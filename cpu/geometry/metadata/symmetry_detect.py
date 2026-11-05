#!/usr/bin/env python3
"""Generic crystallographic symmetry detector for OpenMX inputs.

This is a geometry/transport metadata generator.  It deliberately does not
modify an OpenMX input, project H/S, or enable an SCF symmetry path.  The
output records complete space-group operations (W,t), species-preserving atom
permutations, periodic image shifts, residuals, a tolerance persistence scan,
and finite-group closure diagnostics.

The first implementation targets fully periodic three-dimensional OpenMX
cells (pbc=111).  Slabs and non-periodic clusters need an explicit layer-group
or molecular detector rather than treating vacuum as a third periodic axis.

The command-line control policy is explicit:
  --mode auto   scan tolerances and select the largest stable operation set
  --mode user   validate one requested tolerance and optional point group;
                 a mismatch is an error and never silently becomes P1
                 --operations-file can replace spglib's full group with an
                 explicit closed subgroup
  --mode off    emit raw metadata with no symmetry operations
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import re
import sys
from pathlib import Path

import numpy as np
try:
    import spglib
except ImportError:  # mode=off can emit raw metadata without it
    spglib = None
from scipy.linalg import block_diag
from scipy.spatial import cKDTree

BOHR_TO_ANG = 0.529177210903


def _strip(line: str) -> str:
    return line.split("#", 1)[0].strip()


def _first_value(lines: list[str], key: str) -> str | None:
    key = key.lower()
    for line in lines:
        parts = _strip(line).split()
        if parts and parts[0].lower() == key and len(parts) >= 2:
            return parts[1]
    return None


def _section(lines: list[str], name: str) -> list[str]:
    start = "<" + name.lower()
    end = name.lower() + ">"
    active = False
    out: list[str] = []
    for line in lines:
        text = _strip(line)
        low = text.lower()
        if not active and low.startswith(start):
            active = True
            continue
        if active:
            if low.startswith(end):
                return out
            if text:
                out.append(text)
    raise ValueError(f"OpenMX section <{name} ... {name}> is missing")


def _unit_scale(value: str | None) -> float:
    value = (value or "ang").lower()
    if value.startswith("ang"):
        return 1.0
    if value.startswith("bohr") or value.startswith("au"):
        return BOHR_TO_ANG
    if value.startswith("frac"):
        return 1.0
    raise ValueError(f"Unsupported OpenMX coordinate unit {value!r}")


def parse_openmx(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8").splitlines()
    pbc = "111"
    if len(pbc) != 3 or any(ch not in "01" for ch in pbc):
        raise ValueError("Only pbc=111 is supported by this first detector")

    lattice_unit = _first_value(lines, "Atoms.UnitVectors.Unit")
    lattice_scale = _unit_scale(lattice_unit)
    lattice_rows = []
    for row in _section(lines, "Atoms.UnitVectors"):
        values = row.split()
        if len(values) >= 3:
            lattice_rows.append([float(values[0]), float(values[1]), float(values[2])])
    if len(lattice_rows) != 3:
        raise ValueError("Expected exactly three Atoms.UnitVectors rows")
    cell = np.asarray(lattice_rows, dtype=float) * lattice_scale

    species_basis: dict[str, str] = {}
    for row in _section(lines, "Definition.of.Atomic.Species"):
        values = row.split()
        if len(values) >= 2:
            species_basis[values[0]] = values[1]
    if not species_basis:
        raise ValueError("No Definition.of.Atomic.Species entries found")

    coord_unit = _first_value(lines, "Atoms.SpeciesAndCoordinates.Unit") or "Ang"
    coord_is_fractional = coord_unit.lower().startswith("frac")
    coord_scale = _unit_scale(coord_unit)
    inv_cell = np.linalg.inv(cell)
    atoms = []
    for row in _section(lines, "Atoms.SpeciesAndCoordinates"):
        values = row.split()
        if len(values) < 5:
            continue
        index = int(values[0])
        species = values[1]
        if species not in species_basis:
            raise ValueError(f"Atom {index} uses undefined species {species!r}")
        raw = np.asarray([float(v) for v in values[2:5]], dtype=float)
        if coord_is_fractional:
            frac = raw
            cart = frac @ cell
        else:
            cart = raw * coord_scale
            frac = cart @ inv_cell
        atoms.append({"index": index, "species": species, "basis": species_basis[species],
                      "frac": frac, "cart": cart})
    if not atoms:
        raise ValueError("No atoms found")
    expected = _first_value(lines, "Atoms.Number")
    if expected is not None and int(expected) != len(atoms):
        raise ValueError(f"Atoms.Number={expected} but parsed {len(atoms)} atoms")
    frac = np.vstack([a["frac"] for a in atoms])
    species = [a["species"] for a in atoms]
    labels = {name: i + 1 for i, name in enumerate(sorted(set(species)))}
    numbers = np.asarray([labels[name] for name in species], dtype=int)
    spin = (_first_value(lines, "scf.SpinPolarization") or "off").lower()
    soc = (_first_value(lines, "scf.SpinOrbit.Coupling") or "off").lower()
    if spin == "nc" or soc == "on":
        spin_mode = "spinor"
    elif spin == "on":
        spin_mode = "collinear"
    else:
        spin_mode = "none"
    return {
        "cell": cell,
        "frac": frac,
        "cart": frac @ cell,
        "species": species,
        "species_numbers": numbers,
        "species_basis": species_basis,
        "natoms": len(atoms),
        "pbc": pbc,
        "spin_mode": spin_mode,
    }


def _operation_key(rotation: np.ndarray, translation: np.ndarray) -> tuple:
    t = np.mod(np.asarray(translation, dtype=float), 1.0)
    t[np.abs(t - 1.0) < 1e-10] = 0.0
    return tuple(np.asarray(rotation, dtype=int).reshape(-1).tolist()) + tuple(np.round(t, 8).tolist())


def _species_groups(structure: dict) -> dict:
    groups = {}
    cell = structure["cell"]
    frac = structure["frac"]
    species = structure["species"]
    for name in sorted(set(species)):
        ids = np.asarray([i for i, s in enumerate(species) if s == name], dtype=int)
        # Keep 27 periodic images.  The stored shift is target-image minus
        # target-cell, so frac_target + shift is the transformed source.
        shifts = np.asarray([[i, j, k] for i in (-1, 0, 1)
                              for j in (-1, 0, 1) for k in (-1, 0, 1)], dtype=int)
        target_cell = np.floor(frac[ids]).astype(np.int64)
        wrapped = frac[ids] - target_cell
        image_frac = np.vstack([wrapped + shift for shift in shifts])
        image_cart = image_frac @ cell
        image_target = np.concatenate([ids for _ in shifts])
        image_shift = np.vstack([np.repeat(shift[None, :], len(ids), axis=0)
                                 for shift in shifts])
        groups[name] = {
            "ids": ids,
            "tree": cKDTree(image_cart),
            "image_target": image_target,
            "image_shift": image_shift,
            "image_count": image_cart.shape[0],
            "target_cell": {int(atom): target_cell[k] for k, atom in enumerate(ids)},
        }
    return groups


def _map_operation(structure: dict, rotation: np.ndarray, translation: np.ndarray,
                   map_tol: float) -> dict:
    frac = structure["frac"]
    cell = structure["cell"]
    species = structure["species"]
    groups = _species_groups(structure)
    natoms = structure["natoms"]
    perm = np.full(natoms, -1, dtype=np.int64)
    shifts = np.zeros((natoms, 3), dtype=np.int64)
    residuals = np.zeros(natoms, dtype=float)
    margins = np.full(natoms, np.inf, dtype=float)
    used: set[int] = set()
    transformed = frac @ rotation.T + translation
    transformed_cell = np.floor(transformed).astype(np.int64)
    query_frac = transformed - transformed_cell
    k_query = min(128, max(16, max(v["image_count"] for v in groups.values())))
    for source in range(natoms):
        group = groups[species[source]]
        distances, matches = group["tree"].query(
            query_frac[source] @ cell, k=min(k_query, group["image_count"])
        )
        distances = np.atleast_1d(distances)
        matches = np.atleast_1d(matches)
        best_by_target: dict[int, tuple[float, int]] = {}
        for distance, match in zip(distances, matches):
            match = int(match)
            target = int(group["image_target"][match])
            candidate = (float(distance), match)
            old = best_by_target.get(target)
            if old is None or candidate[0] < old[0]:
                best_by_target[target] = candidate
        ranked = sorted(best_by_target.items(), key=lambda item: item[1][0])
        if not ranked:
            raise ValueError(f"No candidate target for atom {source}")
        nearest_distance = ranked[0][1][0]
        if len(ranked) > 1:
            margins[source] = ranked[1][1][0] - ranked[0][1][0]
        picked = None
        for target, (distance, match) in ranked:
            if target not in used and distance <= map_tol:
                picked = target, distance, match
                break
        if picked is None:
            raise ValueError(
                f"Atom map failed for source {source}: nearest={nearest_distance:.6e} "
                f"Ang, map_tol={map_tol:.6e} Ang"
            )
        target, distance, match = picked
        perm[source] = target
        shifts[source] = (group["image_shift"][match] + transformed_cell[source]
                          - group["target_cell"][int(target)])
        residuals[source] = distance
        used.add(target)
    if len(used) != natoms or not np.array_equal(np.sort(perm), np.arange(natoms)):
        raise ValueError("Atom map is not a bijection")
    return {
        "perm": perm,
        "image_shift": shifts,
        "residual_max_Ang": float(np.max(residuals)),
        "residual_rms_Ang": float(np.sqrt(np.mean(residuals**2))),
        "mapping_margin_min_Ang": float(np.min(margins)),
        "mapping_margin_median_Ang": float(np.median(margins)),
    }


def _pointgroup(rotations: np.ndarray) -> dict:
    try:
        result = spglib.get_pointgroup(rotations)
        # spglib releases have returned both tuple and list containers.
        if isinstance(result, (tuple, list)):
            symbol = str(result[0])
            number = int(result[1]) if len(result) > 1 else None
            transform = result[2].tolist() if len(result) > 2 and hasattr(result[2], "tolist") else (
                result[2] if len(result) > 2 else None
            )
            return {"symbol": symbol, "number": number, "transformation": transform}
        return {"symbol": str(result), "number": None, "transformation": None}
    except Exception as exc:
        return {"symbol": None, "number": None, "transformation": None, "error": str(exc)}


_ROTATION_PROVIDER = None
_ROTATION_PROVIDER_ERROR = None


def _load_orbital_rotation_provider():
    """Load the established OpenMX/tapw real-harmonic s/p/d/f rotations lazily.

    Geometry detection remains usable without this optional package.  In that
    case the metadata says that AO validation was unavailable instead of
    silently treating orbital representations as validated.
    """
    global _ROTATION_PROVIDER, _ROTATION_PROVIDER_ERROR
    if _ROTATION_PROVIDER is not None or _ROTATION_PROVIDER_ERROR is not None:
        return _ROTATION_PROVIDER
    provider_path = Path(__file__).resolve().parents[1] / "openmx_orbitals.py"
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("openmx_local_orbitals", provider_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _ROTATION_PROVIDER = module.get_any_rot_orb_twostep
    except (OSError, ImportError, ValueError) as exc:
        _ROTATION_PROVIDER_ERROR = "repository-local OpenMX orbital provider: " + str(exc)
    return _ROTATION_PROVIDER


def _basis_shells(orbital_spec: str) -> list[tuple[str, int]]:
    text = str(orbital_spec).split("-", 1)[-1]
    shells = [(shell, int(count)) for shell, count in re.findall(r"([spdf])(\d+)", text)]
    if not shells:
        raise ValueError(f"Could not parse orbital basis {orbital_spec!r}")
    return shells


def _nearest_orthogonal(rotation_cart: np.ndarray) -> tuple[np.ndarray, float]:
    """Return the nearest orthogonal matrix and its correction norm.

    A relaxed cell can make the Cartesian form derived from an integer
    fractional operation differ from an exact orthogonal matrix by a small
    amount.  Orbital rotations must use the orthogonal polar factor; the
    correction is retained as a diagnostic rather than hidden.
    """
    matrix = np.asarray(rotation_cart, dtype=float)
    u, _, vt = np.linalg.svd(matrix)
    orthogonal = u @ vt
    return orthogonal, float(np.max(np.abs(matrix - orthogonal)))


def _orbital_rotation_block(orbital_spec: str, rotation_cart: np.ndarray) -> np.ndarray:
    provider = _load_orbital_rotation_provider()
    if provider is None:
        raise RuntimeError(_ROTATION_PROVIDER_ERROR or "orbital rotation provider unavailable")
    blocks = []
    for shell, count in _basis_shells(orbital_spec):
        for _ in range(count):
            blocks.append(np.asarray(provider(shell, rotation_cart), dtype=np.complex128))
    return np.asarray(block_diag(*blocks), dtype=np.complex128)


def _check_orbital_representation(structure: dict, operations: list[dict]) -> dict:
    if structure.get("spin_mode") != "none":
        return {"checked": False, "reason": f"spin mode {structure.get('spin_mode')} requires spin/double-group representation"}
    try:
        by_species = {}
        max_orthogonalization = 0.0
        for species, basis in structure["species_basis"].items():
            by_species[species] = []
            for op in operations:
                orthogonal = op["R_cart_orthogonal"]
                correction = op["R_cart_orthogonalization_correction"]
                max_orthogonalization = max(max_orthogonalization, correction)
                by_species[species].append(_orbital_rotation_block(basis, orthogonal))
        max_unitary = 0.0
        for species, blocks in by_species.items():
            for block in blocks:
                eye = np.eye(block.shape[0], dtype=complex)
                max_unitary = max(max_unitary, float(np.max(np.abs(block.conj().T @ block - eye))))
        keys = {_operation_key(op["W"], op["t"]): i for i, op in enumerate(operations)}
        max_composition = 0.0
        failures = []
        for i, g in enumerate(operations):
            for j, h in enumerate(operations):
                k = keys.get(_operation_key(g["W"] @ h["W"], h["t"] @ g["W"].T + g["t"]))
                if k is None:
                    failures.append({"g": i, "h": j, "reason": "missing_composite"})
                    continue
                for species, blocks in by_species.items():
                    error = np.max(np.abs(blocks[i] @ blocks[j] - blocks[k]))
                    max_composition = max(max_composition, float(error))
                    if error > 1e-10 and len(failures) < 20:
                        failures.append({"g": i, "h": j, "species": species, "error": float(error)})
        return {
            "checked": not failures,
            "provider": "cpu.geometry.openmx_orbitals.get_any_rot_orb_twostep",
            "spin_mode": structure.get("spin_mode"),
            "max_unitarity_error": max_unitary,
            "max_composition_error": max_composition,
            "max_cartesian_orthogonalization_correction": max_orthogonalization,
            "failures": failures[:20],
            "species_orbital_dimensions": {s: int(b[0].shape[0]) for s, b in by_species.items()},
        }
    except Exception as exc:
        return {"checked": False, "reason": repr(exc), "provider_error": _ROTATION_PROVIDER_ERROR}


def _closure(operations: list[dict]) -> dict:
    keys = {_operation_key(op["W"], op["t"]): i for i, op in enumerate(operations)}
    failures = []
    max_shift_variation = 0
    max_perm_mismatch = 0
    for i, g in enumerate(operations):
        for j, h in enumerate(operations):
            # Apply h first, then g.  For row fractional coordinates:
            # s -> s W_h^T+t_h -> s(W_g W_h)^T+t_h W_g^T+t_g.
            wc = g["W"] @ h["W"]
            tc = h["t"] @ g["W"].T + g["t"]
            k = keys.get(_operation_key(wc, tc))
            if k is None:
                failures.append({"g": i, "h": j, "reason": "missing_composite"})
                continue
            expected_perm = g["perm"][h["perm"]]
            if not np.array_equal(expected_perm, operations[k]["perm"]):
                max_perm_mismatch = max(max_perm_mismatch, 1)
                failures.append({"g": i, "h": j, "reason": "permutation_mismatch", "k": k})
            expected_shift = h["image_shift"] @ g["W"].T + g["image_shift"][h["perm"]]
            difference = expected_shift - operations[k]["image_shift"]
            variation = int(np.max(np.abs(difference - difference[0])))
            max_shift_variation = max(max_shift_variation, variation)
            if variation:
                failures.append({"g": i, "h": j, "reason": "image_shift_mismatch", "k": k,
                                 "variation": variation})
    identity = np.eye(3, dtype=int)
    identity_key = _operation_key(identity, np.zeros(3))
    identity_index = keys.get(identity_key)
    identity_ok = identity_index is not None
    if identity_ok:
        identity_ok = np.array_equal(operations[identity_index]["perm"], np.arange(len(operations[identity_index]["perm"])))
    return {
        "closed": not failures and identity_ok,
        "failure_count": len(failures),
        "failures": failures[:20],
        "max_common_shift_variation": max_shift_variation,
        "max_permutation_mismatch": max_perm_mismatch,
        "identity_present": identity_ok,
    }


def _build_operation_set(structure: dict, rotations: list | np.ndarray,
                        translations: list | np.ndarray, symprec: float,
                        map_factor: float, spacegroup: str | None = None,
                        source: str = "explicit") -> dict:
    """Validate and serialize an explicit space-group operation set."""
    rotations = np.asarray(rotations)
    translations = np.asarray(translations, dtype=float)
    if rotations.ndim != 3 or rotations.shape[1:] != (3, 3):
        return {"symprec_Ang": float(symprec), "valid": False,
                "error": "explicit rotations must have shape (n,3,3)"}
    if translations.shape != (rotations.shape[0], 3):
        return {"symprec_Ang": float(symprec), "valid": False,
                "error": "explicit translations must have shape (n,3)"}
    if not np.allclose(rotations, np.rint(rotations), atol=1e-8):
        return {"symprec_Ang": float(symprec), "valid": False,
                "error": "fractional rotations W must be integer matrices"}
    rotations = np.rint(rotations).astype(int)
    determinants = [int(round(np.linalg.det(rotation))) for rotation in rotations]
    if any(abs(det) != 1 for det in determinants):
        return {"symprec_Ang": float(symprec), "valid": False,
                "error": "every W must be unimodular with determinant +/-1"}
    map_tol = float(symprec) * float(map_factor)
    operations = []
    seen = set()
    try:
        for rotation, translation in zip(rotations, translations):
            translation = np.mod(translation, 1.0)
            key = _operation_key(rotation, translation)
            if key in seen:
                return {"symprec_Ang": float(symprec), "valid": False,
                        "error": f"duplicate explicit operation {key}"}
            seen.add(key)
            mapped = _map_operation(structure, rotation, translation, map_tol)
            q = structure["cell"].T @ rotation @ np.linalg.inv(structure["cell"].T)
            q_orthogonal, q_correction = _nearest_orthogonal(q)
            operations.append({
                "W": rotation,
                "t": translation,
                "R_cart": q,
                "R_cart_orthogonal": q_orthogonal,
                "R_cart_orthogonalization_correction": q_correction,
                **mapped,
            })
    except Exception as exc:
        return {"symprec_Ang": float(symprec), "valid": False, "error": str(exc)}
    identity_key = _operation_key(np.eye(3, dtype=int), np.zeros(3))
    operations.sort(key=lambda op: (0 if _operation_key(op["W"], op["t"]) == identity_key else 1,
                                   _operation_key(op["W"], op["t"])))
    closure = _closure(operations)
    if closure["closed"]:
        # One common metric preserves all group products. Independently taking
        # polar factors can destroy closure for a nearly symmetric cell.
        metric = sum(op["R_cart"].T @ op["R_cart"] for op in operations) / len(operations)
        values, vectors = np.linalg.eigh((metric + metric.T) * .5)
        if np.min(values) <= 0:
            raise ValueError("invalid invariant Cartesian metric")
        root = (vectors * np.sqrt(values)) @ vectors.T
        inverse_root = (vectors / np.sqrt(values)) @ vectors.T
        for op in operations:
            orthogonal = root @ op["R_cart"] @ inverse_root
            op["R_cart_orthogonal"] = orthogonal
            op["R_cart_orthogonalization_correction"] = float(np.max(abs(orthogonal-op["R_cart"])))
    orbital_representation = _check_orbital_representation(structure, operations)
    pg = _pointgroup(np.asarray([op["W"] for op in operations], dtype=int))
    max_residual = max(op["residual_max_Ang"] for op in operations)
    min_margin = min(op["mapping_margin_min_Ang"] for op in operations)
    return {
        "symprec_Ang": float(symprec),
        "map_tol_Ang": map_tol,
        "valid": bool(closure["closed"]),
        "operation_count": len(operations),
        "operation_signature": [
            _operation_key(op["W"], op["t"]) for op in operations
        ],
        "space_group": spacegroup,
        "point_group": pg,
        "operation_source": source,
        "cartesian_representation_method": "common_group_invariant_metric; coordinates unchanged",
        "orbital_representation": orbital_representation,
        "max_geometry_residual_Ang": float(max_residual),
        "min_mapping_margin_Ang": float(min_margin),
        "operations": [{
            "W": op["W"].tolist(),
            "t": op["t"].tolist(),
            "R_cart": op["R_cart"].tolist(),
            "R_cart_orthogonal": op["R_cart_orthogonal"].tolist(),
            "R_cart_orthogonalization_correction": op["R_cart_orthogonalization_correction"],
            "atom_permutation": op["perm"].tolist(),
            "image_shift": op["image_shift"].tolist(),
            "residual_max_Ang": op["residual_max_Ang"],
            "residual_rms_Ang": op["residual_rms_Ang"],
            "mapping_margin_min_Ang": op["mapping_margin_min_Ang"],
            "mapping_margin_median_Ang": op["mapping_margin_median_Ang"],
        } for op in operations],
        "closure": closure,
    }


def _read_operations_file(path: Path) -> tuple[list, list]:
    """Read a JSON operation file, accepting W/t or rotation/translation keys."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data, dict):
        entries = data.get("operations")
        if entries is None and isinstance(data.get("selected"), dict):
            entries = data["selected"].get("operations")
    else:
        entries = data
    if not isinstance(entries, list) or not entries:
        raise ValueError("operations file must contain a non-empty operations list")
    rotations = []
    translations = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"operation {index} is not an object")
        rotation = entry.get("W", entry.get("rotation"))
        translation = entry.get("t", entry.get("translation", [0.0, 0.0, 0.0]))
        if rotation is None:
            raise ValueError(f"operation {index} is missing W/rotation")
        rotations.append(rotation)
        translations.append(translation)
    return rotations, translations


def detect_at_tolerance(structure: dict, symprec: float, map_factor: float) -> dict:
    cell = structure["cell"]
    frac = structure["frac"]
    numbers = structure["species_numbers"]
    dataset = spglib.get_symmetry((cell, frac, numbers), symprec=float(symprec))
    if dataset is None:
        return {"symprec_Ang": symprec, "valid": False, "error": "spglib returned None"}
    rotations = np.asarray(dataset["rotations"], dtype=int)
    translations = np.asarray(dataset["translations"], dtype=float)
    try:
        spacegroup = str(spglib.get_spacegroup((cell, frac, numbers), symprec=float(symprec)))
    except Exception:
        spacegroup = None
    return _build_operation_set(
        structure, rotations, translations, symprec, map_factor,
        spacegroup=spacegroup, source="spglib"
    )


def _stable_pair(scan: list[dict]) -> tuple[int, int] | None:
    # There is often a low-symmetry P1 plateau at tight tolerances followed
    # by a higher-symmetry plateau.  Select the stable pair with the largest
    # operation set, then the smallest tolerance on that plateau.  This still
    # records the full ladder so users can see whether the selection is near.
    candidates = []
    for i in range(len(scan) - 1):
        left, right = scan[i], scan[i + 1]
        if left.get("valid") and right.get("valid") and left.get("operation_signature") == right.get("operation_signature"):
            candidates.append((int(left.get("operation_count", 0)), -i, i, i + 1))
    if not candidates:
        return None
    _, _, i, j = max(candidates)
    return i, j


def _normalise_point_group(value: str | None) -> str | None:
    """Map common Schoenflies aliases to spglib symbols."""
    if value is None:
        return None
    text = str(value).strip().lower().replace(" ", "").replace("_", "")
    aliases = {
        "c1": "1", "c2": "2", "c3": "3", "c4": "4", "c6": "6",
        "d2": "222", "d3": "32", "d4": "422", "d6": "622",
        "cs": "m", "ci": "-1", "c2v": "mm2", "c2h": "2/m",
        "c3v": "3m", "c4v": "4mm", "c6v": "6mm",
        "d2h": "mmm", "d3d": "-3m", "d3h": "-6m2",
        "d4h": "4/mmm", "d6h": "6/mmm", "s6": "-3",
    }
    return aliases.get(text, text)


def _raw_selection(reason: str) -> dict:
    return {
        "valid": False,
        "operation_count": 0,
        "point_group": {"symbol": "1", "number": 1, "transformation": None},
        "operations": [],
        "reason": reason,
    }


def run(input_path: Path, output_path: Path, tolerances: list[float],
        strict_residual: float, map_factor: float, mode: str = "auto",
        requested_point_group: str | None = None,
        requested_symprec: float | None = None,
        operations_file: Path | None = None) -> dict:
    mode = str(mode).lower()
    if mode not in {"auto", "user", "off"}:
        raise ValueError(f"Unknown symmetry mode {mode!r}; use auto, user, or off")
    if operations_file is not None and mode != "user":
        raise ValueError("--operations-file is only valid with --mode user")
    structure = parse_openmx(input_path)
    user_group = None
    user_tol = None
    if mode == "off":
        scan = []
        stable = None
        selected_index = None
        selected = _raw_selection("symmetry disabled by user")
        status = "raw"
    elif mode == "user":
        if requested_symprec is None or requested_symprec <= 0.0:
            raise ValueError("user mode requires --symprec > 0")
        user_tol = float(requested_symprec)
        if operations_file is not None:
            rotations, translations = _read_operations_file(operations_file)
            selected = _build_operation_set(
                structure, rotations, translations, user_tol, map_factor,
                source=f"operations_file:{operations_file.name}"
            )
            scan = [selected]
        else:
            selected = detect_at_tolerance(structure, user_tol, map_factor)
            scan = [selected]
        stable = None
        selected_index = 0
        if not selected.get("valid"):
            reason = selected.get("error")
            if reason is None:
                reason = f"operation-set validation failed: {selected.get('closure')}"
            raise ValueError(
                "requested symmetry is geometrically invalid at "
                f"symprec={user_tol:g} Ang: {reason}"
            )
        user_group = _normalise_point_group(requested_point_group)
        detected_group = _normalise_point_group(
            selected.get("point_group", {}).get("symbol")
        )
        if user_group is not None and user_group != detected_group:
            raise ValueError(
                "requested point group does not match the validated geometry: "
                f"requested={requested_point_group!r} ({user_group}), "
                f"detected={selected.get('point_group', {}).get('symbol')!r} "
                f"({detected_group}), operations={selected.get('operation_count')}"
            )
        status = (
            "user_forced_strict"
            if selected.get("max_geometry_residual_Ang", np.inf) <= strict_residual
            else "user_forced_near"
        )
    else:
        scan = [detect_at_tolerance(structure, tol, map_factor) for tol in tolerances]
        stable = _stable_pair(scan)
        if stable is None:
            status = "raw"
            selected = next((item for item in reversed(scan) if item.get("valid")), None)
            selected_index = None if selected is None else scan.index(selected)
        else:
            selected_index = stable[0]
            selected = scan[selected_index]
            if selected["max_geometry_residual_Ang"] <= strict_residual:
                status = "strict_geometry"
            else:
                status = "near_geometry"
        if selected is None:
            selected = _raw_selection("no valid periodic symmetry operation set")
    if mode != "off" and spglib is None:
        raise RuntimeError(
            "spglib is required for auto/user symmetry detection; "
            "install it or run with --mode off"
        )
    metadata = {
        "schema": "openmx-generic-symmetry-metadata-v1",
        "input": str(input_path.resolve()),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": sys.executable,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "spglib_version": getattr(spglib, "__version__", None),
        "host": platform.node(),
        "pbc": structure["pbc"],
        "natoms": structure["natoms"],
        "cell_Ang": structure["cell"].tolist(),
        "species": structure["species"],
        "species_basis": structure["species_basis"],
        "spin_mode": structure["spin_mode"],
        "symmetry_mode": mode,
        "requested_point_group": requested_point_group,
        "requested_point_group_canonical": user_group,
        "requested_symprec_Ang": user_tol,
        "operations_file": str(operations_file.resolve()) if operations_file is not None else None,
        "tolerance_ladder_Ang": tolerances,
        "strict_residual_Ang": strict_residual,
        "map_factor": map_factor,
        "selected_scan_index": selected_index,
        "stable_pair": list(stable) if stable is not None else None,
        "symmetry_status": status,
        "solver_eligibility": "candidate_only",
        "neighbor_closure_checked": False,
        "orbital_representation_checked": bool(selected.get("orbital_representation", {}).get("checked", False)),
        "force_path_supported": False,
        "selected": selected,
        "scan": scan,
        "calculation_category": "geometry detection and transport metadata; no SCF or H/S projection",
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--symprec-ladder", type=float, nargs="+",
                        default=[1e-6, 3e-6, 1e-5, 2e-5, 3e-5, 5e-5, 1e-4])
    parser.add_argument("--strict-residual", type=float, default=5e-6)
    parser.add_argument("--map-factor", type=float, default=1.25)
    parser.add_argument("--mode", choices=("auto", "user", "off"), default="auto",
                        help="auto scans a tolerance ladder; user enforces one tolerance; off emits raw metadata")
    parser.add_argument("--point-group", default=None,
                        help="requested point group in Schoenflies (e.g. D3) or spglib notation (e.g. 32)")
    parser.add_argument("--symprec", type=float, default=None,
                        help="requested geometric tolerance in Angstrom; required in user mode")
    parser.add_argument("--operations-file", type=Path, default=None,
                        help="JSON file containing an explicit closed subgroup of W/t operations")
    args = parser.parse_args()
    metadata = run(args.input.resolve(), args.output.resolve(),
                    args.symprec_ladder, args.strict_residual, args.map_factor,
                    mode=args.mode, requested_point_group=args.point_group,
                    requested_symprec=args.symprec,
                    operations_file=args.operations_file.resolve() if args.operations_file else None)
    selected = metadata["selected"]
    print(json.dumps({
        "input": metadata["input"],
        "status": metadata["symmetry_status"],
        "selected_symprec_Ang": selected.get("symprec_Ang"),
        "operation_count": selected.get("operation_count"),
        "space_group": selected.get("space_group"),
        "point_group": selected.get("point_group"),
        "max_geometry_residual_Ang": selected.get("max_geometry_residual_Ang"),
        "min_mapping_margin_Ang": selected.get("min_mapping_margin_Ang"),
        "orbital_representation": selected.get("orbital_representation"),
        "stable_pair": metadata["stable_pair"],
        "output": str(args.output.resolve()),
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
