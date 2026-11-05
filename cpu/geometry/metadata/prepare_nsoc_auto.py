#!/usr/bin/env python3
"""Build a generic real NSYM v2 sector basis from detected point-group operations.

The detector supplies a finite closed operation set.  For each atom orbit and
radial/angular shell this script builds the real AO representation, forms a
random central element of the group algebra (a class sum), and diagonalizes
that element.  Its eigenspaces are isotypic sectors, so all symmetry commuting
H/S blocks are block diagonal in the resulting orthogonal basis.  D3 keeps
its historical v1 generator and fast path; this script is the generic path.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import struct
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
MAX_SECTORS = 48
MAGIC = 0x4E534433

def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

detector = _module("openmx_generic_symmetry_detect", HERE / "symmetry_detect.py")

def _key(W, t):
    t = np.mod(np.asarray(t, dtype=float), 1.0)
    t[np.abs(t - 1.0) < 1e-10] = 0.0
    return tuple(np.asarray(W, dtype=int).reshape(-1).tolist()) + tuple(np.round(t, 8).tolist())

def _parse_shells(spec):
    return [(l, int(count)) for l, count in detector._basis_shells(spec)]

def _atom_offsets(structure):
    shells = {s: _parse_shells(b) for s, b in structure["species_basis"].items()}
    counts = [sum((2*"spdf".index(l)+1)*nr for l, nr in shells[s])
              for s in structure["species"]]
    offsets = np.r_[0, np.cumsum(counts)].astype(int)
    return shells, offsets, int(offsets[-1])

def _group_tables(operations):
    keys = {_key(op["W"], op["t"]): i for i, op in enumerate(operations)}
    identity = keys.get(_key(np.eye(3, dtype=int), np.zeros(3)))
    if identity is None:
        raise ValueError("operation set has no identity")
    def multiply(i, j):
        g, h = operations[i], operations[j]
        W = np.asarray(g["W"], dtype=int) @ np.asarray(h["W"], dtype=int)
        t = np.asarray(h["t"], dtype=float) @ np.asarray(g["W"], dtype=float).T + np.asarray(g["t"], dtype=float)
        k = keys.get(_key(W, t))
        if k is None:
            raise ValueError(f"operation product ({i},{j}) is absent")
        return k
    inverse = []
    for i in range(len(operations)):
        found = [j for j in range(len(operations)) if multiply(i, j) == identity]
        if len(found) != 1:
            raise ValueError(f"operation {i} has no unique inverse")
        inverse.append(found[0])
    classes = []
    class_of = [-1] * len(operations)
    for i in range(len(operations)):
        if class_of[i] >= 0:
            continue
        members = set()
        for h in range(len(operations)):
            members.add(multiply(multiply(h, i), inverse[h]))
        c = len(classes)
        for j in members:
            class_of[j] = c
        classes.append(sorted(members))
    if any(c < 0 for c in class_of):
        raise ValueError("incomplete conjugacy-class construction")
    return keys, identity, multiply, classes, class_of

def _orbit(atom, perms):
    out = {atom}
    todo = [atom]
    while todo:
        a = todo.pop()
        for p in perms:
            b = int(p[a])
            if b not in out:
                out.add(b)
                todo.append(b)
    return sorted(out)

def _cluster(values, tol):
    groups = []
    for i, value in enumerate(values):
        if not groups or abs(float(value) - groups[-1][-1]) > tol:
            groups.append([float(value)])
        else:
            groups[-1].append(float(value))
    return groups

def _write_v2(output, natoms, n, sizes, columns):
    first = np.zeros(4, dtype="<i4")
    first[:min(4, len(sizes))] = sizes[:4]
    with output.open("wb") as f:
        f.write(struct.pack("<9i", MAGIC, 2, natoms, n, len(sizes), *first))
        if len(sizes) > 4:
            f.write(np.asarray(sizes[4:], dtype="<i4").tobytes())
        nnz = sum(int(np.count_nonzero(np.abs(v) > 1e-13)) for _, v in columns)
        f.write(struct.pack("<q", nnz))
        ident = np.arange(natoms, dtype="<i4")
        f.write(ident.tobytes())
        f.write(np.zeros(3, dtype="<f8").tobytes())
        f.write(ident.tobytes())
        f.write(np.eye(3, dtype="<f8").tobytes())
        f.write(np.zeros(3, dtype="<f8").tobytes())
        ptr = [0]
        rows = []
        vals = []
        for row, value in columns:
            mask = np.abs(value) > 1e-13
            rows.extend(np.asarray(row[mask], dtype="<i4").tolist())
            vals.extend(np.asarray(value[mask], dtype="<f8").tolist())
            ptr.append(len(rows))
        f.write(np.asarray(ptr, dtype="<i8").tobytes())
        f.write(np.asarray(rows, dtype="<i4").tobytes())
        f.write(np.asarray(vals, dtype="<f8").tobytes())
    return nnz

def build(input_path: Path, symmetry_json: Path, output: Path,
          provenance: Path | None = None, eig_tol: float = 2e-8,
          covariance_tol: float = 2e-8):
    structure = detector.parse_openmx(input_path)
    metadata = json.loads(symmetry_json.read_text())
    selected = metadata.get("selected", metadata)
    operations = selected.get("operations", [])
    if not selected.get("valid", False) or len(operations) < 1:
        raise ValueError("symmetry metadata has no validated operation set")
    if structure["spin_mode"] != "none":
        raise ValueError("generic scalar basis is NSOC only; SOC needs a double-group basis")
    shells, offsets, n = _atom_offsets(structure)
    perms = [np.asarray(op["atom_permutation"], dtype=int) for op in operations]
    for p in perms:
        if sorted(p.tolist()) != list(range(structure["natoms"])):
            raise ValueError("operation atom map is not a permutation")
    _, identity, multiply, classes, class_of = _group_tables(operations)
    # Real orbital matrices, cached by species and operation.
    provider = detector._load_orbital_rotation_provider()
    if provider is None:
        raise RuntimeError("OpenMX real-harmonic orbital rotation provider unavailable")
    rotations = {}
    for species, species_shells in shells.items():
        rotations[species] = []
        for op in operations:
            R = np.asarray(op["R_cart_orthogonal"], dtype=float)
            blocks = {}
            for l, nr in species_shells:
                d = 2*"spdf".index(l) + 1
                block = np.asarray(provider(l, R), dtype=np.complex128)
                block = np.real_if_close(block, tol=1000)
                if np.iscomplexobj(block) or block.shape != (d, d):
                    raise ValueError(f"non-real orbital rotation for species={species}, l={l}")
                for radial in range(nr):
                    blocks[(l, radial)] = np.asarray(block, dtype=float)
            rotations[species].append(blocks)
    seen = set()
    sector_columns = {}
    local_checks = []
    for atom in range(structure["natoms"]):
        if atom in seen:
            continue
        orbit = _orbit(atom, perms)
        seen.update(orbit)
        species = structure["species"][atom]
        local_base = 0
        for l, nr in shells[species]:
            d = 2*"spdf".index(l) + 1
            for radial in range(nr):
                rows = np.asarray([offsets[a] + local_base + radial*d + mu
                                   for a in orbit for mu in range(d)], dtype=int)
                size = len(rows)
                reps = []
                where = {a: i for i, a in enumerate(orbit)}
                for gi, op in enumerate(operations):
                    M = np.zeros((size, size), dtype=float)
                    for source in orbit:
                        target = int(perms[gi][source])
                        i0, j0 = where[target] * d, where[source] * d
                        M[i0:i0+d, j0:j0+d] = rotations[species][gi][(l, radial)]
                    reps.append(M)
                central = np.zeros((size, size), dtype=float)
                # Class sums are central in the group algebra and therefore
                # act as scalars on each irreducible representation.
                for ci, members in enumerate(classes):
                    weight = np.sin(0.731 * (ci + 1)) + 0.37 * (ci + 1)
                    for gi in members:
                        central += weight * reps[gi]
                central = (central + central.T) * 0.5
                comm = max(float(np.max(np.abs(central @ M - M @ central)))
                           for M in reps)
                if comm > covariance_tol:
                    raise ValueError(f"central element does not commute: {comm:g}")
                values, vectors = np.linalg.eigh(central)
                groups = _cluster(values, eig_tol * max(1.0, float(np.max(np.abs(values), initial=0))))
                pos = 0
                labels = []
                for group in groups:
                    count = len(group)
                    label = float(np.mean(group))
                    cols = vectors[:, pos:pos+count]
                    pos += count
                    # Verify the eigenspace is invariant under every group op.
                    q = np.eye(size)[:, :0]
                    for M in reps:
                        transformed = cols.T @ M @ cols
                        residual = np.max(np.abs(M @ cols - cols @ transformed), initial=0)
                        if residual > covariance_tol:
                            raise ValueError(f"sector is not invariant: {residual:g}")
                    key = min(sector_columns, key=lambda x: abs(x-label), default=None)
                    if key is None or abs(key-label) > eig_tol * max(1.0, abs(label)):
                        key = label
                        sector_columns[key] = []
                    for col in range(count):
                        sector_columns[key].append((rows, cols[:, col]))
                    labels.append((label, count))
                local_checks.append({"orbit_size": len(orbit), "shell": l,
                                     "radial": radial, "sector_labels": labels})
            local_base += nr * d
    if len(sector_columns) > MAX_SECTORS:
        raise ValueError(f"generated {len(sector_columns)} sectors, maximum is {MAX_SECTORS}")
    ordered = sorted(sector_columns, key=float)
    columns = [item for key in ordered for item in sector_columns[key]]
    sizes = [len(sector_columns[key]) for key in ordered]
    if sum(sizes) != n:
        raise ValueError(f"sector columns span {sum(sizes)} of N={n}")
    output.parent.mkdir(parents=True, exist_ok=True)
    nnz = _write_v2(output, structure["natoms"], n, sizes, columns)
    info = {
        "schema": "openmx-generic-nsoc-sector-basis-v2",
        "input": str(input_path.resolve()),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "symmetry_metadata": str(symmetry_json.resolve()),
        "symmetry_metadata_sha256": hashlib.sha256(symmetry_json.read_bytes()).hexdigest(),
        "operation_count": len(operations),
        "point_group": selected.get("point_group"),
        "sector_count": len(sizes),
        "sector_sizes": sizes,
        "N": n,
        "nnz": nnz,
        "basis_kind": "real_isotypic_central_element",
        "spin_mode": structure["spin_mode"],
        "validation": {
            "group_closed": bool(selected.get("closure", {}).get("closed", False)),
            "orbital_representation_checked": bool(selected.get("orbital_representation", {}).get("checked", False)),
            "central_element_commuting": True,
            "sector_invariance": True,
            "h_s_covariance_checked": False,
        },
        "local_checks": local_checks,
        "final_process_state": "orthogonal sector basis generated; H/S covariance remains a runtime check",
    }
    if provenance:
        provenance.parent.mkdir(parents=True, exist_ok=True)
        provenance.write_text(json.dumps(info, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(info, indent=2, ensure_ascii=False))
    return info

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--symmetry-json", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--provenance", type=Path)
    p.add_argument("--eig-tol", type=float, default=2e-8)
    p.add_argument("--covariance-tol", type=float, default=2e-8)
    a = p.parse_args()
    build(a.input.resolve(), a.symmetry_json.resolve(), a.output.resolve(),
          a.provenance.resolve() if a.provenance else None,
          eig_tol=a.eig_tol, covariance_tol=a.covariance_tol)

if __name__ == "__main__":
    main()
