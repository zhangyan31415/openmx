#!/usr/bin/env python3
"""Geometry-only periodic neighbor support for a validated symmetry action.

No group names, material names or twist angles enter the algorithm.  This
preparation tool emits the existing patched-CPU runtime padding setting; it
does not activate symmetry, alter radial PAOs, or replace native H/S checks.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import re
import socket
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

# Match source/openmx_common.h and Input_std.c, including legacy AU inputs.
BOHR_TO_ANG = 0.529177249
NATIVE_MAX_PADDING_ANG = 0.01


class NeighborSupportError(ValueError):
    """Preparation cannot certify the requested finite neighbor support."""


def _sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def _integer(value, shape, label):
    try:
        a = np.asarray(value, dtype=float)
        if a.shape != shape or not np.isfinite(a).all() or not np.equal(a, np.rint(a)).all():
            raise ValueError()
        if np.max(np.abs(a), initial=0) > np.iinfo(np.int32).max:
            raise ValueError()
        return a.astype(np.int64)
    except (ValueError, TypeError, OverflowError) as exc:
        raise NeighborSupportError(f'{label} must contain finite native int32 integers with shape {shape}') from exc


def _action_key(w, perm, image):
    # A common integer translation cancels from every pair-image difference.
    return (tuple(w.ravel()), perm.tobytes(), (image-image[0]).tobytes())


def _actions(operations, n, radii, species):
    if not operations or len(operations) > 192:
        raise NeighborSupportError('expected 1..192 validated periodic operations')
    result = []
    eye, ids = np.eye(3, dtype=np.int64), np.arange(n, dtype=np.int64)
    for k, op in enumerate(operations):
        w = _integer(op.get('W'), (3, 3), f'operation {k} W')
        a, b, c, d, e, f, g, h, i = map(int, w.ravel())
        if a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g) not in (-1, 1):
            raise NeighborSupportError(f'operation {k} W is not unimodular')
        perm = _integer(op.get('perm', op.get('atom_permutation')), (n,), f'operation {k} permutation')
        image = _integer(op.get('image_shift'), (n, 3), f'operation {k} image_shift')
        if not np.array_equal(np.sort(perm), ids):
            raise NeighborSupportError(f'operation {k} atom permutation is not bijective')
        if species is not None and not np.array_equal(species[perm], species):
            raise NeighborSupportError(f'operation {k} changes atom species')
        if not np.array_equal(radii[perm], radii):
            raise NeighborSupportError(f'operation {k} changes physical PAO radii')
        result.append((w, perm, image))
    keys = {_action_key(*op) for op in result}
    if len(keys) != len(result):
        raise NeighborSupportError('duplicate periodic pair actions')
    if _action_key(eye, ids, np.zeros((n, 3), dtype=np.int64)) not in keys:
        raise NeighborSupportError('operation set is not closed: missing identity')
    for wg, pg, ig in result:
        for wh, ph, ih in result:
            limit = np.iinfo(np.int64).max
            if (3*int(np.max(np.abs(wg)))*int(np.max(np.abs(wh))) > limit or
                    3*int(np.max(np.abs(wg)))*int(np.max(np.abs(ih)))+int(np.max(np.abs(ig))) > limit):
                raise NeighborSupportError('operation composition exceeds safe integer capacity')
            if _action_key(wg@wh, pg[ph], ih@wg.T+ig[ph]) not in keys:
                raise NeighborSupportError('operation atom/image actions are not closed')
    return result


def analyze_support(cell, positions, radii, operations, *, species=None,
                    max_padding=NATIVE_MAX_PADDING_ANG, margin=1e-6,
                    max_image_points=4_000_000, max_candidate_pairs=5_000_000):
    """Return a fixed-point padding certificate for fully periodic geometry.

    Inputs use Angstrom and row-vector cell conventions. Operations must already
    pass the caller's geometry/AO admission; this function independently checks
    species, permutations and complete periodic pair-action group closure.
    Candidate storage is O(neighbor count), never O(AO count squared).
    """
    started = time.perf_counter()
    try:
        cell, xyz, radii = (np.asarray(v, dtype=float) for v in (cell, positions, radii))
    except (ValueError, TypeError) as exc:
        raise NeighborSupportError('invalid numeric geometry') from exc
    if (cell.shape != (3, 3) or xyz.ndim != 2 or xyz.shape[1:] != (3,) or not len(xyz) or
            radii.shape != (len(xyz),) or not all(np.isfinite(v).all() for v in (cell, xyz, radii)) or
            np.any(radii <= 0)):
        raise NeighborSupportError('finite cell/positions and positive PAO radii are required')
    if (not math.isfinite(max_padding) or not 0 <= max_padding <= NATIVE_MAX_PADDING_ANG or
            not math.isfinite(margin) or margin <= 0 or
            not isinstance(max_image_points, int) or max_image_points < 1 or
            not isinstance(max_candidate_pairs, int) or max_candidate_pairs < 1):
        raise NeighborSupportError('invalid padding or preparation resource limits')
    try:
        inverse = np.linalg.inv(cell)
    except np.linalg.LinAlgError as exc:
        raise NeighborSupportError('singular periodic cell') from exc
    if not np.isfinite(inverse).all():
        raise NeighborSupportError('invalid reciprocal cell')
    labels = None if species is None else np.asarray(species)
    if labels is not None and labels.shape != radii.shape:
        raise NeighborSupportError('atom species count differs from geometry')
    actions = _actions(operations, len(xyz), radii, labels)

    frac = xyz@inverse
    if not np.isfinite(frac).all() or np.max(np.abs(frac)) > np.iinfo(np.int32).max:
        raise NeighborSupportError('coordinates exceed safe periodic image capacity')
    wrap = np.floor(frac).astype(np.int64)
    wrapped = (frac-wrap)@cell
    reach = (2*float(radii.max())+max_padding)*np.linalg.norm(inverse, axis=0)
    if not np.isfinite(reach).all():
        raise NeighborSupportError('periodic image bound is not finite')
    # Wrapped atom differences have fractional components in (-1,1).
    bounds = [max(1, math.ceil(float(q)+32*np.finfo(float).eps*max(1., float(q)))) for q in reach]
    image_count = math.prod(2*b+1 for b in bounds)
    image_points = image_count*len(xyz)
    if image_points > max_image_points:
        raise NeighborSupportError(f'periodic image budget exceeded: {image_points} > {max_image_points}; bounds={bounds}')
    shifts = np.asarray(list(itertools.product(*(range(-b, b+1) for b in bounds))), dtype=np.int64)
    points = (wrapped[None, :, :]+(shifts@cell)[:, None, :]).reshape(-1, 3)
    tree = cKDTree(points)
    aa, bb, rr, base = [], [], [], []
    candidate_count = 0
    # Query/transform in bounded chunks, including all reverse pairs and onsite.
    for begin in range(0, len(xyz), 256):
        end = min(len(xyz), begin+256)
        neighbors = tree.query_ball_point(wrapped[begin:end], radii[begin:end]+float(radii.max())+max_padding+1e-10)
        for a, ids in zip(range(begin, end), neighbors):
            ids = np.asarray(ids, dtype=np.int64)
            b, r = ids % len(xyz), shifts[ids//len(xyz)]
            gap = np.linalg.norm(wrapped[b]+r@cell-wrapped[a], axis=1)-radii[a]-radii[b]
            keep = gap <= max_padding
            b, r, gap = b[keep], r[keep], gap[keep]
            candidate_count += len(b)
            if candidate_count > max_candidate_pairs:
                raise NeighborSupportError(f'neighbor candidate budget exceeded: {candidate_count} > {max_candidate_pairs}')
            aa.append(np.full(len(b), a, dtype=np.int64)); bb.append(b); rr.append(r); base.append(gap)
    a, b, images, original_gap = map(np.concatenate, (aa, bb, rr, base))
    low, high = original_gap.copy(), original_gap.copy()
    for w, perm, image in actions:
        shifted = image+wrap[perm]-wrap@w.T
        for begin in range(0, candidate_count, 65536):
            sl = slice(begin, begin+65536)
            ai, bi = a[sl], b[sl]
            mapped = images[sl]@w.T+shifted[bi]-shifted[ai]
            gap = np.linalg.norm(wrapped[perm[bi]]+mapped@cell-wrapped[perm[ai]], axis=1)-radii[ai]-radii[bi]
            low[sl] = np.minimum(low[sl], gap)
            high[sl] = np.maximum(high[sl], gap)
    padding, history = 0., []
    for iteration in range(candidate_count+1):
        active = low <= padding
        needed = max(0., float(np.max(high[active], initial=0.)))
        history.append({'iteration': iteration, 'padding_Ang': padding, 'required_maximum_Ang': needed})
        if needed <= padding:
            break
        if needed > max_padding:
            raise NeighborSupportError(f'closed support requires padding >= {needed:.12g} Ang, above limit {max_padding:.12g}')
        padding = min(max_padding, needed+margin)
    else:
        raise NeighborSupportError('neighbor support fixed point did not converge')
    closed = bool(np.all((low <= padding) == (high <= padding)))
    if not closed:
        raise NeighborSupportError('neighbor support fixed point is not closed')
    return {'schema': 'openmx-neighbor-support-v1', 'closed': True,
            'padding_Ang': padding, 'max_padding_Ang': max_padding, 'margin_Ang': margin,
            'original_blocks': int(np.count_nonzero(original_gap <= 0)),
            'padded_blocks': int(np.count_nonzero(original_gap <= padding)),
            'group_operations': len(actions), 'natoms': len(xyz),
            'periodic_image_bounds': bounds, 'image_points': image_points,
            'candidate_pairs': candidate_count, 'history': history,
            'settings': {'OPENMX_CPU_NEIGHBOR_PAD_ANG': format(padding, '.17g')},
            'geometry_modified': False, 'physical_PAO_radii_modified': False,
            'native_closure_check_required': True,
            'required_native_feature': 'OPENMX_CPU_NEIGHBOR_PAD_ANG' if padding else None,
            'wall_seconds': time.perf_counter()-started}


def _lines(path):
    return [line.split('#', 1)[0].strip() for line in Path(path).read_text().splitlines()]


def _value(lines, key, default=None):
    rows = [line.split()[1:] for line in lines if line.split() and line.split()[0].lower() == key.lower()]
    if len(rows) > 1 or rows and not rows[0]:
        raise NeighborSupportError(f'missing or duplicate {key}')
    return rows[0][0] if rows else default


def _section(lines, name):
    starts = [i for i, line in enumerate(lines) if line.lower() == '<'+name.lower()]
    ends = [i for i, line in enumerate(lines) if line.lower() == name.lower()+'>']
    if len(starts) != 1 or len(ends) != 1 or ends[0] <= starts[0]:
        raise NeighborSupportError(f'invalid OpenMX section {name}')
    return [line.split() for line in lines[starts[0]+1:ends[0]] if line]


def _scale(unit):
    if unit.lower() in ('au', 'bohr'):
        return BOHR_TO_ANG
    if unit.lower() in ('ang', 'angstrom'):
        return 1.
    raise NeighborSupportError(f'unsupported coordinate unit {unit}')


def load_openmx_geometry(input_path, *, data_root=None, runtime_cwd=None):
    """Read real PAO radii; relative DATA.PATH follows explicit launch cwd."""
    input_path = Path(input_path).resolve(); lines = _lines(input_path)
    cwd = input_path.parent if runtime_cwd is None else Path(runtime_cwd).resolve()
    root = Path(data_root) if data_root is not None else Path(_value(lines, 'DATA.PATH', '../DFT_DATA19'))
    root = (cwd/root).resolve() if not root.is_absolute() else root.resolve()
    try:
        cell = np.asarray(_section(lines, 'Atoms.UnitVectors'), dtype=float)*_scale(_value(lines, 'Atoms.UnitVectors.Unit', 'Ang'))
        definitions = _section(lines, 'Definition.of.Atomic.Species')
        atom_rows = _section(lines, 'Atoms.SpeciesAndCoordinates')
        if len({row[0] for row in definitions}) != len(definitions):
            raise NeighborSupportError('duplicate atomic species definition')
        physical, files = {}, []
        for row in definitions:
            species, basis = row[:2]
            name = basis.rsplit('-', 1)[0]
            path = root/'PAO'/(name+'.pao')
            if not path.is_file():
                raise NeighborSupportError(f'PAO data file not found: {path}')
            value = _value(_lines(path), 'radial.cutoff.pao')
            if value is None or not math.isfinite(float(value)) or float(value) <= 0:
                raise NeighborSupportError(f'positive radial.cutoff.pao missing in {path}')
            physical[species] = float(value)*BOHR_TO_ANG
            files.append({'species': species, 'basis': basis, 'path': str(path),
                          'radial_cutoff_Bohr': float(value), 'sha256': _sha(path)})
        count = int(_value(lines, 'Atoms.Number'))
        if count != len(atom_rows) or [int(row[0]) for row in atom_rows] != list(range(1, count+1)):
            raise NeighborSupportError('Atoms.Number or atom ordering mismatch')
        species = [row[1] for row in atom_rows]
        xyz = np.asarray([row[2:5] for row in atom_rows], dtype=float)
        unit = _value(lines, 'Atoms.SpeciesAndCoordinates.Unit', 'Ang')
        xyz = xyz@cell if unit.lower() in ('frac', 'fractional') else xyz*_scale(unit)
        radii = np.asarray([physical[s] for s in species])
    except (IndexError, KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, NeighborSupportError):
            raise
        raise NeighborSupportError(f'invalid OpenMX geometry or PAO metadata: {exc}') from exc
    digest = hashlib.sha256()
    for array in (cell, xyz, radii):
        digest.update(np.asarray(array, dtype='<f8').tobytes())
    digest.update(json.dumps({'species': species, 'PAOs': [
        {key: item[key] for key in ('species', 'basis', 'sha256')} for item in files
    ]}, sort_keys=True).encode())
    return {'cell_Ang': cell, 'positions_Ang': xyz, 'radii_Ang': radii, 'species': species,
            'pao_files': files, 'runtime_cwd': str(cwd), 'data_root': str(root),
            'geometry_basis_sha256': digest.hexdigest()}


def prepare(input_path, operations_path, output_path=None, *, data_root=None,
            runtime_cwd=None, max_padding=NATIVE_MAX_PADDING_ANG, margin=1e-6):
    started = time.perf_counter(); input_path = Path(input_path).resolve(); operations_path = Path(operations_path).resolve()
    try:
        metadata = json.loads(operations_path.read_text())
    except (OSError, ValueError) as exc:
        raise NeighborSupportError(f'cannot read operation metadata: {exc}') from exc
    if not isinstance(metadata, dict):
        raise NeighborSupportError('operation metadata must be a JSON object')
    if metadata.get('input_sha256') and metadata['input_sha256'] != _sha(input_path):
        raise NeighborSupportError('operation metadata input SHA256 differs from current geometry input')
    selected = metadata.get('selected', metadata)
    if not isinstance(selected, dict):
        raise NeighborSupportError('operation metadata has no selected periodic actions')
    if selected.get('valid') is False:
        raise NeighborSupportError('operation metadata failed geometry validation')
    operations = selected.get('operations')
    if not isinstance(operations, list):
        raise NeighborSupportError('operation metadata has no validated atom/image actions')
    geometry = load_openmx_geometry(input_path, data_root=data_root, runtime_cwd=runtime_cwd)
    report = analyze_support(geometry['cell_Ang'], geometry['positions_Ang'], geometry['radii_Ang'], operations,
                             species=geometry['species'], max_padding=max_padding, margin=margin)
    report.update(input=str(input_path), input_sha256=_sha(input_path),
                  operations=str(operations_path), operations_sha256=_sha(operations_path),
                  pao_files=geometry['pao_files'], data_root=geometry['data_root'], runtime_cwd=geometry['runtime_cwd'],
                  geometry_basis_sha256=geometry['geometry_basis_sha256'],
                  Bohr_to_Ang=BOHR_TO_ANG, hostname=socket.gethostname(), environment_prefix=sys.prefix,
                  python_version=platform.python_version(), generator_sha256=_sha(__file__),
                  category='geometry support preparation; no native SCF or energy-equivalence claim',
                  total_wall_seconds=time.perf_counter()-started)
    if output_path is not None:
        _write_json(output_path, report)
    return report


def _write_json(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2)+'\n'); temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path); parser.add_argument('operations', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--settings-output', type=Path)
    parser.add_argument('--data-root', type=Path); parser.add_argument('--runtime-cwd', type=Path)
    parser.add_argument('--max-padding', type=float, default=NATIVE_MAX_PADDING_ANG)
    parser.add_argument('--margin', type=float, default=1e-6)
    args = parser.parse_args()
    try:
        report = prepare(args.input, args.operations, args.output, data_root=args.data_root,
                         runtime_cwd=args.runtime_cwd, max_padding=args.max_padding, margin=args.margin)
    except (NeighborSupportError, OSError) as exc:
        failure = {'schema': 'openmx-neighbor-support-v1', 'closed': False, 'error': str(exc),
                   'settings': {'OPENMX_CPU_NEIGHBOR_PAD_ANG': '0'}, 'native_launch_allowed': False}
        _write_json(args.output, failure)
        if args.settings_output is not None:
            _write_json(args.settings_output, failure['settings'])
        parser.exit(2, f'neighbor support preparation failed: {exc}\n')
    if args.settings_output is not None:
        _write_json(args.settings_output, report['settings'])
    print(json.dumps({k: report[k] for k in ('closed', 'padding_Ang', 'original_blocks', 'padded_blocks', 'group_operations', 'total_wall_seconds')}, indent=2))


if __name__ == '__main__':
    main()
