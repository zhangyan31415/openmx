"""Repository-local OpenMX s/p/d/f real orbital representations of O(3).

Comp2Real data are frozen from OpenMX 4.0 (provenance in comp2real.json).
Cartesian rotations are active, acting on column positions. Coefficient blocks
are active too: p uses (px,py,pz), so D_p(R)=R and D(RS)=D(R)D(S).
No installed OpenMX source or tapw import is needed at runtime.
"""
from functools import lru_cache
import json
from pathlib import Path
import re

import numpy as np
from scipy.linalg import block_diag, expm
from scipy.spatial.transform import Rotation


@lru_cache(maxsize=1)
def _tables():
    document = json.loads(Path(__file__).with_name('comp2real.json').read_text())
    tables = []
    for l, raw in enumerate(document['tables_real_imag']):
        parts = np.asarray(raw, dtype=float)
        c = parts[..., 0] + 1j * parts[..., 1]
        if c.shape != (2*l+1, 2*l+1) or not np.allclose(c @ c.conj().T, np.eye(2*l+1), atol=2e-15, rtol=0):
            raise ValueError('invalid repository Comp2Real table')
        c.setflags(write=False)
        tables.append(c)
    if len(tables) != 4:
        raise ValueError('expected s/p/d/f Comp2Real tables')
    return tuple(tables)


def angular_tables():
    """Return independent complex copies of the four OpenMX Comp2Real tables."""
    return [c.copy() for c in _tables()]


@lru_cache(maxsize=4)
def _generators(l):
    m = np.arange(-l, l+1, dtype=float)
    jp = np.zeros((2*l+1, 2*l+1), dtype=complex)
    for i, v in enumerate(m[:-1]):
        jp[i+1, i] = np.sqrt(l*(l+1)-v*(v+1))
    jm = jp.conj().T
    return (jp+jm)/2, (jp-jm)/(2j), np.diag(m)


def _orthogonal(rotation):
    raw = np.asarray(rotation)
    if np.iscomplexobj(raw):
        raise ValueError('Cartesian operation must be real')
    r = np.asarray(raw, dtype=float)
    if r.shape != (3, 3) or not np.isfinite(r).all():
        raise ValueError('Cartesian operation must be a finite 3x3 matrix')
    if not np.allclose(r.T @ r, np.eye(3), atol=2e-10, rtol=0):
        raise ValueError('Cartesian operation must be orthogonal')
    return r


def get_any_rot_orb_twostep(shell, rotation):
    """Return OpenMX real-harmonic block, matching the historical tapw API.

    Improper operations use inversion parity times the representation of -R.
    No nearest-rotation projection is performed; nonorthogonal inputs fail.
    """
    shell = str(shell).strip()
    if shell not in ('s', 'p', 'd', 'f'):
        raise ValueError(f'unsupported orbital shell: {shell!r}')
    l = 'spdf'.index(shell)
    r = _orthogonal(rotation)
    if l == 0:
        return np.ones((1, 1))
    parity = 1 if np.linalg.det(r) > 0 else -1
    rotvec = Rotation.from_matrix(parity*r).as_rotvec()
    generators = _generators(l)
    exponent = -1j*sum(a*j for a, j in zip(rotvec, generators))
    c = _tables()[l]
    result = (parity**l) * (c.conj() @ expm(exponent) @ c.T)
    if np.max(np.abs(result.imag)) > 5e-12:
        raise ArithmeticError('real harmonic representation acquired imaginary component')
    return result.real


def orbital_rotation_block(orbital_spec, rotation):
    """Block diagonal rotation for e.g. 's3p2d2' or 'Mo7.0-s3p2d2'.

    Orbitals are ordered by l, radial multiplicity, then Comp2Real row.
    """
    spec = str(orbital_spec).strip().rsplit('-', 1)[-1]
    if not re.fullmatch(r'(?:[spdf]\d+)+', spec):
        raise ValueError(f'invalid orbital basis: {orbital_spec!r}')
    parsed = re.findall(r'([spdf])(\d+)', spec)
    indices = ['spdf'.index(s) for s, _ in parsed]
    if indices != sorted(set(indices)):
        raise ValueError('basis shells must occur once in increasing s/p/d/f order')
    r = _orthogonal(rotation)
    blocks = []
    for shell, count in parsed:
        if int(count):
            d = get_any_rot_orb_twostep(shell, r)
            blocks.extend([d]*int(count))
    if not blocks:
        raise ValueError('basis must contain at least one orbital')
    return block_diag(*blocks)
