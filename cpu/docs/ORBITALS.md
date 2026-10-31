# Repository-local OpenMX real orbital rotations

`cpu/geometry/openmx_orbitals.py` replaces the runtime dependence on both an
installed `tapw` package and an external OpenMX `SetPara_DFT.c` source tree.
It requires NumPy and SciPy. No external source path is read at runtime.

## Public API

```python
from cpu.geometry.openmx_orbitals import (
    angular_tables,
    get_any_rot_orb_twostep,
    orbital_rotation_block,
)

c_s, c_p, c_d, c_f = angular_tables()
p_rotation = get_any_rot_orb_twostep('p', R)
ao_rotation = orbital_rotation_block('Mo7.0-s3p2d2', R)
```

`angular_tables()` returns independent complex copies, so caller mutations cannot
corrupt a later call. Complex spherical-harmonic columns are ordered from
`m=-l` to `m=+l`. Real rows follow OpenMX's convention:

| Shell | Real row ordering |
| --- | --- |
| s | s |
| p | px, py, pz |
| d | dz², dx²−y², dxy, dxz, dyz |
| f | fz³, fxz², fyz², fz(x²−y²), fxyz, fx(x²−3y²), fy(3x²−y²) |

The f labels are the standard solid harmonics, not literal monomials. Their
normalized polynomial definitions are independently checked by the unit test.
`orbital_rotation_block` orders shells by l, then radial multiplicity, then
these real rows. It accepts either a full OpenMX basis designation or its
`s3p2d2` suffix; duplicated, out-of-order or unsupported shells fail explicitly.

`R` acts actively on Cartesian column positions. The orbital representation is
active too: `D_p(R)=R` and `D(R1 @ R2)=D(R1) @ D(R2)`. Equivalently,
`f_j(R^-1 x) = sum_i f_i(x) D_ij(R)`. Proper rotations are built with integer-l
angular-momentum generators in the complex basis and transformed using
`C.conj() @ exp(-i theta axis.J) @ C.T`. For improper operations the result is
`(-1)^l D_l(-R)`. Inputs must be finite real orthogonal 3×3 matrices within
`2e-10` absolute orthogonality error; invalid geometry operations fail rather
than being intentionally projected to a new operation. SO(3) conversion can
remove floating-point roundoff within that accepted tolerance.

This module constructs the orbital representation only. It does not certify
that a structure, Hamiltonian, overlap, spin state or numerical grid is invariant
under a proposed operation. It supplies no spinor SU(2) factor or automatic
approximation/projection of H/S.

## Frozen data provenance

`comp2real.json` contains real/imaginary pairs extracted once from OpenMX 4.0's
`source/SetPara_DFT.c`. The original file was
`snapshots/home/software/openmx4.0/source/SetPara_DFT.c`, SHA256:

```
e594f08b4632b19dea5112f3269cd4d4a632ba9ce5344b2dbc521dfa76dac1ad
```

Only literal `Comp2Real[l][row][column] = Complex(real, imag)` assignments for
l=0,1,2,3 were extracted. Their numeric `sqrt(2)` expressions were evaluated
once when generating the table. Runtime loads the JSON, checks dimensions and
unitarity, and never evaluates source text.

## Bounded verification, 2026-10-31

Test host: `node042`, direct SSH. Environment prefix:
`snapshots/home/mambaforge`, Python 3.12.7. OMP/MKL/OpenBLAS threads all 1.
No MPI or production SCF jobs were started by these tests.

```bash
cd snapshots/work/software/2.mp_openmx/openmx-dev
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
snapshots/home/mambaforge/bin/python -m unittest discover \
  -s cpu/tests -p test_orbitals.py -v
```

Six tests pass. They cover Comp2Real ordering/unitarity and copy isolation,
Cartesian p transformations, proper/improper composition, inversion parity,
z-axis rotations, arbitrary C2 axes, explicit s/p/d/f solid-harmonic polynomial
substitution, multi-radial block ordering and invalid input rejection.

An additional one-time reference comparison loaded the established
`snapshots/work/software/1.tapw_code/moirekp-release/tapw/tapw/geometry/rotations.py`
directly through Python's file loader. Four deterministic random O(3) operations
(two proper, two improper), all four shells, gave maximum absolute difference
`9.992007221626409e-16`. This reference is not required for ordinary tests.

Tests were run synchronously and exited successfully; no retained compute
processes. Before implementation, all five initial tests failed because the
repository-local provider did not exist. The z-rotation analytic test was then
corrected to use OpenMX's d-shell ordering (m=2 pair before m=1), independently
confirmed by both the explicit polynomials and the historical provider.
