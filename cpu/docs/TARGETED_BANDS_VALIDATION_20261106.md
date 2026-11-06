# Native targeted bands: validation

The unified executable supports `Band.NumBands` and `Band.ReferenceEnergy` in
ordinary input. The standalone dependency environment is `cpu/environment.yml`;
no moirekp code is imported. Native NSOC and SOC, scalar/complex overlap, sparse
Bloch phases and native `.Band` formatting are covered by the implementation.

## Correctness

On node040, 4 MPI x 1, a 3-atom HgI2 fixture completed independent NSOC and own
one-step SOC with three native path points. The same candidate binary ran both
the unchanged native full-spectrum path and the 20-state path:

| Case | Full dimension | Selected | Maximum energy difference, Ha |
|---|---:|---:|---:|
| NSOC | 57 | 20 | 2.00e-15 |
| SOC | 114 | 20 | 4.00e-15 |

Additional symmetry-enabled runs completed NSOC and projected SOC, including
HS-file output both enabled and disabled and an absolute eV reference. The
projected SOC target stage runs after forces/grid output. Its later H/S export
verified the projected tensor hash and exactly one export handoff. Compared with
the earlier same-projector fixture, energy and forces were unchanged at printed
precision, H differed by at most 5.38e-14, and S by 3.34e-16.

MoTe2 AA, **9.43 degrees, 222 atoms**, also completed one projected SOC step from
its matching converged 64-rank NSOC checkpoint and two targeted native path
points. The matrix dimension is 7696; both points returned 20 states. Maximum
generalized residual was 4.53e-11, S-orthogonality defect 1.51e-11. The independent
solved-versus-exported operator check was below 2.8e-15; final H/S export verified
the identical projected-operator hash. This is one-step SOC, not SOC convergence.

## Timing probe

The initial integration candidate measured the same native NSOC band path on
node040, **64 MPI x 1**, MoTe2 9.43 degrees, 222 atoms, 3848 AOs, Gamma and M.
Both branches used three fixed NSOC steps, the same inputs and candidate binary,
PBE, 150 Ry, 300 K, History30 and StartPulay15. This is a solver probe, not a
converged production workflow. Each branch was measured once.

| Path | Native band-stage wall, s | Whole fixed-step probe wall, s |
|---|---:|---:|
| Full native spectrum | 306.81 | 334.52 |
| SLEPc, nearest 20 | 55.08 | 82.45 |

Observed band-stage speedup: **5.57x**. Maximum selected-energy difference:
4.03e-15 Ha. This pair belongs to the initial candidate binary
`f30c6510c63b4153c8e24f084958eb3ede809443b46151e6aa2a2ecfa1fd856f`, before
the subsequent projected-SOC deferral/compatibility fix. It is not a repeated
performance measurement of the final release, a memory benchmark, or a SOC
speedup claim. Final SOC compatibility and numerical checks are separate above.

## Installation and scope

All **311 translation units** were compiled from an empty output directory with
Intel oneAPI 2022.3 classic, Intel MPI 2021.7.1, MKL 2022.2.1, CORE-AVX2. There were
no external OpenMX/ELPA objects. The test dependency environment was the existing
`moirekp` environment (Python3.11, complex PETSc3.25.2/SLEPc3.25.1, MUMPS); only its
standard numerical dependencies are used. The public environment file provides
those dependencies under the independent name `openmx-cpu`.

Full CPU regressions pass. The fresh source build also passes a projected SOC
installation smoke with 20 states, absolute eV reference and H/S output disabled.
Multi-node targeted-band performance, contracted/LNO bases and overlap bases
requiring stock regularization are not claimed certified by these tests.
