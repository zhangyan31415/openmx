# Targeted native band dispersion

Approved interface: keep ordinary OpenMX launch and `.Band` format. A positive
`Band.NumBands` selects SLEPc; `Band.ReferenceEnergy` is `Fermi` (default) or an
absolute value in eV. Zero/omitted count keeps the upstream full-spectrum path.

The native path snapshots the actual distributed neighbor H/S blocks once,
including SOC imaginary components and periodic translations. Rank-local
workers construct sparse Bloch matrices, solve the generalized Hermitian
problem near the target, and return exactly the requested number of energies.
There is no dense S diagonalization, TAPW truncation, S shift or identity S.
Use the native lower-triangle convention, report Hermiticity defects, and check
overlap conditioning, generalized residuals, S orthogonality and returned count.

SLEPc runs through an automatically invoked isolated process so the OpenMX
oneAPI MPI ABI is not mixed with the dependency environment's PETSc MPI ABI.
The build bundles the helper and records its Python; a standalone OpenMX
environment file contains all Python/PETSc/SLEPc dependencies. No moirekp import
or private path is permitted in published code. Backend capability is checked
before SCF. No user-written helper command, JSON, export or conversion is needed.

Implementation: add numerical/input/packaging tests, sparse snapshot backend,
native dispatch and build resources, then exercise NSOC and SOC band paths on
small fixtures. Compare a few large-angle k points against native full spectra
using the same matrices, record stage and total wall time separately, run CPU
regressions and review before publishing. Preserve all archived results.
