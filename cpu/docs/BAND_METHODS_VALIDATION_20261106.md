# Band method selection: full-spectrum comparisons

`Band.Solver` accepts Auto, LAPACK/GVD, EVR, SLEPc and Legacy. Explicit Auto uses
EVR for all bands and SLEPc for a positive target count. Omitting the method
retains the earlier count-dependent behavior. The methods affect band
postprocessing, not SCF occupations or its eigensolver.

## Matched native full-band result

On **node043, 4 MPI x 16 threads**, MoTe2 AA **9.43 degrees, 222 atoms, 3848 scalar
AOs**, the same candidate binary and input settings computed four path points
from Gamma to M. Each method returned **all 3848 eigenvalues** at every point.
Inputs used three fixed NSOC iterations, PBE150Ry300K, History30, StartPulay15,
criterion1e-7 and D3. This is a solver probe, not converged production SCF.

| Method | Complete band-stage wall, s | Whole three-step probe wall, s | Maximum energy difference from Legacy, Ha |
|---|---:|---:|---:|
| Legacy | 129.393 | 184.068 | reference |
| LAPACK/GVD | 5.276 | 59.302 | 4.14e-12 |
| EVR | 4.300 | 58.208 | 6.23e-12 |

EVR's observed complete band stage is **30.1x** faster than Legacy and 18.5%
shorter than GVD. Each native method was measured once; this is not a repeated
whole-workflow production speedup. Legacy retains S eigendecomposition, full
vectors and native matrix transforms. The dense paths use Cholesky/generalized
drivers, optimized BLAS and eigenvalues only, so the change is not solely a
comparison of two final LAPACK routine names. Legacy links oneAPI MKL2022.2;
the Python dense drivers here use MKL2026.0. Kernel-driver comparisons below
use the same loaded MKL for both methods.

Largest worker lifetime peak RSS was 3.4635 GiB for GVD and 3.4653 GiB for EVR.
This includes that worker's input, BLAS and workspaces. It is not a node physical
peak, must not be summed as physical node memory, and no matching Legacy memory
measurement was collected. Each concurrent dense k-point worker owns matrices.

## Kernel size and threading probe

Identical H/S pairs were compared within each node, with driver order reversed
between two measurements, and 1/4/16 threads. Node041 used a 512-order principal
block; node042 a 1536-order principal block; node043 the complete 3848-order
matrix. The smaller blocks are algorithm-size probes, not independent material
calculations. Actual MKL2026.0 loading was recorded. The following are means of
two **16-thread** kernel calls, excluding file loading and runtime certification:

| Order | GVD, s | Cholesky + HEGST + EVR, s |
|---|---:|---:|
| 512 principal block | 0.0229 | 0.0179 |
| 1536 principal block | 0.1802 | 0.1731 |
| 3848 full matrix | 2.0814 | 1.5033 |

All spectra agreed within 7e-14 Ha across methods/thread counts. The complete
3848 spectrum agreed with the saved native spectrum within 1.78e-13 Ha.
Small absolute differences, initialization and the two-repeat limit prevent
treating every small percentage as a reliable win. Kernel RSS in this first
size sweep is cumulative within a process, so it is not used for method-memory
comparisons. For 57/114-order native fixtures, Legacy avoided worker startup and
was much shorter than external dense workers; kernel time alone cannot choose
the best complete path for very small matrices.

## Correctness and installation

Node041 completed ordinary small NSOC with Legacy/GVD/EVR/SLEPc; node042 completed
projected one-step SOC with those methods plus Auto. Full 57/114-state spectra
agreed within 4.89e-14 Ha, target20 within 3.00e-15 Ha. The existing after-force
projection and hash-verified native H/S handoff also apply to Legacy/full modes.
Dense methods return no wavefunctions and do not claim a runtime eigenpair audit;
finite spectra, native Hermiticity, certified positive overlap and independent
full-spectrum comparisons provide their checks.

A clean build compiled all 311 units, and a fresh projected-SOC installation
smoke confirmed Auto resolves to EVR and outputs all bands. The CPU suite passed
221 tests. The paired timing candidate SHA256 was
`2e1b5629b8dc61c8602695e8b61560d951759a1a64fd8fca3f6fd8bca8879bb8`;
the final change additionally selects the measured EVR method for explicit Auto.
The supplied standalone environment and ordinary OpenMX launch are unchanged.
