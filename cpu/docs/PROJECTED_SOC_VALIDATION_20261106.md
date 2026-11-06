# Projected SOC acceptance and startup report

Explicitly projected general SOC no longer rejects a raw spatial-sector or
real-gauge covariance defect merely for exceeding the raw covariance threshold.
Both measurements remain in the log. Raw, unprojected solves retain their gate.
Nonmagnetic time reversal, finite values, Hermiticity, geometry, basis Gram and
export-versus-solved-operator checks remain mandatory.

The final sampled eigenpair and S-orthogonality acceptance now uses exactly the
Hermitian sector/partner/real projection assembled by the solver, for both H
and S. Both relative residual and S-orthogonality tolerance are `1e-8`, independent
of the raw spatial covariance tolerance. Raw eigenpair diagnostics remain visible.
The independent action uses vector scratch, without a new full dense matrix.

## Verification

Python preparation and native adapter regressions plus the CPU suite: **205 tests
passed**. The numerical regression includes a Hermitian, time-reversal invariant
operator with broken sector covariance and matching nonidentity raw S. It checks
raw rejection, projected acceptance, nonfinite/TR rejection and deliberately
incorrect projected eigenvalues.

Native checks ran on **node040**, Intel oneAPI 2022.3 runtime, using the unified
source-built executable with the changed C units recompiled. All five one-step
SOC runs read existing matching NSOC checkpoints and exited normally:

| Case | Layout | Check |
|---|---|---|
| D2 HgI2 toy geometry, 3 atoms | 4 MPI × 1 | H/S unchanged; energy unchanged |
| C2 HgI2 toy geometry, 3 atoms | 4 MPI × 1 | H/S unchanged; energy unchanged |
| MoTe2 9.43°, 222 atoms, 3848 scalar AOs | 64 MPI × 1 | General SOC with default `1e-5` raw diagnostic threshold |
| Same MoTe2 input and checkpoint | 64 MPI × 1 | Previous binary control, diagnostic threshold `1e-3` |
| Same MoTe2 input and checkpoint | 64 MPI × 1 | New binary, diagnostic threshold `1e-3` |

The previous control requires the looser raw threshold solely to permit that
binary to reach its unchanged projected solve. This is a matched-operator
implementation comparison, not a proposal to relax production physics settings.
The new binary passes with either threshold and produces equivalent results.

For the MoTe2 case, raw H sector covariance is `4.65749e-4`, raw eigenpair residual
is `2.63354e-3`, and raw S-orthogonality defect is `2.64121e-8`. Those raw metrics
describe the difference from the projected problem. Against the actual solved
operator, eigenpair residual is `6.73066e-14`, S-orthogonality defect `3.77821e-15`,
and H/S export relative mismatch below `5e-15`.

Against the previous same-projector control, final energy agrees at printed
precision, full periodic H differs by at most `1.40e-14`, S by zero and exported
DM by less than `2e-14`. This does not establish physical equivalence to an
unprojected calculation. No new NSOC convergence or small-angle study was run.

Observed complete SOC wall time was 58.86 s for the previous control and
59.64/60.04 s for the new checks, about 1.3–2.0% higher in these single runs.
This is a correctness/usability fix, not a claimed speedup or a memory benchmark.

## Startup output

Rank zero prints `CPU_SYMMETRY_SUMMARY` with the selected subgroup, detected
point-group symbol, operation count and maximum geometry residual, followed by
`CPU_SYMMETRY_OPERATION` records containing fractional W/t, Cartesian R and
mapping residual. The convention is printed as `x_prime=W*x+t` for fractional
column vectors. NSOC and SOC each report once during input preparation. Full
atom permutations and periodic image shifts remain in `metadata/selected.symmetry.json`
under the preparation directory named in the log.
