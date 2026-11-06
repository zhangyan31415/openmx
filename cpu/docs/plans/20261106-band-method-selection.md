# Band method selection and full-spectrum comparison

Expose Band.Solver = Auto | LAPACK (GVD) | EVR | SLEPc | Legacy. Omitted solver
preserves the previous count-dependent behavior. Explicit Auto selects EVR for
all bands and SLEPc for a positive target count. Dense methods may also return
a selected subset, but report that the full spectrum was computed.

Reuse native sparse H/S snapshots and projected-SOC deferral. All independent
H/S band methods, including Legacy, can be deferred after force evaluation.
The dense worker uses matching nonorthogonal S, eigenvalues-only LAPACK, and a
safe lower bound on lambda_min(S) from the Frobenius norm of inverse Cholesky.
GVD calls the generalized driver; EVR uses Cholesky, HEGST and standard EVR.

Test analytic complex nonorthogonal problems, parser dispatch, native raw and
projected SOC, and the existing SLEPc path. Compare real-H/S principal blocks of
512/1536 and the complete 3848-AO MoTe2 matrix on node041-043, with paired driver
orders and 1/4/16 threads. Label principal-block kernel probes separately from
physical material calculations; record actual BLAS and memory scope.
