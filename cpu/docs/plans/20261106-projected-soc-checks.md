# Projected SOC checks and startup symmetry report

The public symmetry switch explicitly requests a projected H/S operator. Raw
spatial covariance and real-gauge defects therefore describe the projection,
not the correctness of the projected eigensolve. Keep these measurements in
the log without rejecting projected runs. Unprojected runs retain their gate.
Keep nonfinite, nonmagnetic time-reversal, geometry, basis, Hermiticity and
export checks. Audit projected eigenpairs and S orthogonality against the same
sector/TR/real operator used by assembly, using bounded vector scratch.

At startup, print the selected subgroup, detected point group, operation count,
geometry residual, fractional W/t and Cartesian rotation for each operation.
Print on rank zero once per stage and preserve the complete preparation logs.

Implementation and verification:

1. Add a numerical regression with a Hermitian, time-reversal invariant raw
   operator whose spatial covariance fails; verify projected acceptance, raw
   rejection and projected eigenpair checks, including a deliberately bad pair.
2. Add startup-report tests through the real Python preparation and C adapter.
3. Implement the policy and reuse the existing solved-operator matvec helper.
4. Build on an authorized compute node. Run tiny D2/C2 checks and a large-angle
   MoTe2 generic SOC case from an existing NSOC checkpoint. Compare projected
   H/S and energies with the unchanged projected operator where available.
5. Run the CPU regression suite, inspect native exits and expected artifacts,
   and document the validation scope. No new small-angle production runs.
