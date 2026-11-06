# Readable symmetry output and organized inputs

Put the selected symmetry report near the front of the normal System.Name.out,
with a provisional report available after input initialization and the same
report retained by the native final log merge. Print a readable table, integer
matrix entries where appropriate, and three-digit scientific notation for
residuals/noninteger values. Snap numerical rotation noise for display only;
do not modify metadata or any physical operator.

1. Add failing tests for pretty report formatting, initial .out creation and
   preservation at final log merge. Keep the upstream vendor snapshot unchanged.
2. Add a small report-copy header and a selected Merge_LogFile override; use the
   existing rank-zero bootstrap and input adapter. Off retains stock behavior.
3. Organize copied MoTe2 NSOC/SOC example inputs with restart, MD and symmetry
   beside the other controls before the coordinate block. Verify parsed values
   and every structural block match the archived originals.
4. Build on node040 and run the organized 9.43-degree, 222-atom example with
   64 MPI x 1, NSOC followed by its own one-step SOC. Inspect both .out files
   and compare H/S, density and energies with the existing same-projector run.
5. Run regressions and update the input/workflow documentation. Preserve the
   archived benchmark inputs, logs and hashes.
