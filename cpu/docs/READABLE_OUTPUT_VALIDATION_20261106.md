# Readable symmetry output validation

Symmetry-enabled native input initialization writes a readable report into
`System.Name.out` and standard output. The provisional `.out` explicitly says
the calculation is in progress. Final native log merging retains exactly one
report before the echoed input and physical results. Symmetry-off behavior is
unchanged. The bundled vendor source is not modified; the unified build selects
a log-merger override with only a report-copy hook added.

Residuals use three-decimal scientific notation. Rotation/translation display
values within `5e-12` of an integer are printed as that integer; other values
use scientific notation. This formatting never changes the full-precision
geometry, maps, symmetry bases or solver matrices.

Verification: **210 CPU regression tests passed**, including native-adapter
initial output, final report placement, stock behavior with symmetry off and
display formatting.

An organized copy of the MoTe2 **9.43-degree, 222-atom, 3848-scalar-AO** inputs
was run on **node040, 64 MPI x 1**, using oneAPI 2022.3 and the unified build.
The setup is PBE, 150 Ry, 300 K, Gamma Cluster, criterion `1e-7`, History30,
StartPulay15 and D3 projection. NSOC, SOC and MD/restart controls precede the long
structural blocks. Parsed keyword values and all structural blocks are exactly
equal to the archived inputs; those archives were not rewritten.

The complete fresh NSOC converged in 39 steps. One SOC step read this new
NSOC's own checkpoint. Both exited normally and produced complete native H/S.
Compared with the previous same-projector workflow:

| Stage | Energy difference (Ha) | Maximum full H difference | Maximum S difference | Maximum exported DM difference |
|---|---:|---:|---:|---:|
| NSOC | -1.46e-11 | 6.65e-14 | 4.44e-16 | 1.84e-13 |
| One-step SOC | 3.64e-12 | 1.05e-13 | 4.44e-16 | 1.84e-13 |

Each final `.out` contains one complete six-operation report before the input,
and its report matches standard output. Native scientific output and checkpoint
reading remain valid. This is a presentation/implementation regression, not a
new symmetry-versus-raw accuracy study or a memory/speed benchmark.
