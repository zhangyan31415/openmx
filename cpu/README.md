# OpenMX CPU

[English installation guide](../README.md) · [简体中文安装说明](../README.zh-CN.md)

This directory contains the CPU OpenMX build, symmetry optimizations, tests and workflow tools. The [repository README](../README.md) provides the English installation guide. Use branch `experiments/cpu-d3-20261019`.

After building, add `scf.Symmetry on` to an ordinary input and launch `openmx-cpu` with MPI. Symmetry is off by default; preparation runs internally without user-written JSON. The tested scope is fixed geometry, nonmagnetic, Gamma-point NSOC followed by one SOC step from its own checkpoint. Each material still needs operator and physical-error checks.

Detailed guides are currently in Chinese:

- [Build and dependencies](docs/BUILD.zh-CN.md).
- [Input settings and H/S export](docs/INPUT_SYMMETRY.zh-CN.md).
- [Execution, restart and multi-node layouts](docs/WORKFLOW.zh-CN.md).
- [Native switch validation and SOC limitations](docs/NATIVE_SYMMETRY_VALIDATION_20261106.md).
- [Timing and memory benchmarks](docs/BENCHMARKS_20261031.md), [acceptance criteria](docs/ACCEPTANCE.zh-CN.md).
