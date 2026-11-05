# OpenMX CPU

[中文安装与使用](../README.md#中文) · [Installation and usage in English](../README.md#english)

## 中文

本目录包含 CPU 版 OpenMX 的源码构建、对称性优化、测试与运行工具。完整的中英文安装流程统一放在[仓库首页](../README.md)，请使用 `experiments/cpu-d3-20261019` 分支。

构建后的普通输入只需加入 `scf.Symmetry on`，再用 MPI 启动 `openmx-cpu`；默认关闭。程序内部完成对称性准备，无需单独生成 JSON。当前验证范围是固定结构、非磁性、Gamma 点，NSOC 收敛后从自身检查点接一步 SOC；具体材料仍须通过算子和物理误差检查。

- [构建与依赖](docs/BUILD.zh-CN.md)：编译器版本、源码构建、配套工具。
- [输入设置](docs/INPUT_SYMMETRY.zh-CN.md)：普通对称性开关及 H/S 输出。
- [运行与重启](docs/WORKFLOW.zh-CN.md)：NSOC、一步 SOC 及多节点布局。
- [原生开关验收](docs/NATIVE_SYMMETRY_VALIDATION_20261106.md)：已通过的检查及 SOC 适用性限制。
- [时间与内存基准](docs/BENCHMARKS_20261031.md)、[验收标准](docs/ACCEPTANCE.zh-CN.md)：已测配置与误差范围。

## English

This directory contains the CPU OpenMX build, symmetry optimizations, tests and workflow tools. The [repository README](../README.md#english) provides the complete bilingual installation guide. Use branch `experiments/cpu-d3-20261019`.

After building, add `scf.Symmetry on` to an ordinary input and launch `openmx-cpu` with MPI. Symmetry is off by default; preparation runs internally without user-written JSON. The tested scope is fixed geometry, nonmagnetic, Gamma-point NSOC followed by one SOC step from its own checkpoint. Each material still needs operator and physical-error checks.

Detailed guides are currently in Chinese:

- [Build and dependencies](docs/BUILD.zh-CN.md).
- [Input settings and H/S export](docs/INPUT_SYMMETRY.zh-CN.md).
- [Execution, restart and multi-node layouts](docs/WORKFLOW.zh-CN.md).
- [Native switch validation and SOC limitations](docs/NATIVE_SYMMETRY_VALIDATION_20261106.md).
- [Timing and memory benchmarks](docs/BENCHMARKS_20261031.md), [acceptance criteria](docs/ACCEPTANCE.zh-CN.md).
