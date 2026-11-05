# OpenMX CPU：源码构建与使用

这个目录从源码构建同时包含 NSOC 和 SOC 的 CPU 程序，不依赖历史结果目录中的对象文件或旧二进制。输入文件中加入 `scf.Symmetry on` 即可启用对称性；仍按普通 OpenMX 的方式启动。

当前验证范围是**固定几何、非磁性、Gamma 点，完整 NSOC 收敛后接自身检查点的一步 SOC**。一步 SOC 不等于 SOC 自洽计算。对称性检测成功、数值通过、运行更快分别记录，不能互相代替。

## 安装

完整 CPU 代码目前位于 `experiments/cpu-d3-20261019` 分支；根目录 CMake 的 CPU/GPU 构建尚未统一。先取得正确分支：

```bash
git clone --branch experiments/cpu-d3-20261019 --single-branch \
  https://github.com/zhangyan31415/openmx.git openmx-cpu
cd openmx-cpu
conda create -n openmx-cpu python=3.12 numpy scipy spglib
conda activate openmx-cpu
python cpu/build.py --profile unified \
  --oneapi-prefix /path/to/oneapi \
  --output /path/to/new-build --execute
```

需要外部 Intel oneAPI classic 编译器、Intel MPI 2021.7.1 和 MKL 2022.2.1，具体版本与目录布局见[独立构建](docs/BUILD.zh-CN.md)。当前不是跨平台一键安装器。原子基组和赝势由输入的 `DATA.PATH` 指定。

Python、NumPy、SciPy 和 spglib 用于程序内部调用的对称性准备工具；启用对称性时须保留构建所用的 Python 环境及配套工具。检测没有改成原生 C spglib 实现。普通 MPI 运行不依赖可选的 SSH/systemd 基准控制器。

## 使用

在普通 OpenMX 输入中加入：

```text
scf.Symmetry on
```

省略时默认为 `off`。程序内部识别并验证可用的闭合对称操作，准备对称基、分发元数据并选择求解路径；不需要单独准备 JSON 或手工合并环境开关。日志报告实际选择的群及校验结果。无法满足几何映射、求解支持或 H/S 检查时会明确停止。

加载与构建匹配的 Intel MPI/MKL 运行库后，照常运行，例如单节点 64 MPI×1：

```bash
mpirun -np 64 /path/to/new-build/openmx-cpu input.dat -nt 1
```

NSOC 和 SOC 输入都使用 `scf.Symmetry on`。SOC 必须读取本次 NSOC 自己的收敛检查点，保持相同结构、基组和 MPI 布局。普通 `HS.fileout on` 与对称性开启同时使用时导出相应群下的对称化 H/S。详见[输入说明](docs/INPUT_SYMMETRY.zh-CN.md)与[运行及重启](docs/WORKFLOW.zh-CN.md)。

多节点使用相同可执行文件，通过站点的普通 MPI 或调度器设置主机和绑定。MPI 数变化后应重新运行 NSOC，不能假定已有检查点跨布局兼容。程序没有内部计算倒计时；调度器外部期限仍有效。

## 实现与验证边界

| 群／阶段 | 实现 | 验证边界 |
|---|---|---|
| D3 NSOC | 等价表示块复用、合并返回和按需回变换 | 固定非磁 Gamma 路径 |
| C3 NSOC | 共轭表示配对及公共 DM/EDM 内核 | 不套用 D3 专有算法 |
| D2、C2 等 NSOC | 通用标量表示分块 | 运行时检查 H/S 协变性 |
| 通用 SOC | 自旋表示分块和一致的 H/S 导出 | 目前以自身 NSOC 后一步 SOC 验证 |
| 旧 C3/D3 SOC | 历史专用路径保留用于回归 | 历史最快配置不代表当前自动入口的计时 |

H/S 投影改变离散算子，必须与未对称参考独立比较物理误差；它不意味着 DM/EDM 被导出钩子改写。几何识别成功不能代替数值及性能验收，也不保证任意材料、磁性、任意 k 点、几何优化或多步 SOC 已支持。

- [大角度单／多节点实测](docs/BENCHMARKS_20261031.md)：时间、内存、核时和失败门限。
- [验收方法](docs/ACCEPTANCE.zh-CN.md)：实现等价与物理精度的不同标准。
- [历史输入及导出验证](docs/INPUT_SYMMETRY_VALIDATION_20261102.md)：小体系集成证据，不是当前入口的性能保证。
- [原生开关验收](docs/NATIVE_SYMMETRY_VALIDATION_20261106.md)：直接启动、双节点、整体搬移及当前SOC拒绝边界。
- [NSOC 实现记录](docs/NSOC_INTEGRATION.md)、[历史 SOC 算法记录](docs/LEGACY_SOC.md)：源码来源及算法范围。

历史 MoTe₂ 9.43°、222 原子同节点测试的优化配置，相对 CPU4.0 源码快照参考，完整 NSOC＋SOC 约快 2.05 倍、内存少 3.08 倍；它包含对称投影及收敛轨迹变化，不能都归因于实现优化。其历史最大力差约 12.44 meV/Å，未通过 1 meV/Å 的物理门限。多节点测试也没有达到核时不增长的目标，详见基准报告。
