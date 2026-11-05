# OpenMX CPU：对称性优化

简体中文 · [English](README.md)

这是基于 [OpenMX](https://www.openmx-square.org/) 的 CPU 优化版本，主要改进对称性分块求解、密度矩阵计算及内存使用。同一个可执行文件包含 NSOC 和 SOC 路径，使用普通 OpenMX 输入文件和 MPI 启动方式。

**CPU 源码目前在 `experiments/cpu-d3-20261019` 分支。请按下面命令克隆；只克隆默认分支不会取得 `cpu/` 目录。**

## 1. 获取与安装

已验证的构建环境：Linux、Intel oneAPI 2022.3 classic 编译器、Intel MPI 2021.7.1、MKL 2022.2.1，以及 Python 3.12、NumPy、SciPy、spglib。构建器要求对应的 oneAPI 目录布局，详见[构建说明](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/BUILD.zh-CN.md)。

开启 `scf.Symmetry on` 需要 Python 环境中的 NumPy、SciPy 和 spglib。spglib 由 MPI 主进程在启动阶段调用，用于识别几何对称性，不参与逐步 SCF，也无需链接进 OpenMX 二进制。显式设置 `scf.Symmetry off` 时不调用这套 Python 工具；源码构建本身仍需要 Python。

```bash
git clone --branch experiments/cpu-d3-20261019 --single-branch \
  https://github.com/zhangyan31415/openmx.git openmx-cpu
cd openmx-cpu

conda create -n openmx-cpu python=3.12 numpy scipy spglib
conda activate openmx-cpu

# 改为实际的 oneAPI 安装路径；在计算节点上编译。
OPENMX_ONEAPI=/path/to/oneapi
source "$OPENMX_ONEAPI/setvars.sh"
python cpu/build.py --profile unified \
  --oneapi-prefix "$OPENMX_ONEAPI" \
  --output "$(pwd)/build/cpu" --jobs 8 --execute
```

`build/cpu` 必须是尚不存在的目录。构建会生成 `openmx-cpu`、配套的 `resources/`、日志和来源清单，OpenMX 与 ELPA 均从仓库源码编译。原子基组和赝势另行提供，在输入中通过 `DATA.PATH` 指定。

## 2. 开启对称性并运行

在现有 OpenMX 输入中加入：

```text
scf.Symmetry on
```

默认关闭，也可显式写 `scf.Symmetry off`。开启后，程序自动识别并验证可用子群、准备对称基并选择求解路径，无需用户生成 JSON 或单独运行准备工具。日志报告实际采用的群；映射或算子检查失败时明确停止。

在已加载上述 oneAPI 运行库和 Python 环境的计算节点上，照常运行，例如 **64 MPI × 1 线程**：

```bash
mpirun -np 64 ./build/cpu/openmx-cpu input.dat -nt 1
```

多节点仍使用同一个程序，通过站点的 MPI／调度器设置节点、进程数和绑定。移动安装时，须将二进制和旁边的 `resources/` 一起移动，并保留构建时的 Python 环境。

## 3. NSOC、一步 SOC 与 H/S 输出

目前验证的工作流是**固定结构、非磁性、Gamma 点：NSOC 收敛后，从其自身检查点接一步 SOC**。输入使用 `Cluster` 或 Gamma `Band`，以及 `scf.eigen.lib elpa1` 或 `elpa2`。NSOC 和 SOC 输入均加入 `scf.Symmetry on`，保持结构、基组、物理参数及 MPI 布局一致；SOC 使用原生 C2N 重启。具体输入及重启要求见[运行说明](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/WORKFLOW.zh-CN.md)。一步 SOC 不等于 SOC 自洽收敛。

需要 H/S 时另设普通参数 `HS.fileout on`。开启对称性后，受支持路径导出与所选群投影一致的 H/S。OpenMX 基组非正交，后处理必须使用匹配的 S。

## 4. 实现与验证范围

- D3/C3 包含专用优化；D2/C2 等受支持子群使用通用表示分块，实际路径由日志确认。
- 对称投影可能改变离散算子和 SCF 轨迹，需独立核对能量、力等物理误差。识别成功不保证任意材料都正确或更快。
- 当前直接开关通过了小体系单节点、双节点和安装目录搬移验证；MoTe₂ 9.43°、222 原子的 NSOC 完成，但通用 SOC 未通过默认协变门限。这例不能列为全流程成功。
- 历史加速和内存数据属于各自配置，不能当作当前入口的普遍保证。多节点也不保证核时下降或总内存不变。

[原生开关验收](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/NATIVE_SYMMETRY_VALIDATION_20261106.md) · [时间与内存基准](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/BENCHMARKS_20261031.md) · [数值验收标准](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/ACCEPTANCE.zh-CN.md)

## 代码目录与来源

| 目录 | 内容 |
|---|---|
| `cpu/vendor/` | OpenMX、ELPA 源码及许可 |
| `cpu/overrides/` | 优化求解器与内存内核 |
| `cpu/geometry/` | 对称性识别与表示构造 |
| `cpu/tools/` | 运行、对照与基准工具 |
| `cpu/tests/`, `cpu/docs/` | 测试与详细文档 |

[CPU 文档入口](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/README.md)。OpenMX 和 ELPA 的原作者版权与许可保留在源码及 `cpu/vendor/` 中；来源记录见 `cpu/vendor/ORIGINS.json`。

GPU 部分源自他人的 [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu)，相关说明请参阅原项目。
