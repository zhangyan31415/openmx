# 独立 CPU 源码构建

此目录把原来散布在历史结果目录的 CPU 源码、OpenMX4.0 原生源码和随包 ELPA2018 源码收集到仓库内。默认统一构建编译全部311个翻译单元，不读取历史 OpenMX `.o`、ELPA `.o/.mod`、预编译库或已有 OpenMX 二进制。优化源码的原始版权声明保留，来源及 SHA256 见 `cpu/vendor/ORIGINS.json`，基础目标清单见 `cpu/overrides/soc/SOURCES.json`，统一NSOC/SOC选择清单见 `cpu/overrides/nsoc/UNIFIED_SELECTION.json`。OpenMX GPL 声明和 ELPA LGPLv3 声明保留在各源文件中；`cpu/vendor/COPYING.GPL-3` 和 `cpu/vendor/COPYING.LESSER-3` 提供许可文本。

外部依赖是编译器 Intel oneAPI2022.3 classic、Intel MPI2021.7.1、MKL2022.2.1（LP64、Intel MPI BLACS、线程 MKL；FFTW 兼容头/实现来自 MKL）。不是纯 Python 项目。默认外部前缀取环境变量 `ONEAPI_ROOT`，未设置时为 `/opt/intel/oneapi`；可用 `--oneapi-prefix` 更换同布局安装。运行时加载同一 oneAPI `setvars.sh`。对称性准备工具还需要 NumPy、SciPy 和 spglib；保留构建时 Python 环境及配套工具，原生程序会在启用对称性时调用它们。原子基组/赝势数据由运行输入的 `DATA.PATH` 显式提供，不打包进二进制。

在已授权计算节点上执行：

```bash
"$CONDA_PREFIX/bin/python" /ABS/openmx-dev/cpu/build.py \
  --profile unified --output /ABS/new-empty-build --jobs 15 --execute
```

`--profile` 默认为 `unified`，上面的示例显式写出以便复现。该目标将已验证的scalar NSOC优化、旧SOC及通用SOC路径编译到同一个二进制中；实际运行时由输入及程序内部的准备流程选择适用路径。

不写 `--execute` 则只生成并打印完整构建计划。`--output` 必须是不存在的目录，避免悄悄复用旧对象。C 源码最多 `--jobs` 个并行编译，Fortran 模块串行按依赖顺序构建；使用15时至多16个编译器进程。默认 CPU 目标是 `CORE-AVX2`；`--cpu-target HOST` 仅适合与构建节点兼容的处理器，并会被记录。

三个构建目标如下：

| 目标 | 翻译单元 | 用途 |
|---|---:|---|
| `unified`（默认） | 311 | 完整优化NSOC＋旧SOC/通用SOC，用户常规入口 |
| `soc`（显式可选） | 300 | SOC优化独立变体，未合入统一版scalar NSOC优化，适合已有检查点的SOC验证 |
| `stock` | 285 | 同一份已安装CPU4.0源码快照的原生/ELPA参考，不启用优化替换 |

`stock` 是可重复的同编译器参考，未独立证明等同于未经修改的上游发行包。`soc` 的源计划用于SOC路径部件对照，完整NSOC＋SOC流程应使用默认 `unified`。

`--profile unified` 使用 `cpu/overrides/nsoc/UNIFIED_SELECTION.json`，也可用 `--selection` 指定显式选择清单。该清单只有仓库相对路径：

```json
{
  "replacements": {"DFT.compact.o": "cpu/overrides/nsoc/source/DFT.unified.c"},
  "extra_units": [{"object": "new_helper.o", "source": "cpu/overrides/nsoc/source/new_helper.c", "language": "c", "flags": []}],
  "include_dirs": ["cpu/overrides/nsoc/source"],
  "capabilities": {"nsoc_symmetry": true}
}
```

不能用不存在的对象名做替换，也不能用仓库外部源码或 `..` 路径。共享对象必须通过源码合并/明确选择解决，不能把 NSOC/SOC 两份同名符号同时链接。

输出 `openmx-cpu`、`objects/`、逐对象编译日志、动态依赖日志，以及 `build-manifest.json`。清单包括源码/头文件 SHA256、全部编译链接命令、CPU目标、编译器/MPI版本、二进制 SHA256 和能力标记。`built_pending_native_validation` 仅表明独立构建完成；后续原生数值对照及多节点验收分别记录，不可由编译成功推断。

统一构建还输出 `resources/cpu/`，包含程序自动调用的几何准备工具，并在清单中
记录这些文件的 SHA256 和 Python 解释器。移动构建产物时将 `openmx-cpu` 与
旁边的 `resources/` 一起移动；只复制二进制不足以保证开启对称性后可运行。
不同计算节点须可访问相同的配套文件和 Python 环境。关闭对称性的普通路径
不调用该准备工具。
