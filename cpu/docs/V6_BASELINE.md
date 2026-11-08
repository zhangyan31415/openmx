# v6-clean CPU baseline / CPU 基线

## Release contents / 本次范围

The native CPU sources come from the validated v6-clean snapshot. The build
retains heap allocation of automatic Fortran work arrays, bounded D3 panel
communication with collective resource fallback, demand-based NUMA grid
replicas, and the validated cache/dense-buffer lifetime changes. These are
implementation changes; their eligibility still depends on the solver/input.
They are not a promise of equal speed or memory for arbitrary point groups.

原生 CPU 源码取自已验证的 v6-clean 快照，保留 Fortran 自动数组堆分配、
D3 分批有界通信及集体资源回退、按邻居需求分配的 NUMA 副本，以及缓存和
稠密矩阵生命周期优化。各路径仍有输入和求解器适用条件，并非任意点群均有相同收益。

The independently tested sparse-band Python fixes are included with their tests
and documentation. This does not change native SCF arithmetic. Later bilayer,
registry and direct-DM initial-guess experiments are not part of this candidate.
Existing internal startup hooks remain off; their presence is not a production
endorsement. No local calculation outputs, binary executables or private data trees
are included.

独立验证的稀疏能带 Python 修复随测试与说明保留，不改变原生 SCF 运算。
后续层间/局域堆叠/直接 DM 初猜实验不纳入；已有内部初猜接口保持关闭。
不包含计算输出、可执行文件或个人数据目录。

## Native validation / 原生验证

The original tested executable has SHA256:

```
6ed633854c4d27fe960797dc66c650c633674675c8246a89242ba5227823b610
```

A rebuild can have a different executable hash. This value identifies the
historical validation executable, rather than guaranteeing reproducible linking.

Two clean-up regression pairs for WSe2, 9.43 degrees, 222 atoms retained their
respective 37 and 46 SCF iterations. Maximum implementation differences were:
energy 1.28e-11 Ha; force-vector difference 8.96e-8 meV/angstrom; periodic H
1.25e-12 Ha; S 6.67e-16. These are equivalence checks against their respective
pre-cleanup configurations, **not** fresh comparisons against stock OpenMX.

两组 WSe2 9.43 度、222 原子的清理回归保留各自 37/46 步收敛轨迹。
上面的差值是相对各自清理前配置的实现等价性验证，不是新的原版加速测量。
完整 SOC、多节点以及所有材料并未在本次整理中重新验收。

## Settings and accuracy / 设置与精度

Use the existing build/run instructions and `scf.Symmetry on`. Preserve the
material's basis, physical settings and validated convergence controls. No
cutoff, history length or StartPulay value is universally substituted.

The measured low-memory NSOC profile additionally used these existing internal
runtime controls (record them when reproducing those measurements):

```bash
export OPENMX_CPU_NSOC_NATIVE_RETURN=3
export OPENMX_CPU_VNA_FORCE_STREAM=2
export OPENMX_CPU_GRID_RECOMPUTE=2
export OPENMX_CPU_NSOC_INITIAL_GRID=0
export OPENMX_CPU_NSOC_STARTUP=0
export OPENMX_CPU_NSOC_RHOK_MODE=0
```

These controls do not replace the complete input, runtime binding or eligibility
checks. Do not infer that `scf.Symmetry on` alone reproduces a particular measured
memory peak. See [SCF lifetimes](SCF_MEMORY_LIFETIME.md),
[NUMA replicas](NUMA_GRID_MEMORY.md), and [direct assembly](DIRECT_PANEL_MEMORY.md).

Recent 2.13393-degree, 4326-atom MoTe2/WSe2 tests with 80 Ry and criterion 1e-6
matched the sampled, globally aligned TAPW valence shapes to below 0.5 meV.
Their total-energy and force differences against the historical 150 Ry inputs
failed the user's 0.1 meV/atom and 1 meV/angstrom gates. These input experiments
therefore do not justify changing defaults or claiming structural-force accuracy.
Two-step z-grid-halving tests saved only about 6 percent of the measured physical
peak; they did not establish converged physical accuracy.

80 Ry 与 1e-6 的 2.13 度固定结构能带比较只验证了采样价带形状；总能量与力未
通过原门限。z 网格减半仅做两步试算，峰值内存收益约 6%，未验证收敛精度。
因此均不作为默认生产配置。计时、配置效应与代码收益分别记录。
