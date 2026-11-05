# 2026-11-02 输入对称性与导出接入验证

本轮是功能集成验证，不是生产算例或性能基准。当前分支
`experiments/cpu-d3-20261019`，记录基于 `c6776c34` 的 CPU 输入与导出接入变更。

## 构建与环境

- 源码：`/path/to/openmx`。
- Python 环境：`$CONDA_PREFIX`，Python 3.12.7。
- 构建：gpuh204，oneAPI2022.3，统一 CPU 目标311个源码翻译单元。
- 最终二进制：`snapshots/results/openmx_input_symmetry_20261102/build-final/openmx-cpu`。
- SHA256：`f14eea9dce033faf42717954390edced4efb6af7d4aaf32fb612209402eed7ec`。
- 最终构建清单保留了独立源码构建及输入适配器增量重编译记录；验收时静态源码闭包没有变化。

完整证据目录：`snapshots/results/openmx_input_symmetry_20261102`。

## 通过的原生计算

使用既有 HgI2 三原子测试结构，不代表大角度材料生产结果。结构来自
`results/openmx_unified_cpu_20261031/workflow_probe/raw3atom_two/nsoc/input.dat`。
原有基组 I7.0-s3p2d2 / Hg8.0-s3p2d2、150 Ry、300 K、SCF 1e-7、History35
和 StartPulay15 保留。几何容差1e-5 Å、显式协变限值1e-8。

| 输入／路径 | 原子数 | MPI布局 | 结果 |
|---|---:|---|---|
| User D2，导出对称化 | 3 | gpuh204，4×1 | NSOC＋自身一步SOC通过 |
| User C2，导出对称化 | 3 | gpuh204，4×1 | NSOC＋自身一步SOC通过 |
| Auto，实际选择D2 | 3 | gpuh204，4×1 | NSOC＋自身一步SOC通过 |
| User Ci，含不正操作 | 3 | gpuh204，4×1 | NSOC＋自身一步SOC通过 |
| User D2，双节点 | 3 | gpuh202/gpuh204，各2×1 | NSOC＋自身一步SOC通过 |
| User D2，NSOC各步投影 | 3 | gpuh204，4×1 | NSOC＋自身一步SOC通过 |
| Off，原始控制 | 3 | gpuh203，4×1 | NSOC＋自身一步SOC通过 |

各组使用自己的全新 NSOC 和匹配检查点，没有跨布局复用重启。
最终组命令、输入、节点、环境、布局、正常退出和文件记录分别在
`SMOKE_STATUS_release.json`、`inputs/`、`configured/`、`runs/` 下。
全部正向组与原始控制均正常结束，gpuh202/203/204 已检查无本轮原生进程残留。

## 数值检查

`VERIFICATION_release.json` 独立遍历全部周期邻居块，核对 AO/SU(2) 群变换后的 H/S，
最大元素不变性残差1.78e-15。不是仅比较两个随机向量。
原生程序另用两个确定性算子探针，比较实际求解的 sector/TR/real 投影与
即将导出的实空间群平均；它仍是采样检查，不冒充一般体系逐元素证明。

与同一最终二进制的 Off 控制比较，最大每原子能量差5.44e-8 meV/atom，
最大力差1.95e-5 meV/Å；通过0.1 meV/atom和1 meV/Å物理门限。
电子数 Tr(DS)、导出DM厄米性检查通过。单／双节点 D2 的 NSOC 和 SOC
均通过严格完整 H/S/DM 布局等价比较。EDM未写入scfout，本报告不宣称
由该文件完成了完整EDM比较。

裸 MPI 输入 Auto、但不提供匹配准备配置的负向测试，原生适配器以95退出，
在SCF前明确拒绝，日志见 `UNPREPARED_NATIVE_REJECTION.*`。

## 保留的失败及边界

初次配置误将D3专用 fused/selected-return 用于其他群，已修复并增加回归测试；
失败输出仍保留。另修复实空间元数据生成器假定固定D3二重轴方向的问题，
覆盖六角晶格两种合法C2轴族。

六角三原子 HgI2 额外检查中，D3/C3显式NSOC各步投影均完成；但一步SOC
原始H的相对协变残差分别约1.05e-3和9.10e-4，超过声明的1e-8限值，
明确拒绝继续。没有自动放宽阈值。这些结果不能列为成功的通用SOC生产验证。
对应证据在 `SMOKE_STATUS_hex_projected.json` 及同名 runs 中。

新增数学/C99对照覆盖C2、D2、C3、D3（含s+p轨道和重表示），不正操作
自旋使用轴矢量 det(R)R 的SU(2) lift。完整 CPU 测试137项通过，见
`ALL_TESTS_FINAL.log`。14项运行时限回归确认无内部运行倒计时，旧有限
配置被拒绝；无限内部时限不取消调度器外部限制或客户端连接生命周期。

本轮不声称性能或内存进一步改善，也不把小体系精度通过推广到所有材料。
