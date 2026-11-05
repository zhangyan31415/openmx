# 从输入启用 CPU 对称性

统一 CPU 构建现在提供输入驱动的 `configure` 入口。它读取两个输入中的
`scf.Symmetry.*`，生成匹配的标量／自旋对称基和完整运行配置；用户不需要
手动合并 profile、替换 JSON 占位符或填写求解器环境开关。

适用范围仍是固定几何、非磁性、Gamma 点，完整 NSOC 后接自身检查点的
一步 SOC。识别和求解支持不等于任意材料都通过了物理误差与性能验收。

## 输入

NSOC 和 SOC 两个输入使用相同的对称性设置，例如：

```text
scf.Symmetry.Mode           User
scf.Symmetry.PointGroup     D2
scf.Symmetry.Tolerance      1e-5
scf.Symmetry.HSExport       Symmetrize
scf.Symmetry.CovarianceTolerance 1e-8
```

`Tolerance` 单位是 Å，约束几何映射；`CovarianceTolerance` 是电子算子的
相对协变残差限值，两者不同。上面的 `1e-8` 是小体系集成测试的显式选择，
不是所有材料的推荐值。默认协变限值是 `1e-10`；只允许在明确选择
`HSExport Symmetrize` 时放宽，程序不在失败后自动改限值。

- `Mode Off`：保持原有普通输入和显式历史 profile 的行为。
- `Mode Auto`：在几何容差以内扫描，选择当前求解器支持的闭合子群。省略
  `PointGroup`。它报告实际选择的群；不能把所选子群冒充完整晶体点群。
- `Mode User`：必须指定群，使用给定唯一几何容差检测并验证该子群。
  `OperationsFile` 可给显式操作集合，集合必须与请求的受支持子群完全一致；
  不允许默默删操作或退回 P1。
- `HamiltonianPolicy Free`（默认）：不改写 SCF 各步的原始实空间 H/S。显式
  `HamiltonianPolicy Symmetrize` 另行请求 NSOC 各步投影，且必须同时选择
  `HSExport Symmetrize`；这会改变离散算子，须另行检查物理误差。
- `HSExport Raw`：保留原始实空间 H/S。它不表示“不使用求解器分块”。
- `HSExport Symmetrize`：NSOC 在原生导出时投影实空间 H/S；通用 SOC 在成功
  求解和原生能量／力计算后，使用轨道及 SU(2) 表示投影 H/S。DM 不在导出时
  重写。SOC 额外检查导出算子与求解所用 sector/TR/real 投影算子的一致性，
  不一致则拒绝写出 scfout。不能把投影误差当成舍入误差或自动放宽物理标准。

受支持的准备选项是 C2/D2、C3/D3、C4/D4、C6/D6、Ci、Cs；各自的生产认证
范围以实际验收为准。`GeometryPolicy Strict`、`KPointPolicy Auto`、
`DensityPolicy Free`、`ForcePolicy Free` 是当前支持的策略；其他值和未知
对称性关键字明确报错。

## 配置与启动

复制 [最小请求示例](../tools/examples/input_request.json)，填写二进制及其
构建清单、两个输入、全新结果目录、已授权主机和实际 CPU 布局。示例是
4 MPI×1 的小体系配置；生产核数须按实际硬件调整。无需 `profiles` 或
`symmetry` 字段，它们从输入生成。

```bash
"$CONDA_PREFIX/bin/python" /ABS/REPO/cpu/tools/workflow.py configure \
  /ABS/request.json --output /ABS/new-prepared-config

"$CONDA_PREFIX/bin/python" /ABS/REPO/cpu/tools/workflow.py run \
  /ABS/new-prepared-config/workflow.json --wait
```

`configure` 生成元数据和只读运行计划，不启动 MPI；应在计算节点执行几何
准备。它保留原输入，绑定输入及显式操作文件的 SHA256。修改输入或操作
文件后须重新 configure。`plan` 仍然不写文件、不启动计算。

新 CPU 二进制的原生输入适配器验证 Auto/User 和实际求解器配置一致。
因此直接对带 Auto/User 的输入运行裸 `mpirun openmx-cpu input.dat`，但没有
上述准备步骤，会明确报错；不会忽略请求然后继续原始 SCF。这里的自动
识别由统一 Python 入口实现，未把 Python 几何检测器嵌入 MPI 二进制。
这些变动不修改 HIP/GPU 输入解析器。

## 运行时限

入口、MPI 调用、监控及控制器都没有内部运行倒计时。systemd 服务明确设置
`RuntimeMaxSec=infinity`、启动／停止时限为 infinity，rank scope 也显式无限。
`timeout_seconds` 只能省略、写 `null` 或 `"infinity"`；旧有限值明确拒绝，
不能把旧的86400配置直接照搬。SSH 建连失败和断连检测不作为计算时限。
Slurm 分区或既有 allocation 的外部期限仍由调度器执行，本入口不会绕过它。

`run --wait` 需要保持连接，或在已有持久作业脚本中运行；无限运行时限不等于
客户端断开后服务一定存活。入口不更改用户 Linger 设置。
