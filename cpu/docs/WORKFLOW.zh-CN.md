# 统一 CPU 工作流：单节点和多节点

入口为 `cpu/tools/workflow.py`。它对同一套源码构建的二进制执行完整 NSOC，然后复制这次计算自己的检查点，再运行恰好一步 SOC。两步使用同样的 MPI 数和线程布局；不会复用别的布局留下的重启文件。

当前适用范围为固定几何、非磁性、Gamma 点。几何识别、H/S 协变性、数值结果比较、完整性能是不同的验证；识别出点群不表示已经允许忽略非对角扇区。所有输入物理参数原样保存。输入驱动模式使用 `configure` 读取 `scf.Symmetry.Mode Auto/User` 并生成完整配置，见 [输入驱动说明](INPUT_SYMMETRY.zh-CN.md)。原有显式 profile 方式仍可使用 `Off`。

## 使用顺序

1. 由 `cpu/build.py` 从仓库源码编译，保留 `build-manifest.json`。运行入口检查二进制 SHA256、功能列表、外部对象依赖为空。
2. 准备 NSOC 和一步 SOC 两个输入。NSOC 使用 `scf.restart off`，允许正常收敛；SOC 使用 `nc/on/maxIter 1/c2n`，`scf.restart.filename` 等于本次 NSOC 的 `System.Name`。两者结构、基组、电荷、交换关联、截断、电子温度相同。输入必须给绝对 `DATA.PATH` 和 `HS.fileout on`。
3. 准备同一操作集合的标量和自旋对称基。已有验证资产可以直接使用；新资产用下面的 `prepare`。SOC 矩阵修复和紧凑存储默认关闭；`prepare --project-nsoc` 控制额外的 NSOC 实空间 H/S 投影，默认不加。准备元数据不会自行运行 OpenMX。
4. 完整显式 profile 和布局写入 JSON，先 `plan`，再 `run`。前台 SSH 保持连接，计算节点的用户 systemd 服务监督执行；所有主机另有任务专属 SSH 会话锚点，覆盖完整 NSOC＋SOC。`inspect` 检查服务和科学输出。入口不会申请或取消 Slurm 作业；这些 SSH 主机必须已经得到用户授权或已有资源分配。

先激活提供 numpy/scipy/spglib 的 Python 环境；入口本身使用标准库。下面以 `$CONDA_PREFIX` 表示已激活环境，非交互命令也可直接填写环境的绝对路径。`prepare --environment-prefix` 默认使用当前解释器的 `sys.prefix`：

```bash
"$CONDA_PREFIX/bin/python" /ABS/REPO/cpu/tools/workflow.py prepare \
  --nsoc-input /ABS/case/nsoc.dat --soc-input /ABS/case/soc.dat \
  --output /ABS/metadata-new --group D3 --geometry-tol 0.00005

"$CONDA_PREFIX/bin/python" /ABS/REPO/cpu/tools/workflow.py plan /ABS/config.json
"$CONDA_PREFIX/bin/python" /ABS/REPO/cpu/tools/workflow.py run /ABS/config.json --wait
"$CONDA_PREFIX/bin/python" /ABS/REPO/cpu/tools/workflow.py inspect /ABS/results-new
```

`plan` 只读：不创建目录、不 SSH、不启动任务。`run` 要求结果目录不存在；失败目录保留，不会自动覆盖、重试或静默切换为原始路径。必须加 `--wait` 并保持前台调用，直到本次受监督工作流退出。当前主机没有开启用户 Linger，异步启动会在最后的 SSH 会话退出后被系统结束，因此入口明确拒绝省略 `--wait` 的异步提交。**不承诺断开客户端后继续计算**，也不会擅自改变用户级 Linger 设置。

## JSON 接口 v1

完整例子见 [workflow.json](../tools/examples/workflow.json)。`profiles` 的环境必须来自本次源码构建实际支持并验证过的配置，不应直接混用历史二进制的开关。每一步配置 `name` 与完整 `environment`；入口清除继承的 `OPENMX_CPU_*`、MPI、OFI、OpenMP/BLAS 配置，然后设置所选 profile 和布局。环境中的 MPI transport 由 `runtime` 控制，profile 中不能夹带第二套 MPI 设置。

- `binary`：绝对可执行文件路径；`build_manifest`：同一二进制的源码构建清单。
- `inputs.nsoc/soc`：绝对输入路径；`output`：不存在的新目录，并且在源输入目录外。
- `layout.hosts`：1、2 或4个不同主机；`ranks_per_host`：每节点 MPI 数；`threads`：每进程线程数。
- `layout.cpu_ids_per_host`：每节点显式 CPU ID 列表，长度为每节点 MPI 数 × 线程数。按列表顺序给每个本地 rank 一段 CPU；不假设超线程或 NUMA 编号。需用实际物理核列表。当前 `NSOC_SECTOR_PARALLEL=1` 仅限4节点、总256 MPI×1线程或总64 MPI×4线程；固定总64 MPI×1等诊断布局必须显式关闭扇区独立进程组。原生扇区布局实现还要求每节点CPU ID集合恰好为0–63，且显式提供3个正整数扇区进程数、总和等于MPI总数。`NSOC_ELPA_OMP_THREADS=1` 只支持已启用扇区并行的64 MPI×4线程布局。
- `runtime.mpiexec`：绝对 Intel MPI 启动器；`runtime_prefix`：oneAPI 安装根；`environment_prefix`：Python 环境根；可给绝对 `library_paths`。
- `runtime.transport`：单节点可 `shm`，多节点默认 `shm:ofi`。可显式给 `provider`，实际提供者必须从 MPI 启动日志确认，不能仅根据环境变量归因。
- `symmetry.project/compact`：布尔值，默认 false；必须与显式 generic SOC profile 标志一致。NSOC 全程扇区求解由 `OPENMX_CPU_NSOC_FULL_PROJECT` 单独记录；真正的 NSOC 实空间 H/S 投影由 `OPENMX_CPU_REALSPACE_D3=scf` 表示。二者不能混称。实空间钩子会改写 H/S，但 `AfterDensity` 对 DM/EDM 仅作只读检查，不将投影结果写回。
- `symmetry.soc_export_policy`：默认 `raw`。旧 C3/D3 原生映射会在导出前投影 H/S，必须显式选择 `symmetrize`；通用 GSOC 支持 `raw` 或显式 `symmetrize`；后者要求 `GSOC_PROJECT=1`、`GSOC_EXPORT_SYMMETRIZE=1` 和匹配的 native map，并校验导出与求解投影算子一致。旧 D3 的 TR/REAL 优化要求已启用 D3 原生映射；C3 必须关闭这些 D3 专用开关。
- `monitor_interval_seconds`：默认0.25秒；`timeout_seconds`：默认 null（无内部执行时限），也允许省略或字符串 infinity；旧有限值报错。服务和 MPI 调用均不设置计算倒计时。可选 `memory_limit_bytes_per_host` 为每节点任务 slice 设置硬上限并禁用该 slice 的 swap；默认不设置上限。

`prepare` 输出的 `PREPARED.json` 包含同一闭合子群的 NSOC/SOC 对称基、真实空间操作映射、每个文件 SHA256 和两步环境增量。必须把这些增量合入各自完整 profile；增量本身不是完整生产配置。NSOC 生成器的 `OPENMX_CPU_NSOC_FULL_PROJECT=1` 表示全程使用分块求解路径，与额外实空间 H/S 投影开关分别记录。记录中的 `nsoc_realspace_HS_projection_explicit` 和 `nsoc_export_HS` 分别说明是否投影 H/S 及导出表示；启用 `scf`/`export` 时导出表示为原生标签 `PROJECTED_REALSPACE_H_S`。旧字段 `nsoc_density_projection_explicit` 已弃用并固定为 false，不再宣称显式密度投影。`--group auto` 是此元数据工具的检测选择，不能等同于输入解析器开启。

当前内置库路径对应 oneAPI22.3 / MPI2021.7.1；其他安装可提供 `library_paths`，需重新做 ABI 和实际库检查。环境前缀和输入路径在各节点共享。对称基、映射及其他环境中指定的已有输入文件在 plan 中转换为解析后的绝对路径，再冻结 SHA256；之后改变原始软链接不会让原生程序读取另一个文件。

## 验收与内存口径

每个 MPI rank 在该次任务的专属 systemd scope 内运行，所有同节点 scope 属于同一个 slice。它保留 Hydra 的 PMI 描述符；不会把 MPI rank 当成独立、脱离 MPI 的后台服务。实际主机、MPI rank、本地 rank、线程数、CPU affinity、cgroup 写入 `ranks/rankN.json`。`OPENMX_CPU_VNA_SHARED=1` 当前原生实现只支持每个 MPI 进程一个线程；`OPENMX_CPU_GRID_NUMA=1` 原生检查主线程的亲和性掩码只有一个 CPU；多线程时还需验证整个线程组的 NUMA 局部性，本入口暂按保守策略要求多线程布局显式设为0。VNA共享也必须设为0，`plan` 会提前拒绝不兼容组合，不会擅自替用户关闭。每个 rank 在原生程序启动前按自身 CPU 列表设置逐核 `OMP_PLACES` 和 `OMP_PROC_BIND=close`，清除冲突的 `KMP_*`、`GOMP_CPU_AFFINITY`；原生扇区求解器仍独立检查实际线程掩码。

每个节点独立采样 slice 的 cgroup v2 `memory.current` / `memory.peak`，共享页在节点内只计一次。各节点监控按共同墙钟采样 ID 对齐，报告：

- 每节点独立物理峰值；最大单节点峰值。
- 相同采样 ID、完整覆盖全部主机的同步集群峰值；这是离散采样观测值。即使只有一个节点，它也可能略低于内核 `memory.peak`；表格须区分“采样同步峰值”和“内核单节点峰值”。
- 各节点独立峰值之和，作为保守上界，**不冒充同步总峰值**。
- 节点 `MemAvailable` 和任务 cgroup 的 OOM 事件。

口径包含计算 rank 与小型监控服务，包含 cgroup 计入的缓存；远端 Hydra bootstrap 若位于 slice 外不包含在此任务内存。不会报告进程 RSS 相加作为物理总内存。

NSOC 只有正常退出、原生收敛标志、能量、有效 H/S 和完整同 MPI 数检查点全部具备才继续 SOC；SOC 检查自身 c2n 检查点读取成功及恰好一个 SCF 步。正常退出本身不代表数值与参考一致；完整能量、H/S、DM/EDM 等需要同一套基准验收程序比较。单步 SOC 不宣称 SOC 自洽收敛。

`workflow.json` 给出两步原生完整耗时之和、启动耗时之和、核时与内存。原生 `Elapsed.Time.` 与包含 MPI/systemd 启动的耗时同时保留。核时采用总物理核数 × 原生耗时；选布局时也应检查启动成本。所有终态和失败日志保留。

## 输出目录

临时运行文件默认位于共享结果目录父级的 `.openmx-scratch/`，按本次输出路径的
哈希分开。可用 `runtime.shared_scratch_root` 指定其他绝对共享临时目录。
多节点时，该目录与结果目录一样必须能从所有主机访问；没有固定用户名或
个人站点路径。

```text
results-new/
  config.json, plan.json, SUBMISSION.json, workflow.json
  controller.stdout.log, controller.stderr.log
  nsoc/   input.dat, run.json, stdout.log, stderr.log, *.out, *.scfout, *_rst/
  soc/    input.dat, run.json, stdout.log, stderr.log, *.out, *.scfout, own *_rst/
          RESTART_PROVENANCE.json, environment.json, LOADED_LIBRARY_RESOLUTION.txt
          ranks/rankN.json
          nodes/HOST/{ready.json,samples.jsonl,final.json}
```

失败和超时只停止由这份输出路径生成的任务专属 scope，关闭自身 SSH 会话锚点，不扫描或终止其他 OpenMX 进程。SSH 使用有限连接及心跳超时。输入、二进制和元数据在每个阶段、重启文件复制后以及 MPI 启动前再次校验；复制出的检查点逐文件核对 SHA256。不能更换布局继续用原来的 rank 分布检查点；应在新输出目录完整重跑 NSOC。

## 体网格输出与性能口径

体网格输出也是完整原生耗时的一部分。本轮冻结输入保留了原设置：
MoTe₂ 的 `level.of.fileout=0`，HgI₂ 为1；不能把跨材料的完整时间差全部归因于轨道数或点群。
HgI₂ 7.63°包含91,445,760个网格点，NSOC四份cube约4.816 GB，多节点下
`OutData`明显变慢。若实际工作只需要H/S，可以明确选择原生输入
`level.of.fileout 0`，同时保留`HS.fileout on`；这是输出配置，需在参考与候选中一致。
本轮没有把关闭原有cube输出的收益计入代码加速。
