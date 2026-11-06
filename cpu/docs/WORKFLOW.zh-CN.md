# CPU 运行与重启：NSOC 加一步 SOC

使用同一个源码构建的 `openmx-cpu`，分别运行普通 NSOC 输入和一步 SOC 输入。两步保持相同 MPI／线程布局；SOC 读取这次 NSOC 自己的收敛检查点，不复用其他布局的重启文件。

## 输入开头的设置

把运行控制和电子结构设置放在坐标、晶胞等长数据块之前。下面是 NSOC 输入的控制区；后面接本体系的完整物种、晶胞和坐标，不要将这些开关追加到坐标末尾：

```text
System.Name                    openmx
DATA.PATH                      /path/to/DFT_DATA19
HS.fileout                     on
level.of.fileout                0
MD.Type                        Nomd

scf.XcType                     GGA-PBE
scf.SpinPolarization           off
scf.SpinOrbit.Coupling         off
scf.Symmetry                   on
scf.restart                    off
scf.EigenvalueSolver           Cluster
scf.eigen.lib                  elpa2
scf.maxIter                    100
scf.energycutoff               150
scf.ElectronicTemperature      300
scf.Mixing.Type                rmm-diis
scf.Mixing.History             30
scf.Mixing.StartPulay           15
scf.criterion                  1.0e-7
```

一步 SOC 在同一位置使用 `scf.SpinPolarization NC`、`scf.SpinOrbit.Coupling on`、`scf.restart c2n`、`scf.restart.filename openmx` 和 `scf.maxIter 1`，其余物理参数与 NSOC 保持一致。本例控制区的截断和混合参数是示例，实际材料使用自己的设置。

## 使用顺序

1. 按[构建说明](BUILD.zh-CN.md)编译，保留可执行文件及配套工具。启用对称性需要构建时的 Python 环境及 NumPy、SciPy、spglib。
2. 准备 NSOC 输入，使用正常的 SCF 收敛设置与 `scf.restart off`。加入 `scf.Symmetry on`；`DATA.PATH` 指向实际基组／赝势目录。需要 H/S 时设置 `HS.fileout on`。
3. 运行 NSOC，核对正常退出和实际收敛，不将达到步数上限误算为收敛。
4. 准备匹配的 SOC 输入：非共线、自旋轨道开启、最大迭代数为 1，采用原生 c2n 重启；`scf.restart.filename` 对应这次 NSOC 的 `System.Name`。结构、基组、电荷、交换关联、截断、电子温度保持一致，并同样加入 `scf.Symmetry on`。
5. 在与 NSOC 相同 MPI 布局下运行 SOC，确认读取了自己的 NSOC 检查点且恰好完成一个 SCF 步。一步 SOC 不表示 SOC 自洽收敛。

示例为单节点 64 MPI×1；用站点匹配的 MPI 启动器，在计算资源上执行：

```bash
mpirun -np 64 /path/to/build/openmx-cpu nsoc.dat -nt 1 > nsoc.log 2>&1
# 先核对 NSOC 收敛和检查点，再运行下一行。
mpirun -np 64 /path/to/build/openmx-cpu soc.dat -nt 1 > soc.log 2>&1
```

多节点主机、每节点进程数和线程绑定使用普通 MPI 或调度器配置。每个 MPI 进程的 OpenMP／BLAS 线程数应与 `-nt` 及分配的物理核一致，避免进程和线程嵌套超额。各节点须可访问相同的可执行文件、输入、基组和配套 Python 环境。准备和运行不要求用户手工生成 JSON 或元数据。

## 输出和错误

启用对称性后，普通 `System.Name.out` 前部和标准输出都有可读的对称性摘要：采用的子群、检测点群、操作数量、几何残差及每个操作的旋转／平移表。残差用三位小数科学计数法；矩阵显示将接近 0、±1 的舍入噪声归整，仅影响显示。完整精度仍保留在元数据中。

`.out` 在输入初始化后先写入对称性信息，并明确标记计算尚在进行；正常日志合并时保留该报告，再加入完整输入和物理结果。计算是否成功仍须检查原生正常结束及收敛状态，不能仅凭 `.out` 已存在判断。

程序详细日志报告求解路径和协变性检查。启用对称性且 `HS.fileout on` 时导出相应群下的 H/S；导出投影不代表将 DM/EDM 在导出阶段重写。失败的映射或电子算子检查会明确报错。

输入物理参数不会因启用对称性被自动放宽。但对称约束可能改变离散算子和 SCF 收敛轨迹，因此正常退出、布局等价和相对未对称参考的物理精度须分别检查。参见[验收方法](ACCEPTANCE.zh-CN.md)。

没有内部运行倒计时。外部调度器期限必须遵守；后台运行应由站点的持久批作业管理，普通会话断开后是否存活由站点环境决定。

## 时间与内存比较

完整比较列出 NSOC、一步 SOC、二者之和及核时。对比不同布局时，必须从各自完整 NSOC 的检查点启动 SOC；同时保持基组、精度、输出设置等一致。启动耗时和原生程序计时应分别记录。

物理内存应采用节点 cgroup 或匹配口径的节点监控：共享页在节点内只计一次，不能直接累加进程 RSS。分别报告每节点峰值、相同采样时刻的集群总峰值，以及独立节点峰值之和的保守上界；最后一项不能冒充同步总峰值。节点 `MemAvailable` 及 OOM 事件用于判断是否逼近节点容量。

可选基准工具保留额外的来源、绑定及内存记录，不是普通 OpenMX 使用的必需步骤。历史性能数据在[大角度基准](BENCHMARKS_20261031.md)，不能直接当作新材料或新布局的保证。

## 体网格输出

体网格输出属于完整原生耗时。只需要 H/S 时，可以使用普通输入 `level.of.fileout 0` 并保留 `HS.fileout on`；参考和候选必须保持相同输出配置。历史 HgI₂ 7.63°、678 原子基准含 91,445,760 个网格点及约 4.816 GB 的 NSOC cube 输出，多节点输出成本明显，不能把关闭这些输出的收益当作求解器加速。
