# OpenMX CPU：统一源码与运行入口

这个目录从源码构建一个同时包含 NSOC 和 SOC 的 CPU 程序，并提供同一套单节点、多节点运行与验收工具。它不依赖历史结果目录中的 OpenMX 对象文件、ELPA 对象文件或旧二进制。

当前计算范围是**固定几何、非磁性、Gamma 点，完整 NSOC 收敛后接自身检查点的一步 SOC**。一步 SOC 不等于 SOC 自洽计算。对称性检测成功、求解路径支持、数值通过、运行更快，分别记录，不能互相代替。

## 从这里开始

完整CPU代码目前位于 `experiments/cpu-d3-20261019` 分支，默认 `main`
仍保留GPU实现。新用户先取得正确分支：

```bash
git clone --branch experiments/cpu-d3-20261019 --single-branch \
  https://github.com/zhangyan31415/openmx.git openmx-cpu
cd openmx-cpu
conda create -n openmx-cpu python=3.12 numpy scipy spglib
conda activate openmx-cpu
```

Python环境用于元数据准备，不能替代外部Intel oneAPI classic编译器、MPI和MKL。
现有工作流还依赖Linux、可用的用户systemd服务、cgroup v2和SSH；不是跨平台
一键安装器。原生CPU源码构建和运行控制器的依赖分别见下面文档。

输入驱动的新用法见 [Auto/User 与统一 H/S 导出](docs/INPUT_SYMMETRY.zh-CN.md)：
`configure request.json --output /ABS/new-config` 自动准备匹配元数据和完整
profile，然后 `run /ABS/new-config/workflow.json --wait`。无需手动合并开关。
输入适配器限于统一 CPU 构建；裸二进制不会自行调用 Python 几何检测器。

1. [独立构建](docs/BUILD.zh-CN.md)：使用 `--profile unified`，生成二进制和构建清单。
2. [统一工作流](docs/WORKFLOW.zh-CN.md)：准备两个匹配输入，生成同一闭合子群的元数据，填写明确的运行配置，先 `plan` 再 `run --wait`。
3. [NSOC 源码与路径](docs/NSOC_INTEGRATION.md)、[旧 C3/D3 SOC 路径](docs/LEGACY_SOC.md)：核对所选路径的适用条件和导出语义。
4. [大角度单／多节点实测](docs/BENCHMARKS_20261031.md)、[验收方法](docs/ACCEPTANCE.zh-CN.md)：查看时间、内存、核时和数值边界。

```bash
python /ABS/REPO/cpu/build.py \
  --profile unified --oneapi-prefix /path/to/oneapi \
  --output /ABS/new-build --execute

cp cpu/tools/examples/input_request.json /ABS/request.json
# 按输入驱动说明填写两个输入、build-manifest、主机、CPU和运行库路径。
python /ABS/REPO/cpu/tools/workflow.py configure /ABS/request.json \
  --output /ABS/new-config

python /ABS/REPO/cpu/tools/workflow.py plan /ABS/new-config/workflow.json
python /ABS/REPO/cpu/tools/workflow.py run /ABS/new-config/workflow.json --wait
```

`prepare` 只生成元数据和环境增量，不会自行启用所有优化或启动计算。完整配置示例在 [`tools/examples/workflow.json`](tools/examples/workflow.json)；四条路径的配置模板及合并方法在 [PROFILE_TEMPLATES.md](tools/examples/PROFILE_TEMPLATES.md)。输入的基组、截断、电子温度、SCF 精度和混合历史不会被入口改写；`DATA.PATH` 必须明确给出。

低层 `geometry/prepare.py` 和显式profile留给高级复现。新用户优先使用上面的
`configure`，它自动生成完整配置。当前自动入口选择通用GSOC，不自动切换成
最快历史D3 legacy配置；不能把历史最快计时当成自动入口的测量值。

## 共用入口，按表示选择求解路径

| 群／阶段 | 实现 | 说明 |
|---|---|---|
| D3 NSOC | NSYM v3 | 保留等价表示块复用、可选合并返回和按需回变换 |
| C3 NSOC | NSYM v4 | 保留复共轭表示配对；不套用 D3 专用开关 |
| D2、C2 等 NSOC | 通用标量表示分块 | 按选定闭合子群构造表示，运行时仍检查 H/S 协变性 |
| 通用 SOC | GSOC 自旋表示 | 可显式选择投影和紧凑系数；支持 raw 或原生 symmetrize H/S 导出，含求解／导出算子一致性检查 |
| 旧 C3/D3 SOC | 专用已验证路径 | D3 保留 TR/REAL；C3 禁用 D3 专用 TR/REAL。原生导出投影须显式选择 `symmetrize` |

NSOC 的全程分块求解和实空间 H/S 投影是两个设置。`--project-nsoc` 明确请求在 SCF 中投影实空间 H/S，并导出投影后的 H/S，使用选定闭合操作集合；新群的生产认证仍需单独验收。相应 DM/EDM 钩子只检查，不写回投影密度矩阵。SOC 求解前的矩阵投影与最终导出的 H/S 是否投影也分别声明。几何容差不会自动变成电子结构误差许可。

单节点和多节点使用同一个程序；MPI 数变化后必须重新运行 NSOC，不能直接套用另一布局的重启文件。运行记录包含实际绑定、MPI 提供者、原生耗时、启动耗时、核时、每节点物理内存峰值和同步集群内存峰值。

`run --wait` 为整个流程保留连接，适用于没有用户服务常驻的节点。它不会修改用户的 Linger 设置，也不承诺客户端断开后继续计算。需要无人值守时，应在已有的持久作业脚本内运行此命令。

本目录不宣称任意点群、磁性、任意 k 点、几何优化或多步 SOC 均已优化，也不根据小体系测试承诺多节点线性加速。
