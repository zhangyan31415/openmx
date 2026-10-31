# CPU科学与性能验收

运行检查与科学比较使用独立于历史结果目录的工具：

```sh
python cpu/tools/acceptance.py inspect /absolute/workflow/output
python cpu/tools/acceptance.py compare /absolute/reference/nsoc /absolute/candidate/nsoc --electrons 3404
python cpu/tools/acceptance.py compare /absolute/raw/soc /absolute/projected/soc --electrons 7232 --physical
```

需要NumPy。H/S读取器与验收器都在`cpu/tools/`内，运行时不导入结果目录的脚本。`hs_union.ORIGIN.json`记录已验证读取器的来源和原始SHA256；代码遵循仓库`LICENSE`的GPLv3许可。

默认比较检验实现等价：总能量差不超过`1e-8 Ha`，所有周期H/S块和导出的DM/iDM分量最大差不超过`1e-10`，最大原子力矢量差不超过`1e-3 meV/Å`。检查包括有限数、匹配的结构和基组顺序、电子数以及密度矩阵对角自旋块的反向邻居Hermiticity。SOC的全部实部和虚部Hamiltonian字段都保留，计算使用匹配的非正交S，不能以单位矩阵代替。

`--physical`明确用于原始算符与对称投影算符之间的物理精度比较，门限为`0.1 meV/原子`和最大力差`1 meV/Å`。这种较宽的物理门限不能用于掩盖同一实现应当等价的数值错误。工具完整保留H/S及密度矩阵差异，不把物理通过写成逐元素相等。

完整NSOC必须正常退出并实际收敛。投影路径采用原程序明确打印的`NormRD_projected`、dU和criterion，同时保留未投影NormRD：投影收敛不等于未投影自洽收敛。只达到步数上限却正常退出，不能计为完整验收。一步SOC是自身NSOC检查点上的一次计算，不表示SOC自洽收敛；检查点的MPI和线程布局必须相同。

`.scfout`不包含EDM或完整本征矢。直接比较EDM需要额外的原生`gsoc_density_rank*.bin`小体系审计；验收模块提供`compare_native_audit`，检查6个DM、4个EDM分量及`Tr(DS)=Σf`、`Tr(EDMS)=Σfε`。广义本征对残差和S正交性需要实际本征矢审计，不能由能量相等推断。

性能表必须分别列出完整NSOC、一步SOC、两者之和、物理核时、每节点内存和同步集群总内存。MPI RSS简单求和会重复计算共享页；节点独立峰值相加也不同于同一采样时刻的集群总峰值。缺少节点采样时不生成同步总峰值。阶段rank最大值和rank跨度用于诊断，不作为可直接相加的完整耗时拆分。未完成、失败和外部终止记录单独保留。

仓库测试使用自行生成的完整format12/15小文件，覆盖非正交S布局、全SOC字段保留、密度迹、截断文件拒绝、投影收敛与内存峰值口径：

```sh
python -m unittest discover -s cpu/tests -p test_acceptance.py -v
```
