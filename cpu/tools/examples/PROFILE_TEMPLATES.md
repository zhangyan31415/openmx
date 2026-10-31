# 完整配置模板的使用边界

`profile_templates.json` 收录四条已有实际运行记录的完整算法配置：

| 模板 | NSOC | 一步 SOC | 最终 SOC H/S 导出 |
|---|---|---|---|
| `D3_legacy_fast_candidate` | D3 v3 等价块与直接扇区路径 | 旧 D3 TR/REAL、compact、owner DM、HNL batch | 显式对称化 |
| `C3_gsoc` | C3 v4 分块 | 通用 GSOC compact | raw |
| `D2_gsoc` | 通用 D2 分块 | 通用 GSOC compact | raw |
| `C2_gsoc` | 通用 C2 分块 | 通用 GSOC compact | raw |

它们是复现实验和选择算法路径的模板。程序不会为新材料自动套用这些
投影配置，也没有把它们认证为适用于所有材料的生产默认值。D3 legacy
是复现历史较快 D3 路径的候选；通用 GSOC 仍是准备工具的默认 SOC 路径。

## 和当前元数据合并

先为自己的匹配 NSOC/SOC 输入准备元数据，明确选择同一个群、SOC 路径和
导出策略。D3 legacy 对应 `--soc-route legacy --soc-export-policy symmetrize`；
其他三个模板对应 `--soc-route generic --soc-export-policy raw`。
`--project-nsoc` 是否启用也必须与本次明确选择的 NSOC 投影策略一致。

模板中的 `${PREPARED.files.soc_basis}` 等是 **JSON 字段占位符**，需要替换为
自己 `PREPARED.json` 中的值；没有实现自动替换，也不是 shell 环境变量。
复制完整模板的 `profiles` 后，将当前元数据的环境增量覆盖进去：

```python
profiles = copy.deepcopy(template["profiles"])
profiles["nsoc"]["environment"].update(prepared["NSOC_environment_additions"])
profiles["soc"]["environment"].update(prepared["SOC_environment_additions"])
```

这会补齐本次 basis/map 路径、邻居 padding，以及 legacy 几何门限。
另外，将 `symmetry` 内的占位符替换为相同的实际路径；几何门限保留 JSON
数值类型。最终交给 `workflow.py plan` 前，配置中不能剩余任何占位符。
使用一套已匹配的元数据，避免把 C3/D3 或 legacy/GSOC 的文件混在一起。

完整单次流程应使用自己的全新收敛 NSOC，再由同一工作流从其检查点接
一步 c2n SOC。模板没有改写输入里的基组、电子温度、SCF 门限或 History。
本轮输入设置为：MoTe2/Bi2Se3 的 History=30、StartPulay=15；HgI2 的
History=35、NSOC StartPulay=15、一步 SOC StartPulay=50。SCF criterion 均为
`1e-7`。这些是记录的测试条件，不是新的材料无关收敛建议。

## 单节点与多节点

算法配置和 MPI 布局分别填写。模板列出 1 节点 64×1、2 节点每台 64×1、
4 节点每台 64×1，以及 4 节点每台 16×4 四种布局形状。总并行核数分别为
64、128、256、256。为每台机器显式填写实际获授权的 CPU IDs；单节点使用
shm，多节点使用工作流支持的 OFI 传输，并运行 plan 检查。

这些模板将私有 `OPENMX_CPU_NSOC_SECTOR_PARALLEL` 保持为 `0`，省略了
`OPENMX_CPU_NSOC_SECTOR_RANKS`。特别不要把旧的 `1,1,1` 或 64-rank 私有
分组复制给总共 128/256 ranks 的布局。D3 的直接扇区求解、等价块复用
与这种私有 MPI 分组是不同开关。

单节点 64×1 已有运行记录。布局示例自身不构成多节点加速和数值验收
证据；应以最终报告中完成的配对次数、完整耗时、核时和总内存为准。

每进程多线程时，必须显式设置 `OPENMX_CPU_VNA_SHARED=0` 和
`OPENMX_CPU_VNA_RAW0=0`：该共享缓存的原生保护条件要求一线程。
私有扇区并行仅接受四节点总 `256 MPI×1` 或 `64 MPI×4`；固定总
`64 MPI×1` 的四节点诊断必须关闭该开关。入口会在启动前检查这些组合。

## 已知物理边界

通用 SOC 的默认协变门限在模板里保留为 `3e-5`，没有遇到拒绝后自动放宽
的逻辑。MoTe2 和 Bi2Se3 的单独诊断曾显式使用 `1e-3` 完成运行，但它们
相对 raw 参考的最大力差仍未通过 `1 meV/Å` 门限：MoTe2 约 12.37、
Bi2Se3 约 1.10 meV/Å。Bi2Se3 该结构还显式使用了 `3e-4 Å` 的几何门限；
这个值同样没有推广为模板默认值。

NSOC 的 raw 参考对照也分别得到约 12.44 和 1.04 meV/Å 的力差。因此
模板明确保留 `raw_physical_force_gate=FAIL`；这是实际离散算子的投影误差
边界。legacy D3 不消费 GSOC 的协变门限，不能通过设置它改变 legacy
路径，也不能把通用 GSOC 的 SOC 力差当成 legacy SOC 的实测结果。

对于 legacy 的对称化 H/S，直接与 raw H/S 比较不应要求逐项为零；要用
同一群投影 raw 参考后比较。能量、力、电子数和完整流程一致性仍各自
验收。模板中的源配置 SHA、构建 SHA、测试布局和证据记录名称用于追溯；
运行不导入结果目录内的代码。
