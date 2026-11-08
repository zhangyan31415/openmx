# Changes since c260a3e2 / 相对已发布 CPU 分支的改动

Baseline: `c260a3e225c9eb74ded33c7d4c35bc793c8335f9`, the GitHub
`experiments/cpu-d3-20261019` head when this update was prepared.
These changes do not describe a comparison against the unrelated GPU main branch.

| Area | Change |
|---|---|
| Direct D3 assembly | Bound send/receive batches; budget planning and collective resource fallback preserve D3 solving instead of a fixed 64 MiB failure. |
| ELPA build | Put automatic Fortran work arrays on the heap to avoid relying on an 8 MiB process stack. |
| NSOC storage | Retain validated native coefficient handling, packed-overlap page reclamation, dense/sparse buffer lifetimes and final history release. |
| Integration caches | Retain demand-based NUMA replicas, VNA value/derivative lifetime controls, force streaming and grid recomputation paths; eligibility and explicit runtime controls still matter. |
| Reciprocal mixing | Retain the clean low-memory compatibility path and the native LU scratch-size fix; no default switch to an experimental density guess or adaptive polishing. |
| Sparse bands | Remove mandatory expensive eigenpair/overlap audits and redundant Hermitian rebuilding; retain opt-in diagnostics and the tested sparse-LU backsolve refinement. No change to the generalized nonorthogonal eigenproblem. |
| User documentation | Bilingual README for develop, memory-path notes, v6 validation scope, cutoff band figures and small numerical datasets. |

本次以 v6-clean 冻结源码为原生基线，包含上表的构建、通信和内存修改。
这些修改并不是本次重新实现的算法，也不代表对任意材料都有同样加速。
后续层间/堆叠密度初猜实验没有合入；实验开关保持关闭。
完整的原生数值证据来自该冻结版本的已有回归，本次发布整理另执行 Python
测试及构建计划检查，没有重新运行整套 SCF/SOC 或多节点性能矩阵。

See [baseline/settings](V6_BASELINE.md) and
[cutoff results](benchmarks/20261109-cutoff/README.md).
