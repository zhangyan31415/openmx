# Unified standalone CPU build and multi-node acceptance

User authorization: unify verified paths, standalone CPU source build, common single/multi-node acceptance and several large-angle benches; direct SSH node040--043. User explicitly forbids commits until completion. Do not push. Preserve unrelated dirty files and existing worktrees.

Hardware observed: node040--043 each64 physical cores, Intel Xeon Gold6530, ~251GiB memory. Shared environment snapshots/home/mambaforge; native compiler/runtime oneAPI22.3. No resource scans or re-approval gates under the direct-SSH authorization.

## Objective and boundaries

Produce one CPU source target with the verified scalar NSOC symmetry/memory/fused-return path and legacy/general SOC paths; no external OpenMX .o files or prebuilt ELPA objects may be needed by a clean build. External compiler/MPI/BLAS/FFTW development libraries and explicit DATA.PATH are allowed. Build source manifests identify provenance, original copyright and exact selected implementations.

One documented workflow prepares geometry/group metadata, runs complete NSOC then one SOC from its own same-MPI checkpoint, on1/2/4hosts. Existing strict feature/covariance checks remain. Generic support and performance acceptance are separate; no arbitrary magnetic/k-point/MD promise.

## Work packages and ownership

A. cpu/build.py, cpu/vendor/, cpu/overrides/soc/, cpu/build/ documentation: recover exact tested SOC source + stock CPU/ELPA sources and include closure; compile every object freshly in a new build directory. No cached parent objects. Prove an empty-output clean build and dependency list with no result-tree source/object dependencies.

B. cpu/overrides/nsoc/, cpu/overrides/merge_stages.py and tests: recover best verified dc095 NSOC scalar/full-projection, D3 equivalent-block reuse, fused selected backtransform and safe memory paths; make deterministic checked native hook merge against A's SOC DFT. Preserve D3 fast path and arbitrary scalar sector fallback; C3/D2/C2 selection must use supplied subgroup operations, not silently redetect D3. Root owns final integration of source selection/overlap fixes.

C. cpu/tools/, cpu/tests/test_workflow*: common prepare/plan/run/inspect CLI; explicit source build binary capabilities, group/geometry/covariance contract, same-own-checkpoint NSOC->SOC,1/2/4host layouts, numeric options preserved, full environment/provenance, isolated outputs, per-node cgroup memory and simultaneous total peaks. No detached unsupervised processes. API should accept config JSON and new output path. Node040--043 only for this acceptance.

D. results/openmx_unified_cpu_20261031/bench/: prepare original large-angle datasets and source/hash refs, layout configs and strict science validator; no formal timed launch until root freezes binary/profile. Measure full NSOC+own SOC, ranks/threads, actual binding, per-node/cumulative physical memory, core-hours; successful expected artifacts, density trace and residual/Hermiticity checks where available. Preserve failures and anomalies.

## Sequence and checkpoints

1. Read-only source and provenance census. Freeze exact selected files before combining.
2. Mathematical/unit tests for new merge/launch behavior before implementation. Compile source-level stage integrations independently, then unified target. No commits.
3. Tiny nonmagnetic Gamma c2n smoke, raw/legacy/general routes as applicable; compare energy, force, full periodic H/S and native DM/EDM/electron metrics. Multi-node tiny cases test collectives/owner layouts, not speed.
4. Frozen original large-angle benches: MoTe2 9.43/7.34 where useful; HgI2 7.63; D3/C3 and D2/C2 through explicit validated groups. C3/C2 may be explicitly documented closed subgroups of the same structure. Keep criterion1e-7, original cutoff/temperature, History30 for MoTe2 and actual historical HgI2 settings unless a separate matched baseline is prepared. No reduced-state/energy/force bypass.
5. Matched layout comparisons:1x64x1,2x64x1,4x64x1 and4x16x4 (64/128/256 physical cores); the mixed layout is tested only if numerical guards and actual thread ownership pass. Fixed-total64 layout may be used as a diagnostic, not a speedup claim. Initial small component/probe before full runs to avoid known high-thread failures.
6. Use at least two independent full repeats for headline layout/candidate comparisons; baseline is clearly named (already-optimized historical versus raw4.0). Keep hardware/layout/code effects separate. Do not require an invented speedup threshold to finish; report honest tradeoffs. User's practical objectives: total CPU core-hours not markedly greater and total memory toward<=1.3x; failing layouts remain nonrecommended.
7. Review source portability and common workflow; rerun necessary checks after fixes. Complete bench artifacts and documentation first. Only then targeted final commit, preserving unrelated input-parser/old-run changes. No push unless separately requested.

## Completion requirements

- Clean CPU build uses repository sources and declared external libraries, zero parent OpenMX/ELPA objects.
- Common NSOC+SOC workflow works from own checkpoint on accepted single/multi-node layouts; incompatible modes fail clearly.
- Native numerical comparisons pass stated implementation gates or separately approved projection physical gates. A detected geometry alone does not license dropping H/S off-block elements.
- Report tables include group, N, hardware, layout, NSOC/SOC/full wall, core-hours, maximum per-node and synchronized total physical peaks, convergence and H/S/energy errors, failures and acceptance labels.
- No claim of all-point-group optimality, ideal scaling or several-fold whole-task memory unless measured.
- All task services/processes final states inspected. Commit only after all required work is complete.
