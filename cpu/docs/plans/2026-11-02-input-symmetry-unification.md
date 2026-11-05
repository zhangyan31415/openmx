# Input symmetry and native export implementation plan

Goal: make Auto/User actionable through the CPU entry point, unify explicit raw/symmetrize SOC export, and remove internal calculation deadlines.

Architecture: keep validated scalar and spinor solvers. The preparation entry reads input symmetry settings, validates one common operation set, produces complete profiles automatically, and freezes their provenance. A native input adapter rejects Auto/User when matching prepared metadata is absent. The owner-block projector serves both legacy and generic SOC after native energy/force evaluation. No replicated dense matrix is introduced for export.

Authorization: the user requested these three changes on 2026-11-02; proceed without an additional permission round. Preserve all unrelated dirty files; no commit/push until implementation and acceptance are complete.

1. Reproduce time-limit defaults; test unlimited plans and long-lived remote calls. Remove internal run/service deadlines; reject stale finite timeout configurations. Keep connection-failure and explicit-stop cleanup separate from calculation duration.
2. Test input symmetry contract: Off compatibility, Auto supported subgroup selection, User exact subgroup/tolerance and operations validation, unsupported policies rejected. Add an input-driven configure command which produces complete runnable JSON instead of requiring manual template merging.
3. Add native CPU input adapter to verify selected profiles consume Auto/User requests. Do not modify the unrelated CPU/HIP parser worktree edits or claim HIP support.
4. Generalize owner-block native SOC export to selected closed operation sets; retain conjugation, reverse-block Hermiticity, SU(2), and improper-operation axial spin lifts. Add post-solve checksum/order guards. Raw export stays raw.
5. Run Python and native component regressions, independently rebuild the unified CPU target on an authorized compute node, then execute small actual NSOC + one-SOC checks for representative groups. Inspect H/S symmetry, electron count, energy/force preservation and final state. No production angle study.
6. Update public quickstart, export semantics, supported scope, unlimited runtime contract and verification evidence. Report any unvalidated scope explicitly.
