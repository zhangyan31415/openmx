# Geometry tolerance and missing-symmetry diagnostics

```text
scf.Symmetry             on
scf.Symmetry.Tolerance   1.0e-3
```

`Tolerance` is a geometry matching upper bound in Å, not an SCF stopping threshold or an H/S covariance tolerance. With only `scf.Symmetry on`, the bound is5e-5Å. The supported range is finite positive values up to0.02Å. Larger values are not automatically more physically appropriate.

The automatic preparation scan now includes the exact user upper bound, with intervening probes above the previous1e-4Å cap. An irregular endpoint is included too. A sufficiently stable, validated operation set is selected; the endpoint does not force acceptance of a particular point group. Atomic coordinates and lattice vectors are not projected or rewritten. Setting a geometric tolerance alone does not explicitly request the covariance-controlled generic SOC route.

## When no usable subgroup is selected

The program exits before native SCF rather than silently changing the requested symmetry or using a larger tolerance. There is no interactive prompt inside an MPI job. The console prints the scan limit, a diagnostic classification and the diagnostic-file path. `symmetry-diagnostic.txt` and `SYMMETRY_DIAGNOSTIC.json` in the preparation directory retain details; no successful launch environment is written.

| Diagnostic | Meaning | Next action |
|---|---|---|
| `larger_tolerance_candidate` | A separate, bounded diagnostic-only scan found a validated supported subgroup above the requested limit. | Inspect the residual and whether the approximation is appropriate; explicitly increase `scf.Symmetry.Tolerance` and rerun if desired. |
| `identity_only_in_scanned_range` | The solve scan validated only the identity and the additional probes supplied no usable candidate. This does not prove symmetry is absent at every possible tolerance. | Inspect the structure; use `scf.Symmetry off` for an unrestricted calculation when no useful symmetry is expected. |
| `operation_validation_failed` | An operation or atom/periodic/orbital mapping failed validation. | Inspect the actual error and geometry; do not equate this with an asymmetric material or blindly loosen the tolerance. |
| `unsupported_subgroup` | Nontrivial operations were found, but none of the supported subgroups matches the request. | Check supported groups / the requested subgroup, or turn symmetry off. |
| `unstable_or_unusable_selection` | Candidates occur in the scan but the selected set is unstable or unusable. | Inspect the tolerance platform, residuals and mapping checks. |

Diagnostic-only probes use up to3×,10× and30× the upper bound, capped at0.02Å; these values do **not** become solve tolerances. A diagnostic candidate is not a proof of physical accuracy, a stable plateau, or solver eligibility after all later checks. At the0.02Å cap no larger probe is made. Missing dependencies or parser/runtime failures retain their own errors rather than being called “no symmetry”.

A finite-precision geometry check cannot in general distinguish an actually asymmetric structure from a perturbed symmetric one without an accuracy policy supplied by the user. We therefore report evidence and the scanned range, not a categorical physical diagnosis. User mode and explicit operation sets retain their stricter contracts.

The fix is in the Python preparation resources. Newly built installations include it; previously frozen benchmark binaries/resources are not modified. Existing installations must update their preparation resources through a fresh build/install to use the fix.

## Verification

`cpu/tests/test_symmetry_scan.py` checks exact bound coverage, custom/invalid ladders, diagnostic-only candidates, identity-only structures, validation failures and unsupported groups. Real three-atom fixtures test P1 failure and a perturbed D3 structure requiring a tolerance above the old cap. Successful and failing preparation leave the original input unchanged, and a failure does not create a native-launch environment. Existing bootstrap and input-contract tests are also run.

中文要点：只在用户允许的范围内选择有效对称群；范围外的探测仅用于提示。找到较宽容差候选时由用户决定是否增大；只找到恒等操作时提示可关闭对称性；映射失败、群不支持、平台不稳定分别报出，不伪称“材料没有对称性”。

## Visible solver selection

The ordinary `.out` summary now includes the prepared solver path, why it was selected, the geometry scan upper bound and the actual detection tolerance. Standard eligible D3/C3 cells use the specialized path; supported other cases can use general spinor blocks. Explicitly writing the public default covariance value has the same route as omitting it. A genuinely non-default advanced covariance contract retains the generic native checks and a visible reason; it is not silently ignored. This normalization does not certify arbitrary large-system memory reductions.

## Periodic coordinate representatives / 周期坐标代表

When OpenMX wraps periodic fractional coordinates into the reference cell,
native NSOC/SOC projector setup rebases the metadata image shifts to the
actual atom representatives before neighbor routes are built. This changes
only the image bookkeeping; the input file and physical geometry remain
unchanged. A displacement that is not an integer lattice translation still
fails the fixed-geometry check. This correction does not enlarge the requested
symmetry tolerance.

中文：无需手工折回周期等价坐标。NSOC/SOC 会在建立近邻映射前统一周期像
索引，保留输入和物理几何；真正的非整晶格位移仍会报错，不会通过放宽容差
掩盖坐标不一致。
