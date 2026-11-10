# Native band solvers / 原生能带求解方法

Install and build once using the repository's environment and normal oneAPI
toolchain. No moirekp package is used or required:

```bash
conda env create -f cpu/environment.yml
conda activate openmx-cpu
source /path/to/oneapi/setvars.sh
python cpu/build.py --profile unified --oneapi-prefix /path/to/oneapi \
  --output "$PWD/build/cpu" --jobs 8 --execute
mpirun -np 64 ./build/cpu/openmx-cpu input.dat -nt 1
```

The executable carries its companion helper in `resources/` and records the
build-time Python. Keep both the resource directory and dependency environment.
The user does not run a converter, Python helper or separate export command.
Complex PETSc/SLEPc and an available sequential sparse LU factorizer are checked
before SCF. Missing dependencies cause an explicit failure, not a full-spectrum
fallback. Ordinary SCF and the original full-band mode do not require SLEPc.

## Input / 输入

`Band.Solver` selects the band postprocessor independently of the SCF solver:

| Method | Meaning |
|---|---|
| `Auto` | EVR for all bands; SLEPc for a positive target count |
| `LAPACK` or `GVD` | Generalized LAPACK GVD, eigenvalues only |
| `EVR` | Cholesky + HEGST + standard LAPACK EVR, eigenvalues only |
| `SLEPc` | Sparse shift-invert; requires a positive target count |
| `Legacy` | Original native full-spectrum implementation; count must be zero |

For all bands, use `Band.Solver EVR` and `Band.NumBands 0`. Dense GVD/EVR may
also return a positive count near the reference, but still compute the full
spectrum internally; the audit records this. If `Band.Solver` is omitted,
count0 retains Legacy and positive count retains SLEPc. Explicit Auto uses the
table above. SLEPc needs PETSc; dense methods only need NumPy/SciPy.

中文：全带推荐显式写 `Band.Solver Auto` 或 `EVR`，并设置 `Band.NumBands 0`。
少量目标态可选 `SLEPc`。方法只控制能带后处理，不改变 SCF 求解器。


Put the controls near the front of the ordinary input:

```text
Band.dispersion       on
Band.Solver           Auto
Band.NumBands         20
Band.ReferenceEnergy  Fermi
Band.Nkpath           1
<Band.kpath
  61  0 0 0  0.5 0 0  G M
Band.kpath>
```

- Positive `Band.NumBands`: request that many energies per k point and spin
  channel. They are selected by distance from the reference and sorted by energy.
- `Band.ReferenceEnergy Fermi` (default): use the current SCF chemical potential.
  A numerical value, for example `-4.85`, is an **absolute energy in eV**.
- Zero or omitted `Band.NumBands`: return the full spectrum with the selected method.
- Each path segment needs at least two samples. Band count must be less than
  the matrix dimension minus one. SLEPc is asked for exactly the requested
  count (`search_nev = Band.NumBands`), with no extra-state padding or minimum
  search count. The Krylov working subspace uses `min(N,max(4*nev,20))` vectors
  to improve coverage of near-degenerate partners; this does not increase the
  requested or returned band count. The audit records its actual `search_ncv`.
  MUMPS keeps its parallel factorization and uses the conventional backsolve
  after factorization, avoiding a verified silent error in its tree-parallel
  solve path. BLAS thread settings are retained.
  When explicitly enabled, residual and orthogonality checks validate returned
  eigenpairs, not completeness of the nearest-state subset at a degenerate
  boundary. Full-spectrum selection differences are reported separately.
- This does not change SCF occupations or limit the number of occupied states
  used to form density. It does not guarantee equal numbers above and below the
  reference. In SOC mode the count includes both members of any returned pair.

中文：正数条数表示求参考能量附近最近的 N 个态，不是完整求解后裁剪。
默认参考费米能；数值参考用 eV 绝对能量。自旋极化的两个通道各返回 N 条；SOC
按完整 spinor 态计数。保留原有路径定义和普通 MPI 启动方式。

## Operator and parallelism / 算子和并行

The worker methods snapshot the actual distributed finite-range H/S blocks once,
including periodic translations and all SOC components. Independent workers
construct sparse Bloch H/S using the same positive Bloch phase and lower-triangle
Hermitian convention as the native full-band path. SLEPc solves the generalized
problem H(k)c = S(k)c epsilon by Krylov-Schur with shift-invert. No TAPW basis
truncation, identity-overlap replacement or overlap diagonal shift is used.

Hermiticity is enforced once by construction: scalar blocks and spin-diagonal
SOC blocks are built from the native lower triangle and conjugate counterparts;
the spin-off-diagonal block is retained and paired with its conjugate transpose.
Solver entry does not recheck and rebuild these CSR matrices. Compatible contiguous
CSR indices and values are borrowed for PETSc; incompatible dtypes are converted.
The input arrays remain alive until PETSc objects are destroyed. Strict audit mode
can measure Hermiticity without rewriting the supplied operator.

原生装配按下三角及其共轭关系一次生成厄米 H/S。求解入口直接复用 CSR，不重复计算
厄米误差或重建矩阵；类型兼容时不额外复制数组。SOC 的自旋非对角块不单独强制厄米。

Workers are assigned to the first min(MPI ranks, k points) ranks, each using
the native `-nt` thread count. Their MPI/loader environment is isolated from the
parent oneAPI process. Files and the build-time Python must be visible on all
nodes. This is k-point parallelism, not one eigensolve distributed over all
OpenMX ranks. LU fill-in can still require substantial memory.

The worker sets `MKL_DYNAMIC=TRUE` while retaining the requested thread ceiling.
It does not inherit a parent dense-SCF policy that disables dynamic MKL teams.
In the tested MKL2026/PETSc3.25.2/SLEPc3.25.1/MUMPS environment, disabling that
policy produced invalid eigenpairs for the 13000-order Si example; enabling it
passed independent residual, overlap-orthogonality and full-spectrum checks.
The native parent environment is unchanged.

SLEPc has no matrix-order threshold for overlap preprocessing. It keeps S sparse
and checks finite entries and a valid positive diagonal, without an extra dense
Cholesky factorization or minimum-overlap eigenvalue solve. Positive-definite S
is the generalized Hermitian solver's input requirement; these necessary input
checks do not certify its full spectrum, and the audit explicitly records that
no global positivity certificate was computed. Convergence, returned-state count
and finite-real-eigenvalue checks remain enabled. An ill-conditioned or indefinite
basis can still make the solver fail.

SLEPc 对所有矩阵规模使用同一条稀疏路径，不设置按阶数切换的门槛，也不额外求
S 的最小本征值。保留输入合法性、收敛状态、返回数量和有限实数检查。正定的 S
仍是广义厄米问题的输入要求。

The expensive post-solve generalized-residual and S-orthogonality audit is
**off by default**. Normal energy-only runs do not extract eigenvectors into
Python or recompute Hc/Sc. Developers can enable it with `OPENMX_BAND_AUDIT=1`
in the native launch environment, or `--audit` when invoking the helper directly.
The JSON records `eigenpair_audit_enabled`; residual and Sorth are null when the
audit was not performed, never reported as a successful zero-error check.
SLEPc still computes its internal Ritz vectors during iteration.

默认不提取波函数来复算广义残差和 S 正交性。调试时可设置
`OPENMX_BAND_AUDIT=1` 开启；关闭时 JSON 中这两项为 null，并记录开关状态。

Dense methods retain matching nonorthogonal S and compute no eigenvectors. They
check S by Cholesky and the rigorous lower bound
`lambda_min(S) >= 1 / ||inverse(L)||_F^2`; an ambiguous bound uses an extremal
overlap solve. No full S eigendecomposition is required. Dense audits do not
claim eigenpair residuals when eigenvectors were not computed. Every concurrent
k worker holds its own dense matrices and workspace; peak worker RSS is not a
node physical-memory measurement.

For symmetry-enabled one-step SOC, the supported generic projected solver is
selected conservatively. Its band stage is deferred until after native force,
total-energy and grid-output consumers. The normal one-shot real-space projection
and solved/exported-operator validation then run before sparse snapshots. The
later H/S file export accepts one hash-verified handoff of that same operator.
This also works without requesting an H/S file. Raw force consumers are not
exposed to prematurely projected tensors. The specialized legacy SOC memory
profile is not claimed certified with this new postprocessor.

## Outputs and checks / 输出与检查

- `System.Name.Band` retains the ordinary OpenMX layout; only its band count is
  reduced. Its eigenvalues and chemical potential remain in **Hartree**.
- `System.Name.Band.info` and the front of the final `.out` record the solver,
  selected count, reference energy, parallelism and work directory.
- `.openmx-bands-*` contains native block snapshots, worker logs and per-k JSON
  audits. It is created by the program; users do not prepare these files.
- All methods check finite values; Hermiticity diagnostics are opt-in. Dense methods
  also validate overlap conditioning. SLEPc performs sparse input checks and
  optionally validates generalized eigenpair residuals and S orthogonality; it does not
  certify the complete overlap spectrum. The worker does not apply the stock
  overlap-penalty regularization: detected invalid inputs or failed solves are
  reported explicitly.
- The targeted output is an energy window selection, not a guarantee of a fixed
  global band index or continuous subspace when that window changes along k.

当前限普通未收缩 AO 基，不支持 LNO。数值验收、阶段计时和完整 SCF 收敛分别报告；
只跑固定几步的求解器对照不能当作完整生产 SCF 的性能结果。

See [native checks and timing scope](TARGETED_BANDS_VALIDATION_20261106.md).

See [full-spectrum method comparisons](BAND_METHODS_VALIDATION_20261106.md).

### Floating-point roundoff and failures / 舍入误差与报错

SLEPc's basis orthogonalization uses the matching nonorthogonal S inner
product. Its imaginary-part roundoff guard is set to `256 * eps` for FP64
(about `5.7e-14`), independently of the eigenpair convergence tolerance.
This handles cancellation in the stock two-atom diamond example without
changing H, S, the target, or the requested band count. S must still satisfy
the positive-definite input contract. Materially complex overlap diagonals
are rejected; only FP64 roundoff is tolerated. Optional audit mode remains off
by default. Solver failures retain the original PETSc diagnostic, with any
cleanup errors reported separately; there is no silent dense fallback.

中文：内积检查允许浮点舍入范围内的微小虚部，不修改 H/S，也不改变目标
能量、带数或本征求解收敛精度。求解失败时保留首个数值错误，清理对象时
的后续错误单独记录，不再覆盖首个错误，也不会悄悄切换成全带求解。
