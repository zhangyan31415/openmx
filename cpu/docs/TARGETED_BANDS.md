# Targeted native band dispersion / 原生目标能带

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

Put the controls near the front of the ordinary input:

```text
Band.dispersion       on
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
- Zero or omitted `Band.NumBands`: original full-spectrum band implementation.
- Each path segment needs at least two samples. Band count must be less than
  the matrix dimension minus one; a small buffer of extra converged states is
  requested internally to select the nearest states.
- This does not change SCF occupations or limit the number of occupied states
  used to form density. It does not guarantee equal numbers above and below the
  reference. In SOC mode the count includes both members of any returned pair.

中文：正数条数表示求参考能量附近最近的 N 个态，不是完整求解后裁剪。
默认参考费米能；数值参考用 eV 绝对能量。自旋极化的两个通道各返回 N 条；SOC
按完整 spinor 态计数。保留原有路径定义和普通 MPI 启动方式。

## Operator and parallelism / 算子和并行

The native code snapshots the actual distributed finite-range H/S blocks once,
including periodic translations and all SOC components. Independent workers
construct sparse Bloch H/S using the same positive Bloch phase and lower-triangle
Hermitian convention as the native full-band path. SLEPc solves the generalized
problem H(k)c = S(k)c epsilon by Krylov-Schur with shift-invert. No TAPW basis
truncation, identity-overlap replacement or overlap diagonal shift is used.

Workers are assigned to the first min(MPI ranks, k points) ranks, each using
the native `-nt` thread count. Their MPI/loader environment is isolated from the
parent oneAPI process. Files and the build-time Python must be visible on all
nodes. This is k-point parallelism, not one eigensolve distributed over all
OpenMX ranks. LU fill-in can still require substantial memory.

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
- Audits check finite values, raw Hermiticity defects, an extremal overlap
  eigenvalue, generalized eigenpair residuals and S orthogonality. Bases needing
  the stock overlap-penalty regularization are explicitly rejected in this mode.
- The targeted output is an energy window selection, not a guarantee of a fixed
  global band index or continuous subspace when that window changes along k.

当前限普通未收缩 AO 基，不支持 LNO。数值验收、阶段计时和完整 SCF 收敛分别报告；
只跑固定几步的求解器对照不能当作完整生产 SCF 的性能结果。

See [native checks and timing scope](TARGETED_BANDS_VALIDATION_20261106.md).
