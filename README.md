# OpenMX CPU: Symmetry optimization

[简体中文](README.zh-CN.md) · English

This CPU version of [OpenMX](https://www.openmx-square.org/) adds symmetry-based block solvers, density-matrix optimizations and memory improvements. One executable contains both NSOC and SOC paths and uses ordinary OpenMX inputs and MPI launches.

**The CPU source is currently on `experiments/cpu-d3-20261019`. Clone that branch explicitly; the default branch does not contain `cpu/`.**

## 1. Get and build

The tested build uses Linux, Intel oneAPI 2022.3 classic compilers, Intel MPI 2021.7.1 and MKL 2022.2.1. The supplied Python environment includes NumPy, SciPy, spglib and complex PETSc/SLEPc for targeted bands. The builder expects this oneAPI installation layout; see the [build guide](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/BUILD.zh-CN.md) for details.

Enabling `scf.Symmetry on` requires NumPy, SciPy and spglib in the Python environment. The MPI root process uses spglib during initialization to detect geometric symmetry; it is not called at every SCF iteration or linked into the OpenMX executable. Explicit `scf.Symmetry off` does not invoke these Python tools. The source builder itself still requires Python.

```bash
git clone --branch experiments/cpu-d3-20261019 --single-branch \
  https://github.com/zhangyan31415/openmx.git openmx-cpu
cd openmx-cpu

conda env create -f cpu/environment.yml
conda activate openmx-cpu

# Set the actual oneAPI installation path; compile on a compute node.
OPENMX_ONEAPI=/path/to/oneapi
source "$OPENMX_ONEAPI/setvars.sh"
python cpu/build.py --profile unified \
  --oneapi-prefix "$OPENMX_ONEAPI" \
  --output "$(pwd)/build/cpu" --jobs 8 --execute
```

The output directory must not already exist. It will contain `openmx-cpu`, its companion `resources/`, logs and a build manifest. OpenMX and bundled ELPA are compiled from repository sources. Supply basis sets and pseudopotentials separately through the input's `DATA.PATH`.

## 2. Enable symmetry and run

Place the symmetry setting beside the other SCF controls, before long coordinate blocks. For a fresh NSOC calculation:

```text
MD.Type                        Nomd
scf.SpinPolarization           off
scf.SpinOrbit.Coupling         off
scf.Symmetry                   on
scf.restart                    off
```

Symmetry is off by default; `scf.Symmetry off` disables it explicitly. When enabled, the program detects and validates a supported subgroup, prepares its symmetry basis, and automatically enables the verified solver and memory optimizations compatible with the input, thread count and rank binding. No separate preparation command or user-written JSON is needed. At startup, rank zero writes a readable symmetry table to standard output and the front of `System.Name.out`: selected group, operation count, geometry residual and each rotation/translation. Residuals use scientific notation; rounding matrix noise affects display only. Invalid mappings or failed solver checks stop the run; raw spatial covariance defects are diagnostic when the general SOC solver explicitly projects H/S.

With the matching oneAPI runtime and Python environment available on the compute node, launch normally, for example **64 MPI ranks × 1 thread**:

```bash
mpirun -np 64 ./build/cpu/openmx-cpu input.dat -nt 1
```

Multi-node runs use the same executable with the site's MPI or scheduler settings for hosts, processes and binding. Keep the executable and adjacent `resources/` together when moving an installation, and retain the build-time Python environment.

## 3. NSOC, one-step SOC and H/S export

The tested workflow is **fixed geometry, nonmagnetic, Gamma point: converge NSOC, then run one SOC step from that run's own checkpoint**. Use `Cluster` or Gamma `Band`, with `scf.eigen.lib elpa1` or `elpa2`. Add `scf.Symmetry on` to both inputs. Keep geometry, basis, physical settings and MPI layout consistent, and use native C2N restart for SOC. See the [run guide](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/WORKFLOW.zh-CN.md) for input and restart requirements. One SOC step is not a self-consistent SOC calculation.

Request H/S files with the ordinary `HS.fileout on` setting. Supported symmetry paths export H/S consistent with the selected group projection. OpenMX orbitals are nonorthogonal; post-processing requires the matching S.

For band dispersion, keep the usual `Band.Nkpath` and `<Band.kpath` definitions and add:

```text
Band.dispersion       on
Band.NumBands         20
Band.ReferenceEnergy  Fermi
```

This computes the 20 states nearest the reference at each k point using SLEPc, directly from native sparse H/S. `Fermi` is the default reference; a number specifies an absolute energy in **eV**. The normal `.Band` format and its **Hartree** energy units are retained. Zero or omitted `Band.NumBands` keeps the original full-spectrum path. Launch with the same `mpirun` command; no separate export, converter or moirekp installation is required. See the [band guide](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/TARGETED_BANDS.md).

## 4. Implementation and validation scope

- D3/C3 include specialized optimizations; other supported subgroups, including D2/C2, use general representation blocks. Check the log for the actual path.
- Symmetry projection can change the discretized operator and SCF trajectory. Energy, forces and other physical errors need separate validation; detection alone does not establish accuracy or speedup.
- The current default profile passed complete NSOC plus own one-step SOC for MoTe₂ at 9.43° (222 atoms) and 7.34° (366 atoms). Standard oriented C3/D3 slabs select the verified specialized SOC path; other supported cases use the general path and its operator checks.
- Historical timing and memory results apply to their recorded configurations. More nodes do not guarantee lower core-hours or constant total memory.

[Native switch validation](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/NATIVE_SYMMETRY_VALIDATION_20261106.md) · [Timing and memory benchmarks](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/BENCHMARKS_20261031.md) · [Acceptance criteria](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/ACCEPTANCE.zh-CN.md) (detailed reports are currently in Chinese).

## 5. Reference vs latest benchmarks

Fresh same-node measurements: **node040, 64 MPI × 1**, same compiler options and physical inputs. Both versions complete NSOC followed by one SOC step from their own checkpoints. Each version/size was measured once.

![Reference and latest CPU time and physical memory versus atom count](https://raw.githubusercontent.com/zhangyan31415/openmx/experiments/cpu-d3-20261019/cpu/docs/benchmarks/20261106/latest_vs_reference.png)

| MoTe₂ | Reference → latest total time | Reference → latest physical peak |
|---|---:|---:|
| 9.43°, 222 atoms | 184.86 → **86.42 s** | 73.59 → **25.22 GiB** |
| 7.34°, 366 atoms | 552.02 → **147.25 s** | 110.99 → **35.35 GiB** |

The reference is the bundled CPU4.0 source snapshot, built without optimized replacements. Symmetry constraints and SCF trajectories change together with the implementation; this is not a claim of identical raw operators. See the [full report](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/benchmarks/20261106/README.md), [data](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/benchmarks/20261106/measurements.csv) and [physical differences](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/docs/benchmarks/20261106/physical_differences.json), including force errors exceeding the earlier physical gate. The two-point lines are not an extrapolated scaling law.

## Source layout and credits

| Directory | Contents |
|---|---|
| `cpu/vendor/` | OpenMX and ELPA sources and licenses |
| `cpu/overrides/` | Optimized solvers and memory kernels |
| `cpu/geometry/` | Symmetry detection and representation construction |
| `cpu/tools/` | Workflow, comparison and benchmark tools |
| `cpu/tests/`, `cpu/docs/` | Tests and detailed documentation |

[CPU documentation](https://github.com/zhangyan31415/openmx/blob/experiments/cpu-d3-20261019/cpu/README.md). Original OpenMX and ELPA copyright notices and licenses are retained in the sources and `cpu/vendor/`; provenance is recorded in `cpu/vendor/ORIGINS.json`.

The GPU implementation comes from [dc1394/openmx4.0_gpu](https://github.com/dc1394/openmx4.0_gpu); please refer to that project.
