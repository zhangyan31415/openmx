# Force3 compatibility with NVHPC 24.1

The OpenMP Force3 implementation used a thread-private array of four potential
buffer pointers. On the tested NVHPC 24.1 toolchain, stores through that array
were compiled as writes into the thread's stack. The resulting corruption could
appear later as a null `work_dObs.RF` dereference during orbital evaluation.

The source workaround replaces the pointer array with four scalar pointers.
It preserves OpenMP, GPU selection, grid layout, synchronization and the
2/2/4/4 noncollinear potential weights. Spin modes 0 and 1 keep their original
4 and 2 weights. Components absent from a spin mode are not accessed.

## Independent compiler reproducer

`tests/nvc_omp_pointer_array_bug.c` contains no OpenMX, MPI or GPU calls. Its
shared grid buffers have disjoint per-thread writes and are checked after the
parallel region. These are exactly representable integer-valued test data.

```sh
nvc -O3 -Mnovect -fopenmp tests/nvc_omp_pointer_array_bug.c -o repro
OMP_NUM_THREADS=1 ./repro 3
nvc -O3 -Mnovect -fopenmp -DUSE_SCALAR tests/nvc_omp_pointer_array_bug.c -o scalar
OMP_NUM_THREADS=1 ./scalar 3
```

On the tested AMD EPYC 7282 system, the NVHPC 24.1 array variant crashes for
spin 0, 1 and 3. Its scalar and no-OpenMP variants pass. NVHPC 24.7 and GCC
11.5 pass the array tests. The scalar variant also passes all three spin modes
with eight OpenMP threads under both NVIDIA compiler versions.

Adding `volatile` to the pointer array did not repair the full Force3 case.
Do not replace this workaround with that qualifier or disable Force OpenMP as
an assumed equivalent performance fix.

The reproducer establishes a compiler-dependent OpenMP code-generation issue;
it does not identify the particular internal compiler pass. The observed
failure persists with an effective trailing `-Mnovect`.

## Full-application regression

`tests/run_force3_nvhpc.py` runs a saved 24- or 96-atom MoTe2 SOC Gamma input
against an independent saved CPU result. Supply absolute paths and a new output
directory, inside an authorized GPU allocation after loading the intended
compiler/MPI environment:

```sh
python tests/run_force3_nvhpc.py \
  --binary /absolute/path/to/openmx \
  --input /absolute/path/to/input.dat \
  --reference /absolute/path/to/cpu_reference \
  --work /absolute/path/to/new_output_directory \
  --atoms 24 --threads 8 --force3-gpu 0
```

The case uses `mote2_<atoms>` as its OpenMX system name and must write its
noncollinear eigenvalues, forces and HS export. The test requires normal
completion, a completed Force3 profile, finite values, matching eigenvalue and
force counts, and maximum energy/eigenvalue/force differences no greater than
1e-7 atomic units. The HS check establishes that an export exists, not matrix
parity. `--all-force-cpu` additionally checks that Force.c launched no device
kernels; `--force3-gpu 1` requires an actual Force3 GPU launch. Logs include
commands, input/binary hashes, loaded runtime paths, GPU samples and thread
CPU masks. OpenMP binding operates within the CPU set inherited from Slurm.

The server-specific evidence directory is
`snapshots/openmx_force3_fix_20261018`.
Its report and summary distinguish expected failing controls from the fixed
application runs. The source structure is bulk `mp_poscar/mp-602.vasp`, using
Mo7.0-s3p2d2 and Te7.0-s3p2d2f1, SOC, Gamma, native FP64 and three SCF iterations.
These checks establish smoke-test parity, not SCF convergence, twisted-system
validation, performance scaling or validation of thousands of atoms.

The diagnostic A100 executables link the frozen upstream objects at
`d7719488d298f8ddef31610df800bf118fb7988b` with only Force.o replaced. They do
not incorporate the separate local rank-distributed DM/EDM acceleration. The
source fix in this repository can coexist with that code, but the diagnostic
executables must not be described as combined acceleration benchmarks.

The original external server's failing input and original NVHPC 24.1 CUDA/MPI
bundle were not available. This local comparison uses CUDA 12.8 and common
OpenMPI 4.1.5, with each compiler's matching OpenMP/OpenACC runtime.
