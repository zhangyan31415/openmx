# Legacy C3/D3 SOC adapter

`cpu/geometry/legacy_soc.py` prepares the **existing verified** sparse spinor
basis consumed by `cpu_c3.recv.c`, plus the matching native real-space SOC map.
The ordinary unified GSOC path remains the general default. This adapter is an
explicit route for comparing/reusing the earlier optimized legacy D3 path.
It does not introduce a new solver, truncate states, or modify coordinates.

The common preparation entry point selects it explicitly:

```bash
"$CONDA_PREFIX/bin/python" cpu/geometry/prepare.py \
 --nsoc-input /absolute/NSOC.dat --soc-input /absolute/SOC.dat \
 --output /absolute/new-metadata --group D3 \
 --soc-route legacy --soc-export-policy symmetrize
```

Omitting the last two options retains the existing generic GSOC/raw-export
behavior. `--export-policy` is an alias for `--soc-export-policy`. Common
`PREPARED.json` now records `SOC_route`, `SOC_export_policy`, `SOC_export_HS`
and hashes all selected SOC artifacts; legacy metadata also includes the
`soc_native_map`. These declarations do not activate an external process.

## API

```python
result = legacy_soc.prepare(
    input_soc="/absolute/SOC.dat",
    symmetry_json="/absolute/selected.symmetry.json",
    output_dir="/absolute/new-directory",
    group="D3",                    # D3 or C3 only
    export_policy="symmetrize",     # raw is the default
)
```

The caller must load `legacy_soc.py` from this repository, use the selected
closed subgroup produced by the unified preparation workflow, and merge
`result["env_additions"]` into the complete **clean** SOC runtime profile. These
additions alone are not a complete launch configuration. In particular, the
verified D3 antiunitary path needs `OPENMX_CPU_SOC_COMPACT=1`, distributed DM
settings, and its own successfully converged nonmagnetic NSOC checkpoint with
the same restart rank layout. Unset inherited `OPENMX_CPU_*` flags before
applying the selected profile. Never combine legacy and GSOC activation.

The complete portable D3 algorithm profile is frozen in
`cpu/geometry/legacy_helpers/d3_fast_profile.json`. Its source is the verified
2026-10-25 `select_on1_soc` run with binary SHA256
`9ba1dd686da8abef7d45a76e19dcd91aa49ac6e52a4a18f9e50004e8a4f59acc`;
the JSON records the source-profile hash. It contains no result-tree paths and
does not import the original run at runtime. Load its `environment`, then
overlay the current prepared `env_additions`. Geometry/basis/native-map paths
come from that preparation, and the workflow sets MPI transport/threading and
binding independently. Do not copy historical node lists, thread counts or
transport settings into this algorithm profile.

The full D3 profile restores `SOC_COMPACT=1`, TR/REAL, owner-sparse DM,
local H/S and overlap, shared force DM, and DM SORT/TRIM/CYCLIC, with BLAS3
density and HNL batching enabled. New GSOC selectors and scalar NSOC solver
selectors are explicitly disabled for the SOC stage. Using only the basis/map
additions with a generic GSOC memory profile would lose this tested path.
The portable profile is for D3 only; it is not a C3 antiunitary profile.

Returned fields include `basis`, `map`, `files`, `env_additions`, `export_HS`,
SHA256 hashes, orbital dimensions, geometry residual/bound and algebra checks.
`LEGACY_SOC.json` records the same information. Existing output directories are
rejected. Only local NumPy/SciPy and the repository geometry detector/orbital
tables are imported; installed OpenMX source and tapw are not used.

## Export policy is explicit

| Policy | Native H/S export | Legacy TR/real optimization |
|---|---|---|
| `raw` (default) | Raw H/S | Disabled |
| `symmetrize`, D3 | Internally symmetrized spinor H/S | TR and real enabled |
| `symmetrize`, C3 | Internally symmetrized spinor H/S | TR/real disabled |

In the preserved legacy implementation, `SOC_TR=1` requires the active native
C3/D3 validation hook; that hook also projects H/S before export. There is no
independent raw-export switch for that fast D3 configuration. The adapter does
not silently change raw exports to symmetrized ones. The legacy antiunitary
implementation only supports the four-sector D3 basis; setting `SOC_TR=1` on
the three-sector C3 basis would abort. General GSOC supplies the separately
implemented C3 time-reversal route.

The legacy solver itself solves projected sectors even with `raw` exports.
Consequently raw H/S exports describe the original assembled operators, not
an assertion that those operators are exactly covariant. Geometry and basis
preparation checks do not replace native H/S covariance/physical-error checks.

## Scope

Fixed geometry, nonmagnetic Gamma NSOC charge followed by exactly one c2n
Cluster SOC step, in the established hexagonal slab orientation. Unsupported
rotated cells must use the general GSOC adapter. Runtime guards retain the
restrictions on contraction, Hubbard terms, fields, and auxiliary states.
The adapter revalidates integer atom/image actions and group closure against
the SOC geometry, then uses the exact historical +120-degree z convention.
It preserves the selected source C2 gauge. It does not redetect or enlarge the
requested subgroup and independently checks the nominal-rotation residual.

## Verification

On node041 with `$CONDA_PREFIX` and one BLAS/OpenMP thread:

```bash
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
 "$CONDA_PREFIX/bin/python" -m unittest discover \
 -s cpu/tests -p test_legacy_soc.py -v
```

Seven adapter tests check C3 cycle/fixed atoms, D3 partner and +/-i gauges, malformed
permutations, binary headers, explicit export policy, and one-step/output
boundaries, plus the portable full D3 profile. These are metadata tests, not
native SCF performance measurements.
Four additional `test_legacy_soc_integration.py` tests verify common artifact
selection, unchanged generic defaults, explicit raw/projected export semantics,
and rejection of unsupported routes/groups. All eleven pass on node042 with the
same one-thread runtime; node041 was not used for this integration check.

`LEGACY_SOC_MATCHED_GAUGE.json` records a further comparison with the exact
MoTe2 7.34-degree basis used by the verified 2026-10-25 fused-return runs:
the regenerated basis is byte-identical, SHA256
`f605c9ff3a8ee254d9680bdf81c7725ad636a20433b12dacdc701e7cbca7468e`.
The matching SOC map has identical integer actions, spin and orbital rotation
blocks; parsed Cartesian coordinates differ by at most 7.11e-15 Ang from
floating point coordinate conversion. The map is therefore numerically
equivalent but not byte-identical. No native SCF was launched for this check.

The preserved mathematical source and original SHA256 are in
`cpu/geometry/legacy_helpers/ORIGINS.json`.
