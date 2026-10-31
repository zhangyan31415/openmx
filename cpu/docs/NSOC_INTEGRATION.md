# Unified scalar and spinor CPU source integration

The NSOC sources in `cpu/overrides/nsoc` recover the tested CPU lineage ending in
`nsoc-selected-dc095f2c7b5b09cb` (25 October 2026). They retain full native SCF,
forces, checkpoint/export behavior, owner-tiled DM/EDM, shared grids/fields/VNA,
compact imaginary history, representative real history, scratch reuse and
optional D3 fused/selected eigenvector return. This integration does not change
the input SCF criterion, mixing-history length or returned-state rule.

## Source provenance

`SOURCES.json` records original source and imported-file SHA256 values. The
reconstructed scalar solver translation has exactly the SHA256 recorded by the
verified selected-return build. `selected_left.c` only changes its former
absolute include of the frozen sparse helper to a repository-local include.
The compatibility superset of `cpu_nsoc_sym.h` adds scratch API declarations to
the older header without changing existing declarations.

`recover_sources.py` is a one-time historical recovery recipe, not a build
requirement. Normal builds use the checked-in source files and stock sources
under `cpu/vendor/openmx4.0/source`; they do not consume result-tree objects,
source trees or historical binaries.

The v3 D3 and v4 C3 basis constructors are in `metadata/`, with exact origins in
`metadata/ORIGINS.json`. They now resolve the group/geometry helper from the
repository. Angular representations and accepted operation subsets remain
subject to the geometry helper's validation.

## Merge and build

Generate the deterministic merged sources/selection after changing either
archived source family:

```sh
python cpu/overrides/merge_stages.py --unified-repo "$PWD" \
  --output cpu/overrides/nsoc/UNIFIED_SELECTION.json
```

`merge_stages.py` maps each NSOC edit through the original stock-to-SOC line
alignment and applies the edits in reverse order. It rejects overlapping edits
instead of matching a repeated allocation block elsewhere in the file.
The scalar and spinor allocation/history validity predicates are kept in their
respective spin paths. All conflicts are checked before generated sources are
written.

The unified selection retains the SOC batched-force implementation, which
already contains the NSOC native shared-field force path. Scalar H-grid BLAS3
is enabled only on its original scalar eligibility conditions; spinor assembly
uses its native fallback. The common grid and VNA helpers are byte-identical between these recovered
source families. Shared-field C bodies also match, but their original private
headers encode incompatible scalar/spinor layouts. The unified build therefore
uses a new explicit common header and collectively configures two or four spin
components before allocation. Every selected consumer includes that same
header. Scalar payload retains 8 fields and SOC retains 14; configuration
disagreement is rejected without changing the live layout. See
`validation/GRID_FIELDS_FIX.json` for the failed native case and regression.

## Route contracts

| Basis | Scalar solve | Applicable extra optimizations |
|---|---|---|
| D3 v3 | Real Ex/Ey equivalent blocks plus A1/A2 | Direct-sector assembly, D3 representative reuse, fused/selected return when explicitly enabled |
| C3 v4 | Conjugate complex pair plus real invariant sector | Common owner DM/EDM and sparse transforms; D3-only switches must remain disabled |
| Generic v2 | Validated real isotypic sectors | Common scalar sector solver; no claim of D3 equivalent-block reuse |
| Off | Native scalar path | Only independently compatible memory/runtime options |

The experimental scalar solver supports fixed-geometry, nonmagnetic Gamma
Cluster calculations with an uncontracted basis. Enabling `FULL_PROJECT` is an
explicit physical symmetry projection, distinct from implementation-equivalent
fused/selected return. Unsupported combinations fail rather than silently
falling back to unoptimized NSOC.

Representative real-history compression retains the requested History30 and
remains restricted to that verified value. The separately recovered 28 October
H-history capacity header accepts any positive input history while compacting
only unused Hamiltonian history to two slots; it does not shorten DM history.
The original History30-only header remains unchanged in the archival closure.
For History35, representative real-history compression must remain disabled. SOC and NSOC invocations must use their
own generated environments, since scalar-only flags intentionally reject SOC.

## Optional XC workspace reuse

The unified build includes the verified `Set_XC_Grid.workspace.c` patch and
`cpu_xc_workspace.h`, both recovered from the 28 October generic-memory work.
The generated source matches that archived source exactly; its helper header
SHA256 matches the previous native build. Runtime `OPENMX_CPU_XC_WORKSPACE=0`
remains the default. Mode 1 reuses private gradient storage. Mode 2 additionally
aliases exactly equal finite scalar spin-density inputs/core corrections;
unequal inputs are rejected. External derivative outputs keep native storage.
The XC arithmetic loops are unchanged.

Previous component validation covered 144 cases with zero grid-value difference
in accepted cases. The earlier HgI2 7.63° native run completed 32 NSOC steps with
History35, StartPulay15, H-history capacity reuse and XC mode2. This is prior
lineage evidence, not a substitute for validating the newly unified binary.
`MEMORY_ORIGINS.json` records the exact source and run provenance.

## Validation status

The merger/generator unit checks and compile-only checks validate source
integration, not timing or scientific convergence. Native NSOC and one-step SOC
benchmarks on node040–043 are required before this candidate is committed.
