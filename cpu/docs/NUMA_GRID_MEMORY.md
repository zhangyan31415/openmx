# Orbital-grid NUMA replicas

`scf.Symmetry on` prepares the existing fixed-geometry CPU cache profile.
NUMA replicas are a locality optimization: canonical orbital values are still
stored once per physical node, while remote-domain values can be read from a
local replica. No changed basis, grid, symmetry tolerance or numerical precision
is involved in the replica optimization.

The previous implementation copied every other NUMA domain's atoms into every
domain. The current implementation forms the union of first-neighbor atoms
needed by the centers owned by each physical-node/NUMA domain, using the existing
FNAN/natn graph. Only needed atoms whose owners belong to another local NUMA
domain are replicated. Off-node ghosts retain the existing native evaluation.

Replica payload is packed with per-atom offsets. There are no canonical-layout
holes for owner atoms or unused remote atoms. Row reads for an uncached atom
return its canonical shared row, so a future reader outside the planned graph
remains numerically correct. Publishing copies only the planned atoms, once
owner writes are complete; the same behavior applies to both orbital-grid
publication sites, including the derivative/force preparation site.

`CPU_GRID_NUMA` reports:

- `canonical_bytes`: the canonical node-shared orbital grid capacity;
- `replica_copy_payload_bytes`: bytes actually selected for this domain;
- `replica_mapped_capacity_bytes`: the compact replica capacity;
- `all_remote_copy_bytes`: the old policy's remote payload for this domain;
- `replica_cached_atoms` and `policy=first_neighbor_demand`.

Summing mapped capacities describes allocation capacity, not a synchronous
physical-resident peak. In particular, the old full-size mappings contained
untouched owner-domain holes; their removal must not be claimed as physical
memory saved. Use copy payload and a measured task PSS/cgroup peak separately.

The implementation lives in `cpu/overrides/soc/source/cpu_grid_shared.c`.
The NSOC helper source/header forward to the same implementation/API. The
standalone unified build already selects the common helper for both scalar and
spinor calculations; the duplicate scalar implementation has been removed.

The MPI test `cpu/overrides/nsoc/tests/grid_numa_demand.c` exercises sparse,
self-only and fully connected neighbor graphs, canonical fallback, exact values
and pointer routing, two publication epochs, destruction/reinitialization, and
physical-node identity with repeated NUMA IDs. It requires ranks pinned across
two physical NUMA domains per host. Every test reader finishes before the next
publication's owner writes begin.

This change does not compress VNA projector/derivative caches, mixing history,
or dense solver matrices. Larger savings from those objects require separate
lifetime/representation changes and validation.
