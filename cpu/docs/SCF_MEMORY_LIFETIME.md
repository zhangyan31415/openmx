# Scalar SCF memory lifetime

The fixed-geometry CPU workflow now has two exact memory-lifetime changes.
They do not change orbital basis, integration grid, returned-state policy,
SCF tolerance, arithmetic precision, Hamiltonian export or native force output.

## Packed overlap factors

After symmetry-sector Cholesky factorization, the overlap allocation contains
only the live packed factors in its prefix. Complete pages beyond that prefix
are returned to Linux using MADV_DONTNEED. The allocation address and virtual
capacity remain valid, allowing a subsequent explicit raw-S rebuild to write
those pages again. Boundary pages are preserved to protect both the live prefix
and allocator/neighbor data. Unsupported platforms or failed advice continue
without reclamation. No memory saving is inferred from the advised byte count;
resident/PSS/cgroup measurements remain separate.

## VNA value/derivative stages

The VNA shared cache keeps its four-component FP32 layout. In the selected
memory-saving path, the SCF phase computes/stores component0 only. Gradient
components1–3 are generated with the original quadrature and real-harmonic
conversion immediately before Force. Raw inter-node publication transfers only
component0; weighted force publication transfers the three gradient components.
A collective readiness check forbids weighting incomplete derivatives.

Before derivative generation the now-dead Hs/Ss/Cs/EVec1 allocations are freed.
Sparse H/S tensors, density and energy-density matrices, occupations, eigenvalue
lists and communicator infrastructure remain available to native energy/force
and H/S output. End-of-DFT cleanup is safe with the cleared pointers.

The path is selected automatically only for fixed, first-MD-step, scalar NSOC
Cluster calculations already eligible for SCF-buffer release, with VNA sharing
active and without contraction/LNO/partial-charge/temperature-control or
state-based postprocessing/stress requirements. The selection additionally
requires **every physical node's conservative dead-dense capacity** to cover its
deferred derivative payload. The dense count uses actual eigen-scratch capacity
and allocated Hs/Cs sizes; S is excluded because most of its pages were already
reclaimed. This is a memory-benefit selection rule, not an OOM guarantee.

Small cases retain the joint native value+gradient construction. Other solver,
spin and postprocessing paths retain their previous VNA lifetime. Scalar and
spinor builds share the same projector/cache implementation; no material,
angle or point-group special case is used.

The two-stage construction repeats some common intermediate integral work.
The final derivative values are generated once; they are not approximated or
recovered from rounded matrix outputs. Additional work occurs once per fixed
calculation, rather than once per SCF iteration. It can increase runtime; consult
matched measurements instead of assuming a zero-cost memory reduction.

Diagnostics:

- `CPU_NSOC_OVERLAP_RECLAIM`: allocation/live-factor/advised-tail byte counts;
- `CPU_VNA_LIFETIME_POLICY`: recorded dense capacity, derivative payload and
  collective selection;
- `CPU_VNA_DEFERRED`: value-only SCF publication contract;
- `CPU_NSOC_DENSE_RELEASED`: dead solver arrays freed before the gradient pass;
- `CPU_VNA_DEFERRED_COMPLETE`: native gradient-pass time.

The MPI state test covers both deferred and legacy publications, reinitialization,
physical-node communication and rejection of premature weighting. Native
acceptance must compare energy, complete periodic H/S and forces: this change
alters when derivative arrays are constructed, so an implementation-equivalence
force check is required even when forces are not a production timing gate.
