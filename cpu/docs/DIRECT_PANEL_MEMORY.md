# Direct D3 assembly memory

The scalar Gamma D3 direct assembler computes representative atom-block
contributions to the symmetry sectors. It reuses the existing eigenvector
scratch for its sector output. Skipping full-AO assembly does **not** remove
all solver matrices: the current solver still allocates distributed Hs, Ss,
Cs and the eigenvector storage needed by its later stages.

## Bounded communication

The assembler plans contiguous batches of local panels. It exchanges proposed
counts and reduces a batch until both its sends and the receiving ranks' inputs
fit their budgets. Empty source ranks still participate. If many one-panel
senders overload one receiver, that batch uses one sender at a time. A single
panel that cannot fit is a collective resource failure.

Each panel's BLAS contractions execute once per H/S application. After each
exchange the received values are accumulated into the existing sector output;
the same send/receive arrays are reused for the next batch. Small messages use
one batch. Batch count and actual capacities are printed in
`CPU_NSOC_DIRECT_PLAN`. Floating-point summation order can differ between
batches; the generalized eigenproblem and symmetry operations are unchanged.

The default combined send/receive capacity is at most 8 MiB per MPI rank.
This is a cap, not a preallocation. The planner further limits it by its available
allocation budget. Basis metadata, routes, the packed atom-block input, solver
matrices, grids and mixing history are **not** included in this 8 MiB figure.
All routing metadata is currently retained. Metadata construction can therefore
still exhaust the overall planner budget; the change does not promise constant
memory for the entire solver or eliminate all large-count limits.

The total planner budget uses one eighth of measured node/cgroup headroom,
divided among local ranks. Resource/count failure can collectively fall back to
the established projected-D3 assembly and paired solver. Invalid basis or field
metadata remains an error. The fallback keeps compact Work1 and owner DM/EDM
configuration; it is not a switch back to stock OpenMX.

## Developer diagnostics

`OPENMX_CPU_NSOC_PANEL_BUFFER_BYTES` overrides the combined numerical buffer
cap. `OPENMX_CPU_NSOC_PANEL_BUDGET_BYTES` can lower the whole planner budget.
These two resource controls survive native symmetry preparation; inherited
solver-selector flags continue to be cleared. No extra public input is required.
A cap of zero is useful only for testing resource failure/fallback.

Tests must inspect the actual branch. A requested environment variable alone is
not proof: require `value_batches > 1` for streaming, verify the printed capacity
bound, and require `CPU_NSOC_DIRECT_FALLBACK` for fallback tests. Earlier test
records whose variables were cleared by native bootstrap cannot establish those
branches were exercised.

`cpu/overrides/nsoc/tests/panel_streaming.c` compares repeated H/S applications
against an independently constructed dense D3-sector oracle, with uneven input
ownership, empty ranks, and reversed BLACS rank mappings. It takes buffer bytes
and grid row count as optional arguments. Compile it with `rep_panel.c`, MPI and
MKL, and run it on a compute node with 1, 2, 4, 8 and 64 MPI ranks. The default
512-byte cap deliberately creates many batches; it is not a performance setting.
