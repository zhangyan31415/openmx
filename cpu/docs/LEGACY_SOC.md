# Historical C3/D3 SOC implementation

This document records algorithm and validation provenance. Normal symmetry use
is documented in the [input guide](INPUT_SYMMETRY.zh-CN.md): add
`scf.Symmetry on` and run OpenMX normally.

The preserved C3/D3 SOC implementation consumes a sparse spinor basis and a
matching native real-space map. It provides an earlier specialized path beside
the general spinor solver. It does not truncate states or modify coordinates.
The mathematical source and original hashes are recorded in
`cpu/geometry/legacy_helpers/ORIGINS.json`.

## Algorithms and output

The historical D3 optimization combines equivalent representation blocks,
time-reversal and real reductions, owner-sparse density matrices, compact
coefficients, local H/S and overlap storage, shared force density matrices,
BLAS-3 density construction and batched nonlocal-force contraction. C3 does
not inherit the four-sector D3 time-reversal implementation; the general
spinor framework implements its separate C3 route.

The fastest saved D3 configuration projected native spinor H/S before export.
Earlier variants also supported original assembled H/S export even when the
solver worked with projected sectors; those exports must not be interpreted
as exactly covariant operators. Current simple input behavior and operator
consistency checks are described in the input guide. Geometry and basis
preparation checks do not replace H/S covariance or physical-error checks.

The complete saved D3 algorithm configuration originates from the verified
2026-10-25 `select_on1_soc` run, binary SHA256
`9ba1dd686da8abef7d45a76e19dcd91aa49ac6e52a4a18f9e50004e8a4f59acc`.
The historical record retains the original configuration and source hashes.
It is provenance for a particular measured path, not a guarantee that a
current automatically selected path has the same timing.

## Scope

The preserved specialized path was validated for fixed geometry,
nonmagnetic Gamma NSOC charge followed by one c2n Cluster SOC step, in the
established hexagonal slab orientation. Other orientations require compatible
general representation handling. Runtime guards retain restrictions on
contraction, Hubbard terms, fields and auxiliary states.

Integer atom/image actions and group closure are revalidated against the SOC
geometry. The historical rotation convention and source C2 gauge are retained;
nominal-rotation residuals are checked independently. A valid geometry does
not by itself certify the electronic operator or final physical accuracy.

## Saved validation

Metadata tests covered C3 cycle/fixed atoms, D3 partner and spinor gauges,
malformed permutations, binary headers, export semantics, one-step boundaries
and the preserved D3 configuration. Common integration tests checked artifact
selection and rejection of incompatible routes. These checks passed; they are
metadata/integration tests, not native SCF performance measurements.

`LEGACY_SOC_MATCHED_GAUGE.json` records a comparison with the exact MoTe₂
7.34° / 366-atom basis used by the verified 2026-10-25 fused-return runs.
The regenerated basis is byte-identical, SHA256
`f605c9ff3a8ee254d9680bdf81c7725ad636a20433b12dacdc701e7cbca7468e`.
The matching SOC map has identical integer actions and spin/orbital rotation
blocks; Cartesian coordinates differ by at most 7.11e-15 Å due to coordinate
conversion. The maps are numerically equivalent but not byte-identical.
No native SCF was launched for this reconstruction check.

Physical-accuracy and performance limitations remain in the
[benchmark report](BENCHMARKS_20261031.md). In particular, historical
projected-versus-original force failures must be retained when interpreting
implementation/layout agreement.
