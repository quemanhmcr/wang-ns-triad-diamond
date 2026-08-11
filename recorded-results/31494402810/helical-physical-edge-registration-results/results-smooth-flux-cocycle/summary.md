# Smooth log-scale flux cocycle

The theorem identities are analytic; the stress test only validates implementation.

- gamma*: `0.492815285342135`
- filter-independent all-log Mellin moment: exact for every conservative finite interaction
- common-midgap smooth-tail identity: exact on a spectral moat
- no-cancellation polarization identity: exact
- stress blocks: `50000`
- test filter half-width delta: `0.050` log units
- test shell half-width: `0.060` log units
- minimum common-moat margin seen: `0.087018431`
- worst midgap equality residual: `1.421e-14`
- worst polarization residual: `6.661e-16`

The key PDE-facing formula is

`2 int_tau^infty Pi_delta(t) dt = sum_e T_e log(q_e/p_e)`

when tau is the positive-transfer-weighted mean of the edge midgaps and a common
smooth transition moat separates all top parents from all children.
