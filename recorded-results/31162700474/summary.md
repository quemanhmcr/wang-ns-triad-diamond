# Localized SGS / pressure-work ledger

The identities are analytic; random traces only validate implementation.

- random finite chains: `50000`
- worst weighted summation-by-parts residual: `1.421e-14`
- worst positive-work depletion margin: `7.414e-03`
- pressure-cancellation branches checked: `27037`

For a positive raw SGS flux `S`, either combined work retains at least `S/2`,
or pressure boundary work has magnitude at least `S/2`. In the latter branch,

`int_A (|u|^3+|p|^(3/2)) >= S/(2 ||grad chi||_inf)`.
