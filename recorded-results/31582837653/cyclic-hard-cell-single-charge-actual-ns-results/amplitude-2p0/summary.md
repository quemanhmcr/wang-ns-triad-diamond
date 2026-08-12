# Actual Galerkin NS cyclic hard-cell single-charge audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_CYCLIC_HARD_CELL_SINGLE_CHARGE_AUDIT__ACTUAL_DW_MINUS_ROWS_TO_CANONICAL_DW_PLUS_FATE_COLUMNS__COARSE_SELF_LOOPS_ZERO_DEPTH__ONE_DONOR_TWO_RECIPIENT_AND_PHASE_REVERSED_TWO_DONOR_ONE_RECIPIENT**.

The probe evolves the repository's real 2/3-dealiased incompressible Fourier--Galerkin Navier--Stokes system. At each physical snapshot it reads the actual three-root signed helical triad work, forms the certified cyclic donor kernel, and only then pushes donor/recipient roots through deterministic hard-cell maps.

- common cutoff: `7`
- phase sign: `+1`
- FFT representations: `24`
- maximum total-work representation native residual: `0.000e+00`
- maximum good-recipient representation native residual: `0.000e+00`
- maximum bad-recipient representation native residual: `0.000e+00`

## resolution 24
- steps/snapshots: `24` / `4`
- worst row/column native residuals: `5.735e-17` / `3.829e-17`
- worst fate-partition native residual: `0.000e+00`
- worst restricted-donor pushforward native residual: `5.735e-17`
- worst fine/coarse total/fate native residuals: `0.000e+00` / `0.000e+00`
- NS energy-balance residual: `2.906e-12`
- initial donor/recipient counts: `1` / `2`
- initial good/bad recipient masses: `68.0888945495` / `21.7191397006`
- initial overlapping recipient charges: `0`
- initial maximal-coarsening self-loop fraction: `1`

Changing the sign of the whole divergence-free initial field is used only to expose the opposite physical cubic-work sign pattern at t=0: the + branch begins one-donor/two-recipient, while the - branch begins two-donor/one-recipient. Both are then evolved by the same Navier--Stokes equations.

Coarse self-loops remain real same-time redistribution. They create no event depth and no scale progress. Geometry-bad recipient work remains the existing transfer-loss causal sublaw, not vanished PDE energy. No between-time deposit/withdrawal matching and no global-regularity claim are made.
