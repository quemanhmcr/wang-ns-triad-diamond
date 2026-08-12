# Actual Galerkin NS helical-mode energy continuity audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HELICAL_MODE_ENERGY_CONTINUITY__GROSS_CANONICAL_EDGE_DW_PLUS_DW_MINUS__MODE_STOCK_PLUS_VISCOSITY__SAME_CUTOFF_CROSS_FFT_REPRESENTATION__NO_REHAHN_OR_TEMPORAL_MATCHING**.

At every RK4 output time the probe reconstructs the full unordered canonical child-edge ledger from the actual evolving Fourier coefficients, restricts it to one physical child helicity, and integrates gross dW+ and dW- against the actual helical modal energy and viscous dissipation.

- common cutoff: `7`
- child/helicity: `(7, 6, 5)` / `-1`
- phase sign: `+1`
- FFT representations: `24, 28`
- representation native energy-throughput scale: `7.44536678347e-05`
- initial/final energy representation native residuals: `6.231e-28` / `6.444e-14`
- integrated positive/negative work representation native residuals: `2.239e-14` / `1.209e-16`
- viscous dissipation representation native residual: `9.386e-17`

## resolution 24
- steps: `96`
- initial/final helical modal energy: `6.9118235465e-32` / `7.38660111133e-05`
- integrated gross positive/negative work: `7.44536161723e-05` / `3.92012387067e-07`
- viscous dissipation: `1.9564433437e-07`
- interval continuity native residual: `6.939e-07`
- worst instantaneous signed reconstruction native residual: `2.128e-15`
- global NS energy-balance residual: `4.400e-13`
- steps with positive/negative gross edge work: `97` / `97`

## resolution 28
- steps: `96`
- initial/final helical modal energy: `2.27291407636e-32` / `7.38660111133e-05`
- integrated gross positive/negative work: `7.44536161723e-05` / `3.92012387067e-07`
- viscous dissipation: `1.9564433437e-07`
- interval continuity native residual: `6.939e-07`
- worst instantaneous signed reconstruction native residual: `7.739e-15`
- global NS energy-balance residual: `4.392e-13`
- steps with positive/negative gross edge work: `97` / `97`

The positive and negative terms are the gross canonical edge Hahn marginals for this physical helical child mode, not a fresh Hahn split of its summed modal net work.  The identity is aggregate stock/flow continuity only: it creates no FIFO/LIFO matching between earlier deposits and later withdrawals and does not bound total gross nonlinear transfer variation.
