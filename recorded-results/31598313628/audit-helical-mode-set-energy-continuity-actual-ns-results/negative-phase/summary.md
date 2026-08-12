# Actual Galerkin NS helical-mode energy continuity audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HELICAL_MODE_ENERGY_CONTINUITY__GROSS_CANONICAL_EDGE_DW_PLUS_DW_MINUS__MODE_STOCK_PLUS_VISCOSITY__SAME_CUTOFF_CROSS_FFT_REPRESENTATION__NO_REHAHN_OR_TEMPORAL_MATCHING**.

At every RK4 output time the probe reconstructs the full unordered canonical child-edge ledger from the actual evolving Fourier coefficients, restricts it to one physical child helicity, and integrates gross dW+ and dW- against the actual helical modal energy and viscous dissipation.

- common cutoff: `7`
- child/helicity: `(7, 6, 5)` / `+1`
- phase sign: `-1`
- FFT representations: `24, 28`
- representation native energy-throughput scale: `1.69255786931`
- initial/final energy representation native residuals: `2.624e-16` / `3.936e-16`
- integrated positive/negative work representation native residuals: `7.687e-19` / `4.100e-18`
- viscous dissipation representation native residual: `4.100e-18`

## resolution 24
- steps: `96`
- initial/final helical modal energy: `1.69` / `1.6529287465`
- integrated gross positive/negative work: `0.00255786931112` / `0.0263913748228`
- viscous dissipation: `0.0132377478819`
- interval continuity native residual: `6.485e-11`
- worst instantaneous signed reconstruction native residual: `8.834e-16`
- global NS energy-balance residual: `6.435e-12`
- steps with positive/negative gross edge work: `97` / `97`

## resolution 28
- steps: `96`
- initial/final helical modal energy: `1.69` / `1.6529287465`
- integrated gross positive/negative work: `0.00255786931112` / `0.0263913748228`
- viscous dissipation: `0.0132377478819`
- interval continuity native residual: `6.485e-11`
- worst instantaneous signed reconstruction native residual: `1.033e-15`
- global NS energy-balance residual: `6.434e-12`
- steps with positive/negative gross edge work: `97` / `97`

The positive and negative terms are the gross canonical edge Hahn marginals for this physical helical child mode, not a fresh Hahn split of its summed modal net work.  The identity is aggregate stock/flow continuity only: it creates no FIFO/LIFO matching between earlier deposits and later withdrawals and does not bound total gross nonlinear transfer variation.
