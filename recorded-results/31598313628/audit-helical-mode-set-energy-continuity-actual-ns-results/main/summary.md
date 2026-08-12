# Actual Galerkin NS helical-mode energy continuity audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HELICAL_MODE_ENERGY_CONTINUITY__GROSS_CANONICAL_EDGE_DW_PLUS_DW_MINUS__MODE_STOCK_PLUS_VISCOSITY__SAME_CUTOFF_CROSS_FFT_REPRESENTATION__NO_REHAHN_OR_TEMPORAL_MATCHING**.

At every RK4 output time the probe reconstructs the full unordered canonical child-edge ledger from the actual evolving Fourier coefficients, restricts it to one physical child helicity, and integrates gross dW+ and dW- against the actual helical modal energy and viscous dissipation.

- common cutoff: `7`
- child/helicity: `(7, 6, 5)` / `+1`
- phase sign: `+1`
- FFT representations: `24, 28`
- representation native energy-throughput scale: `1.71637024221`
- initial/final energy representation native residuals: `2.587e-16` / `1.035e-15`
- integrated positive/negative work representation native residuals: `0.000e+00` / `1.263e-18`
- viscous dissipation representation native residual: `2.021e-18`

## resolution 24
- steps: `96`
- initial/final helical modal energy: `1.69` / `1.70031618647`
- integrated gross positive/negative work: `0.0263702421445` / `0.00262817336548`
- viscous dissipation: `0.0134258823783`
- interval continuity native residual: `4.051e-11`
- worst instantaneous signed reconstruction native residual: `7.367e-16`
- global NS energy-balance residual: `4.400e-13`
- steps with positive/negative gross edge work: `97` / `97`

## resolution 28
- steps: `96`
- initial/final helical modal energy: `1.69` / `1.70031618647`
- integrated gross positive/negative work: `0.0263702421445` / `0.00262817336548`
- viscous dissipation: `0.0134258823783`
- interval continuity native residual: `4.051e-11`
- worst instantaneous signed reconstruction native residual: `8.825e-16`
- global NS energy-balance residual: `4.392e-13`
- steps with positive/negative gross edge work: `97` / `97`

The positive and negative terms are the gross canonical edge Hahn marginals for this physical helical child mode, not a fresh Hahn split of its summed modal net work.  The identity is aggregate stock/flow continuity only: it creates no FIFO/LIFO matching between earlier deposits and later withdrawals and does not bound total gross nonlinear transfer variation.
