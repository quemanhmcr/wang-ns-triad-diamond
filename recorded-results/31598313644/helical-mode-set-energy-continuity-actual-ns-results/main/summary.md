# Actual Galerkin NS helical-mode energy continuity audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HELICAL_MODE_ENERGY_CONTINUITY__GROSS_CANONICAL_EDGE_DW_PLUS_DW_MINUS__MODE_STOCK_PLUS_VISCOSITY__SAME_CUTOFF_CROSS_FFT_REPRESENTATION__NO_REHAHN_OR_TEMPORAL_MATCHING**.

At every RK4 output time the probe reconstructs the full unordered canonical child-edge ledger from the actual evolving Fourier coefficients, restricts it to one physical child helicity, and integrates gross dW+ and dW- against the actual helical modal energy and viscous dissipation.

- common cutoff: `7`
- child/helicity: `(7, 6, 5)` / `+1`
- phase sign: `+1`
- FFT representations: `24, 28`
- representation native energy-throughput scale: `1.01000884888`
- initial/final energy representation native residuals: `0.000e+00` / `8.794e-16`
- integrated positive/negative work representation native residuals: `1.718e-18` / `0.000e+00`
- viscous dissipation representation native residual: `0.000e+00`

## resolution 24
- steps: `64`
- initial/final helical modal energy: `1` / `1.00240956249`
- integrated gross positive/negative work: `0.0100088488666` / `0.000991283875857`
- viscous dissipation: `0.00660800252052`
- interval continuity native residual: `1.563e-11`
- worst instantaneous signed reconstruction native residual: `8.081e-16`
- global NS energy-balance residual: `7.410e-13`
- steps with positive/negative gross edge work: `65` / `65`

## resolution 28
- steps: `64`
- initial/final helical modal energy: `1` / `1.00240956249`
- integrated gross positive/negative work: `0.0100088488666` / `0.000991283875857`
- viscous dissipation: `0.00660800252052`
- interval continuity native residual: `1.563e-11`
- worst instantaneous signed reconstruction native residual: `6.471e-16`
- global NS energy-balance residual: `7.409e-13`
- steps with positive/negative gross edge work: `65` / `65`

The positive and negative terms are the gross canonical edge Hahn marginals for this physical helical child mode, not a fresh Hahn split of its summed modal net work.  The identity is aggregate stock/flow continuity only: it creates no FIFO/LIFO matching between earlier deposits and later withdrawals and does not bound total gross nonlinear transfer variation.
