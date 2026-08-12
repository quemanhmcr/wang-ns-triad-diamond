# Actual Galerkin NS helical-mode energy continuity audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HELICAL_MODE_ENERGY_CONTINUITY__GROSS_CANONICAL_EDGE_DW_PLUS_DW_MINUS__MODE_STOCK_PLUS_VISCOSITY__SAME_CUTOFF_CROSS_FFT_REPRESENTATION__NO_REHAHN_OR_TEMPORAL_MATCHING**.

At every RK4 output time the probe reconstructs the full unordered canonical child-edge ledger from the actual evolving Fourier coefficients, restricts it to one physical child helicity, and integrates gross dW+ and dW- against the actual helical modal energy and viscous dissipation.

- common cutoff: `7`
- child/helicity: `(7, 6, 5)` / `-1`
- phase sign: `+1`
- FFT representations: `24, 28`
- representation native energy-throughput scale: `1.81542776069e-05`
- initial/final energy representation native residuals: `1.478e-28` / `3.079e-14`
- integrated positive/negative work representation native residuals: `4.890e-14` / `2.595e-16`
- viscous dissipation representation native residual: `1.910e-16`

## resolution 24
- steps: `64`
- initial/final helical modal energy: `1.28590080344e-32` / `1.80191453799e-05`
- integrated gross positive/negative work: `1.81542566415e-05` / `9.53908587922e-08`
- viscous dissipation: `3.97413682063e-08`
- interval continuity native residual: `1.155e-06`
- worst instantaneous signed reconstruction native residual: `1.890e-15`
- global NS energy-balance residual: `7.410e-13`
- steps with positive/negative gross edge work: `65` / `65`

## resolution 28
- steps: `64`
- initial/final helical modal energy: `1.55429512854e-32` / `1.80191453799e-05`
- integrated gross positive/negative work: `1.81542566415e-05` / `9.53908587923e-08`
- viscous dissipation: `3.97413682063e-08`
- interval continuity native residual: `1.155e-06`
- worst instantaneous signed reconstruction native residual: `1.461e-15`
- global NS energy-balance residual: `7.409e-13`
- steps with positive/negative gross edge work: `65` / `65`

The positive and negative terms are the gross canonical edge Hahn marginals for this physical helical child mode, not a fresh Hahn split of its summed modal net work.  The identity is aggregate stock/flow continuity only: it creates no FIFO/LIFO matching between earlier deposits and later withdrawals and does not bound total gross nonlinear transfer variation.
