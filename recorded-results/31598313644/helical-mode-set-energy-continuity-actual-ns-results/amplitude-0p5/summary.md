# Actual Galerkin NS helical-mode energy continuity audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HELICAL_MODE_ENERGY_CONTINUITY__GROSS_CANONICAL_EDGE_DW_PLUS_DW_MINUS__MODE_STOCK_PLUS_VISCOSITY__SAME_CUTOFF_CROSS_FFT_REPRESENTATION__NO_REHAHN_OR_TEMPORAL_MATCHING**.

At every RK4 output time the probe reconstructs the full unordered canonical child-edge ledger from the actual evolving Fourier coefficients, restricts it to one physical child helicity, and integrates gross dW+ and dW- against the actual helical modal energy and viscous dissipation.

- common cutoff: `7`
- child/helicity: `(7, 6, 5)` / `+1`
- phase sign: `+1`
- FFT representations: `24`
- representation native energy-throughput scale: `0.250939063538`
- initial/final energy representation native residuals: `0.000e+00` / `0.000e+00`
- integrated positive/negative work representation native residuals: `0.000e+00` / `0.000e+00`
- viscous dissipation representation native residual: `0.000e+00`

## resolution 24
- steps: `48`
- initial/final helical modal energy: `0.25` / `0.249610168205`
- integrated gross positive/negative work: `0.000939063538141` / `9.23590817659e-05`
- viscous dissipation: `0.00123653625125`
- interval continuity native residual: `1.989e-12`
- worst instantaneous signed reconstruction native residual: `4.845e-16`
- global NS energy-balance residual: `1.059e-12`
- steps with positive/negative gross edge work: `49` / `49`

The positive and negative terms are the gross canonical edge Hahn marginals for this physical helical child mode, not a fresh Hahn split of its summed modal net work.  The identity is aggregate stock/flow continuity only: it creates no FIFO/LIFO matching between earlier deposits and later withdrawals and does not bound total gross nonlinear transfer variation.
