# Actual Fourier--Galerkin Navier--Stokes referee: material ownership is downstream of service

Status: **ACTUAL_DEALIASED_FOURIER_GALERKIN_NS_INCREMENT_SERVICE__SAME_PHYSICAL_LAW_DIFFERENT_MATERIAL_READINGS__MATERIAL_INTERFACE_MASS_IS_NOT_GENERATION**.

The referee evolves the same dealiased incompressible Galerkin Navier--Stokes state used by the physical helical audits. At genuine trajectory times it forms one actual positive velocity-increment square law at the fixed physical displacement r=(pi/2,0,0). The state and service law are then held literally fixed while two geometric old-material readings are applied: an x-slab, whose endpoints cross the displacement boundary, and a y-slab, whose endpoints do not. The OO/ON/NN masses change, but the underlying service and Navier--Stokes trajectory do not.

- maximum cross-FFT representation service residual: `2.065e-16`
- minimum OO/ON/NN partition-change fraction at fixed service: `1`

## resolution 24
- global NS energy-balance residual: `2.963e-12`
- minimum positive x-boundary interface service: `8.59536445504`
- maximum material-rereading service residual: `0.000e+00`
- maximum global nonlinear-work relative rate: `1.111e-18`
- maximum divergence / initial L2: `4.498e-16`
- t=0: total service `17.24`, x-ON `8.62`, y-ON `0`, partition L1 change `17.24`
- t=0.0003125: total service `17.2246252688`, x-ON `8.61231265863`, y-ON `0`, partition L1 change `17.2246253173`
- t=0.00065625: total service `17.2076892231`, x-ON `8.60384471823`, y-ON `0`, partition L1 change `17.2076894365`
- t=0.001: total service `17.1907284158`, x-ON `8.59536445504`, y-ON `0`, partition L1 change `17.1907289101`

## resolution 28
- global NS energy-balance residual: `2.963e-12`
- minimum positive x-boundary interface service: `8.59536445973`
- maximum material-rereading service residual: `0.000e+00`
- maximum global nonlinear-work relative rate: `1.174e-18`
- maximum divergence / initial L2: `4.664e-16`
- t=0: total service `17.24`, x-ON `8.62`, y-ON `0`, partition L1 change `17.24`
- t=0.0003125: total service `17.2246252688`, x-ON `8.61231265909`, y-ON `0`, partition L1 change `17.2246253182`
- t=0.00065625: total service `17.2076892231`, x-ON `8.60384472025`, y-ON `0`, partition L1 change `17.2076894405`
- t=0.001: total service `17.1907284158`, x-ON `8.59536445973`, y-ON `0`, partition L1 change `17.1907289195`

This referee does not invent a material PDE, does not infer a causal supplier for the snapshot increment law, and does not claim that the chosen geometric slabs are persistent material packets. Its point is narrower and exact: for an actual NS field, changing only the material reading can change interface/fresh provenance while the positive physical service law is unchanged. Therefore material ownership alone cannot be its generator. No global-regularity claim is made.
