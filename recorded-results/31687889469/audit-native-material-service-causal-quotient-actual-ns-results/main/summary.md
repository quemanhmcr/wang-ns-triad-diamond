# Actual Fourier--Galerkin Navier--Stokes referee: material ownership is downstream of service

Status: **ACTUAL_DEALIASED_FOURIER_GALERKIN_NS_INCREMENT_SERVICE__SAME_PHYSICAL_LAW_DIFFERENT_MATERIAL_READINGS__MATERIAL_INTERFACE_MASS_IS_NOT_GENERATION**.

The referee evolves the same dealiased incompressible Galerkin Navier--Stokes state used by the physical helical audits. At genuine trajectory times it forms one actual positive velocity-increment square law at the fixed physical displacement r=(pi/2,0,0). The state and service law are then held literally fixed while two geometric old-material readings are applied: an x-slab, whose endpoints cross the displacement boundary, and a y-slab, whose endpoints do not. The OO/ON/NN masses change, but the underlying service and Navier--Stokes trajectory do not.

- maximum cross-FFT representation service residual: `2.074e-16`
- minimum OO/ON/NN partition-change fraction at fixed service: `1`

## resolution 24
- global NS energy-balance residual: `2.614e-11`
- minimum positive x-boundary interface service: `2.1415665097`
- maximum material-rereading service residual: `0.000e+00`
- maximum global nonlinear-work relative rate: `6.941e-19`
- maximum divergence / initial L2: `3.807e-16`
- t=0: total service `4.31`, x-ON `2.155`, y-ON `0`, partition L1 change `4.31`
- t=0.0005: total service `4.30101941243`, x-ON `2.15050971009`, y-ON `0`, partition L1 change `4.30101942018`
- t=0.001: total service `4.29206376401`, x-ON `2.14603189743`, y-ON `0`, partition L1 change `4.29206379486`
- t=0.0015: total service `4.28313295029`, x-ON `2.1415665097`, y-ON `0`, partition L1 change `4.2831330194`

## resolution 28
- global NS energy-balance residual: `2.614e-11`
- minimum positive x-boundary interface service: `2.14156651036`
- maximum material-rereading service residual: `0.000e+00`
- maximum global nonlinear-work relative rate: `9.910e-19`
- maximum divergence / initial L2: `4.477e-16`
- t=0: total service `4.31`, x-ON `2.155`, y-ON `0`, partition L1 change `4.31`
- t=0.0005: total service `4.30101941243`, x-ON `2.15050971016`, y-ON `0`, partition L1 change `4.30101942032`
- t=0.001: total service `4.29206376401`, x-ON `2.14603189772`, y-ON `0`, partition L1 change `4.29206379544`
- t=0.0015: total service `4.28313295029`, x-ON `2.14156651036`, y-ON `0`, partition L1 change `4.28313302071`

This referee does not invent a material PDE, does not infer a causal supplier for the snapshot increment law, and does not claim that the chosen geometric slabs are persistent material packets. Its point is narrower and exact: for an actual NS field, changing only the material reading can change interface/fresh provenance while the positive physical service law is unchanged. Therefore material ownership alone cannot be its generator. No global-regularity claim is made.
