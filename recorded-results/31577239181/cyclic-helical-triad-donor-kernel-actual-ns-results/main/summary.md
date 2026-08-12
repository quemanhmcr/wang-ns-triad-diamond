# Actual Galerkin NS cyclic helical-triad donor-kernel audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_CYCLIC_HELICAL_TRIAD_DONOR_KERNEL_AUDIT__ACTUAL_MODAL_COEFFICIENTS__THREE_ROOT_ENERGY_CONSERVATION__DW_MINUS_TO_DW_PLUS_MARGINALS__SIGNED_GOOD_NONFORWARD_SIDE_RECIPIENT**.

The probe evolves the repository's real 2/3-dealiased incompressible Fourier--Galerkin Navier--Stokes system.  At each physical snapshot it reads the actual evolving helical coefficients of one closed triad and registers all three cyclic child-work roots before Hahn.

- common cutoff: `7`
- FFT representations: `24, 28`
- maximum selected-triad positive-work representation native residual: `5.751e-17`
- maximum selected root-work representation native residual: `7.668e-17`

## resolution 24
- steps/snapshots: `32` / `5`
- worst cyclic energy-conservation native residual: `8.148e-17`
- worst donor/recipient measure marginal native residuals: `5.751e-17` / `1.917e-17`
- NS energy-balance residual: `2.963e-12`
- initial signed-good efficiency: `0.999996646308`
- initial side/child work ratio: `0.318982116604`
- initial child/donor work ratio: `0.758160393088`
- initial side/donor work ratio: `0.241839606912`
- initial side forward ratio/J: `0.610513941468` / `0`
- initial donor/recipient counts: `1` / `2`

## resolution 28
- steps/snapshots: `32` / `5`
- worst cyclic energy-conservation native residual: `7.179e-17`
- worst donor/recipient measure marginal native residuals: `3.834e-17` / `3.829e-17`
- NS energy-balance residual: `2.963e-12`
- initial signed-good efficiency: `0.999996646308`
- initial side/child work ratio: `0.318982116604`
- initial child/donor work ratio: `0.758160393088`
- initial side/donor work ratio: `0.241839606912`
- initial side forward ratio/J: `0.610513941468` / `0`
- initial donor/recipient counts: `1` / `2`

The side recipient is actual positive nonforward NS work at the same triad/time.  The probe does not reinterpret it as dissipation or a reset.  The donor kernel adds provenance to canonical dW+ and creates no new event.  No global-regularity claim is made.
