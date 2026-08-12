# Actual Galerkin NS hard-tail true upward supply audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HARD_TAIL_TRUE_UPWARD_SUPPLY__FULL_TAIL_STOCK_WORK_VISCOSITY_SEPARATE_FROM_SELECTED_CYCLIC_SUPPORT__PURE_UV_FIRST_SHELL_AND_DEEP_RESOLVED_CONTACT__NO_INTERFACE_OVERCLAIM**.

The full-tail stock/work/viscosity trajectory remains the certified radial PDE reading.  Selected cyclic triads are separate sub-observables used only to test parent-support geometry; they are never substituted for the full tail law.

- main FFT representations: `24, 28`
- pure-support cross-FFT native residual: `2.087e-16`
- deep FFT representations: `20, 24`
- deep initial/final cross-FFT native residuals: `1.968e-17` / `3.936e-17`
- full closed-triad tail cross-FFT native residual: `4.995e-16`

## pure upward support, resolution 24
- upward / pure-UV / resolved-contact work: `8.51111181869` / `8.51111181869` / `0`
- pure/contact atoms: `1` / `0`
- max parent/shell ratio: `0.40019526484`
- donor/shell ratio range: `0.40019526484` .. `0.40019526484`

## pure upward support, resolution 28
- upward / pure-UV / resolved-contact work: `8.51111181869` / `8.51111181869` / `0`
- pure/contact atoms: `1` / `0`
- max parent/shell ratio: `0.40019526484`
- donor/shell ratio range: `0.40019526484` .. `0.40019526484`

## deep resolved contact, resolution 20
- initial upward / deep / pure work: `0.11459078671` / `0.0748690219794` / `0.039721764731`
- final upward / deep work: `0.11458218965` / `0.074863405005`
- evolved snapshots with deep upward work: `25`
- max donor excess above M/4: `0.000e+00`
- global NS energy-balance residual: `2.212e-16`

## deep resolved contact, resolution 24
- initial upward / deep / pure work: `0.11459078671` / `0.0748690219794` / `0.039721764731`
- final upward / deep work: `0.11458218965` / `0.074863405005`
- evolved snapshots with deep upward work: `25`
- max donor excess above M/4: `0.000e+00`
- global NS energy-balance residual: `8.809e-19`

## full six-mode Galerkin tail ledger, resolution 20
- initial/final tail energy: `0.666666666667` / `0.666683875767`
- integrated full upward/downward crossing: `8.02119775288e-05` / `2.00708881372e-09`
- normalized tail dissipation: `0.00105001450028`
- inherited/upward common work versus threshold: `0.666666666667` / `8.02119775288e-05` versus `3.15004350083e-05`
- inherited/upward owner flags: `True` / `True`
- tail continuity residual: `9.991e-16`
- worst full cyclic/direct tail-work residual: `5.904e-17`
- worst full Phi_up-Phi_down/signed-tail residual: `3.936e-17`
- maximum internal high-tail circulation: `5.61855735194e-05`
- global Galerkin energy-balance residual: `0.000e+00`

## full six-mode Galerkin tail ledger, resolution 24
- initial/final tail energy: `0.666666666667` / `0.666683875767`
- integrated full upward/downward crossing: `8.02119775288e-05` / `2.00708881372e-09`
- normalized tail dissipation: `0.00105001450028`
- inherited/upward common work versus threshold: `0.666666666667` / `8.02119775288e-05` versus `3.15004350083e-05`
- inherited/upward owner flags: `True` / `True`
- tail continuity residual: `4.995e-16`
- worst full cyclic/direct tail-work residual: `7.872e-17`
- worst full Phi_up-Phi_down/signed-tail residual: `2.952e-17`
- maximum internal high-tail circulation: `5.61855735195e-05`
- global Galerkin energy-balance residual: `0.000e+00`

Deep resolved-scale parent contact is a Fourier-support statement only.  The probe does not call it a smooth-cutoff interface owner.  High-to-high circulation is not reintroduced as tail supply, and no recipient-shell causal reweighting is used.
