# Actual Galerkin NS hard-tail true upward supply audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HARD_TAIL_TRUE_UPWARD_SUPPLY__FULL_TAIL_STOCK_WORK_VISCOSITY_SEPARATE_FROM_SELECTED_CYCLIC_SUPPORT__PURE_UV_FIRST_SHELL_AND_DEEP_RESOLVED_CONTACT__NO_INTERFACE_OVERCLAIM**.

The full-tail stock/work/viscosity trajectory remains the certified radial PDE reading.  Selected cyclic triads are separate sub-observables used only to test parent-support geometry; they are never substituted for the full tail law.

- main FFT representations: `24, 28`
- pure-support cross-FFT native residual: `1.900e-16`
- deep FFT representations: `20, 24`
- deep initial/final cross-FFT native residuals: `2.688e-17` / `5.375e-17`
- full closed-triad tail cross-FFT native residual: `1.970e-16`

## pure upward support, resolution 24
- upward / pure-UV / resolved-contact work: `18.6989126657` / `18.6989126657` / `0`
- pure/contact atoms: `1` / `0`
- max parent/shell ratio: `0.40019526484`
- donor/shell ratio range: `0.40019526484` .. `0.40019526484`

## pure upward support, resolution 28
- upward / pure-UV / resolved-contact work: `18.6989126657` / `18.6989126657` / `0`
- pure/contact atoms: `1` / `0`
- max parent/shell ratio: `0.40019526484`
- donor/shell ratio range: `0.40019526484` .. `0.40019526484`

## deep resolved contact, resolution 20
- initial upward / deep / pure work: `0.251755958403` / `0.164487241289` / `0.0872687171139`
- final upward / deep work: `0.251731771878` / `0.164471438784`
- evolved snapshots with deep upward work: `33`
- max donor excess above M/4: `0.000e+00`
- global NS energy-balance residual: `2.967e-17`

## deep resolved contact, resolution 24
- initial upward / deep / pure work: `0.251755958403` / `0.164487241289` / `0.0872687171139`
- final upward / deep work: `0.251731771878` / `0.164471438784`
- evolved snapshots with deep upward work: `33`
- max donor excess above M/4: `0.000e+00`
- global NS energy-balance residual: `3.645e-16`

## full six-mode Galerkin tail ledger, resolution 20
- initial/final tail energy: `1.12666666667` / `1.12675633878`
- integrated full upward/downward crossing: `0.000226577519868` / `9.47567173902e-09`
- normalized tail dissipation: `0.00228159886111`
- inherited/upward common work versus threshold: `1.12666666667` / `0.000226577519868` versus `6.84479658334e-05`
- inherited/upward owner flags: `True` / `True`
- tail continuity residual: `1.576e-15`
- worst full cyclic/direct tail-work residual: `8.956e-17`
- worst full Phi_up-Phi_down/signed-tail residual: `3.583e-17`
- maximum internal high-tail circulation: `0.000206304406891`
- global Galerkin energy-balance residual: `0.000e+00`

## full six-mode Galerkin tail ledger, resolution 24
- initial/final tail energy: `1.12666666667` / `1.12675633878`
- integrated full upward/downward crossing: `0.000226577519868` / `9.47567173902e-09`
- normalized tail dissipation: `0.00228159886111`
- inherited/upward common work versus threshold: `1.12666666667` / `0.000226577519868` versus `6.84479658334e-05`
- inherited/upward owner flags: `True` / `True`
- tail continuity residual: `1.773e-15`
- worst full cyclic/direct tail-work residual: `8.958e-17`
- worst full Phi_up-Phi_down/signed-tail residual: `5.375e-17`
- maximum internal high-tail circulation: `0.000206304406891`
- global Galerkin energy-balance residual: `0.000e+00`

Deep resolved-scale parent contact is a Fourier-support statement only.  The probe does not call it a smooth-cutoff interface owner.  High-to-high circulation is not reintroduced as tail supply, and no recipient-shell causal reweighting is used.
