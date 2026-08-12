# Actual Galerkin NS radial spectral crossing audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_RADIAL_SPECTRAL_CROSSING__ACTUAL_TAIL_STOCK_SIGNED_NONLINEAR_WORK_VISCOSITY__SELECTED_CYCLIC_UPWARD_DOWNWARD_CROSSING__SAME_CUTOFF_CROSS_FFT**.

The probe evolves the actual 2/3-dealiased incompressible Fourier--Galerkin Navier--Stokes system.  It reads full radial-tail stock, signed nonlinear tail work and viscosity at every RK4 output time.  Independently, it reads one actual evolving closed helical triad from the same state and restricts its already-certified cyclic donor flow by the same Fourier radius; the selected triad is not substituted for the full tail law.

- common cutoff: `7`
- radial boundary: `6.5`
- phase sign: `1`
- FFT representations: `24`
- tail representation native scale: `2.38295949971`
- selected crossing representation native scale: `0.00840205574797`
- max tail energy/work/viscosity representation residual: `0.000e+00`
- max selected up/down representation residual: `0.000e+00`

## resolution 24
- steps: `48`
- initial/final tail energy: `2.37` / `2.37068617015`
- integrated signed tail nonlinear work: `0.0129594995715`
- tail viscous dissipation: `0.0122733295585`
- tail interval balance native residual: `5.739e-11`
- integrated selected upward/downward crossing: `0.00637010588863` / `0`
- initial selected upward/downward crossing: `8.51111181869` / `0`
- worst selected radial divergence residual: `0.000e+00`
- worst selected layer-cake residual: `9.967e-17`
- global NS energy-balance residual: `5.613e-13`

The full-tail reading and selected-triad reading are separate physical observables on the same evolved state.  No net-tail Hahn law is minted, no selected triad is claimed to exhaust the full tail, and no radial crossing is promoted to a finite traffic budget or recursive event count.
