# Actual Galerkin NS radial spectral crossing audit

Status: **EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_RADIAL_SPECTRAL_CROSSING__ACTUAL_TAIL_STOCK_SIGNED_NONLINEAR_WORK_VISCOSITY__SELECTED_CYCLIC_UPWARD_DOWNWARD_CROSSING__SAME_CUTOFF_CROSS_FFT**.

The probe evolves the actual 2/3-dealiased incompressible Fourier--Galerkin Navier--Stokes system.  It reads full radial-tail stock, signed nonlinear tail work and viscosity at every RK4 output time.  Independently, it reads one actual evolving closed helical triad from the same state and restricts its already-certified cyclic donor flow by the same Fourier radius; the selected triad is not substituted for the full tail law.

- common cutoff: `7`
- radial boundary: `8`
- phase sign: `-1`
- FFT representations: `24, 28`
- tail representation native scale: `2.38714082193`
- selected crossing representation native scale: `0.0111918064899`
- max tail energy/work/viscosity representation residual: `5.581e-16`
- max selected up/down representation residual: `1.550e-16`

## resolution 24
- steps: `64`
- initial/final tail energy: `2.37` / `2.33660511857`
- integrated signed tail nonlinear work: `-0.0171408219284`
- tail viscous dissipation: `0.0162540594788`
- tail interval balance native residual: `1.006e-11`
- integrated selected upward/downward crossing: `0` / `0.00848518440777`
- initial selected upward/downward crossing: `0` / `8.51111181869`
- worst selected radial divergence residual: `0.000e+00`
- worst selected layer-cake residual: `9.994e-17`
- global NS energy-balance residual: `6.692e-12`

## resolution 28
- steps: `64`
- initial/final tail energy: `2.37` / `2.33660511857`
- integrated signed tail nonlinear work: `-0.0171408219284`
- tail viscous dissipation: `0.0162540594788`
- tail interval balance native residual: `1.006e-11`
- integrated selected upward/downward crossing: `0` / `0.00848518440777`
- initial selected upward/downward crossing: `0` / `8.51111181869`
- worst selected radial divergence residual: `0.000e+00`
- worst selected layer-cake residual: `9.996e-17`
- global NS energy-balance residual: `6.692e-12`

The full-tail reading and selected-triad reading are separate physical observables on the same evolved state.  No net-tail Hahn law is minted, no selected triad is claimed to exhaust the full tail, and no radial crossing is promoted to a finite traffic budget or recursive event count.
