# Multi-child continuum edge measure on one actual Galerkin Navier--Stokes orbit

Resolution `24`, common cutoff `5`, child modes `((5, 1, 0), (3, 4, 1), (1, 5, 2))`, `32` RK4 steps, `3` sampled PDE times.

Every child is reconstructed from all retained ordered parents, then unordered parent orbits and all eight helicity sectors. The outer child family is aggregated **before** the Hahn identity is checked on the same evolved PDE state.

- child registrations: `9`;
- minimum unordered pairs/helical edges per child: `269` / `2152`;
- worst child source/work/progress/Hahn residuals: `2.838e-16`, `3.953e-16`, `6.956e-16`, `3.781e-16`;
- worst outer-child work/modal/Hahn/progress residuals: `2.115e-16`, `1.057e-16`, `1.057e-16`, `4.104e-16`;
- minimum outer-child positive-Hahn dominance margin: `7.890e-03`;
- child snapshots with positive nonforward physical work: `9`;
- NS energy-balance residual: `6.822e-10`;
- maximum normalized divergence: `4.817e-17`.

This is a finite Galerkin outer-child aggregation falsifier of the companion analytic joint Radon theorem `dLambda_unord=(1/16) dz d(q_#dr)`.  The PDE probe does not replace that proof; it attacks its normalization, parent quotient, helicity reconstruction, and Hahn ordering on one evolved Navier--Stokes orbit.
