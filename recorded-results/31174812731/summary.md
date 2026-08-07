# Exact relative-polarization transport

The identities are exact; random checks validate the implementation.

- exact parent wedge law: `Wdot = U^T J (D1-D2) V`
- exact polarization numerator law includes only relative parent generator and child generator
- pointwise capacity bound: `|Pdot| <= 2 sqrt(||D1-D2||_F^2+||D3||_F^2)||U||||V||||Z||`
- random differential/time-ordered checks: `50000`
- worst wedge RHS residual: `9.653e-09`
- worst polarization RHS residual: `1.556e-08`
- minimum pointwise bound margin: `1.760e-01`
- worst arbitrary common time-ordered wedge residual: `5.661e-15`
- minimum nonlinear-forcing bound margin: `1.038e+00`
- hyperbolic common-gauge countermodel: `||M-I||=2.980e+03`, `cond(M)=8.886e+06`, wedge residual `1.462e-10`

Thus the full time-ordered parent observable does not require a Magnus expansion:
common `SL(2)` motion cancels pointwise in the symplectic wedge.  Euclidean
propagator distance from the identity is not a physical polarization defect.
