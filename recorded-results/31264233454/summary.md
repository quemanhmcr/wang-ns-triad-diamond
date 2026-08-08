# Coherent affine projection and deformation-dissipation collision

Status: **EXACT_COHERENT_AFFINE_PROJECTION_AND_DEFORMATION_DISSIPATION_COLLISION__AVERAGED_TRANSPORTER_SOURCE_CALCULUS_REMAINS**.

The center-jet countermodel indicates that the affine gauge should be selected by the whole coherent eddy.  In intrinsic Gaussian coordinates `W(z)=L^-1 V(X+Lz)`, define

`vbar=E_gamma W`, `Abar=E_gamma grad W`.

Gaussian integration by parts gives `Abar=E_gamma[W tensor z]`, so

`R=W-vbar-Abar z`

has no Hermite degree 0 or 1.  The Ornstein--Uhlenbeck spectral gap therefore gives exactly

`E|R|^2 <= K_coh^2/2`, where `K_coh^2=E||grad W-Abar||^2`.

Creation/annihilation operators give the weighted bound

`E |z|^2 |R|^2 <= 7 K_coh^2`.

For the coherent Gaussian carrier `psi=g exp(i q.z)`, scalar residual advection plus vector amplitude mismatch obey

`||R_nonaff psi||/||psi|| <= [1+|q|/sqrt(2)+sqrt(7)/2] K_coh`.

Thus the **entire** spatial non-affine Gaussian-core forcing, including all higher Hermite degrees, is controlled by one physical deformation-variance observable.  There is no need to create H4/H5/... master currencies.

Large coherent deformation is already critical dissipation.  The Gaussian density peak, `cond(L)<=567/500`, shell uncertainty `r_g>=2/(3N)`, Cauchy in time and `D_V=N int||grad V||_2^2 dt` give

`I_K^2 <= C_coh c D_V`, `I_K=int K_coh dt`,

with

`C_coh=(567/500)^2 (2 pi)^(-3/2) (3/2)^3 = 0.275568824559`.

Hence `D_V >= I_K^2/(C_coh c)`.  This is a scale-critical dissipation branch, **not** a uniform finite reset count.

On a scale-matched radius branch `s=N r_g<=s0`, the intrinsic carrier is uniformly bounded by

`|q| <= cond(L)^(2/3) s0 (|k|/N)`,

so small `I_K` makes the whole non-affine Gaussian-core low--high residual a controlled perturbation.  Large radius remains the existing sticky affine-radius ancestry branch.

Stress: `50000` spectral/radius/dissipation states
- minimum OU spectral-gap margin: `0.000e+00`
- minimum clean position-weight margin: `-1.776e-15`
- minimum dissipation-collision margin: `9.355e-08`
- minimum radius/aspect axis margin: `6.502e-07`
- clean collision constant: `0.275568824559`

This removes the need for a separate high-Hermite curvature currency at the Gaussian-core forcing level.  The remaining conceptual bridge is now the **source calculus of the coherent averaged affine jet** `Abar(t)`: if it replaces the center jet as the common material transporter, its corotational time variation must be derived from resolved Navier--Stokes and routed to the existing pressure/SGS/viscous/service currencies without double charging.  No global-regularity claim is made.
