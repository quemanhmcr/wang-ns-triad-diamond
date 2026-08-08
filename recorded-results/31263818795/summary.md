# Resolved-role Egorov identity and coherent-deformation redesign

Status: **EXACT_AFFINE_SUBTRACTED_EGOROV_IDENTITY_AND_CENTER_JET_COUNTERMODEL__COHERENT_VARIANCE_REDESIGN**.

For a smooth moving Fourier cell with convolution kernel `K_N`, after transporting the multiplier by the affine jet `A=grad V(X)`, the complete scalar advection commutator is exactly

`(dot Q+[V.grad,Q])f = int K_N(y)[V(x)-V(x-y)-A y].grad f(x-y) dy`.

No higher Taylor expansion is hidden here.  If `||Hess V||<=H`, the exact integral-Taylor bound is

`|V(x)-V(x-y)-A y| <= H(|x-X||y|+|y|^2/2)`.

Including the vector amplitude mismatch `Q[(grad V-A)f]`, fixed dimensionless kernel moments and packet moments give

`||R_Eg f||/||f|| <= (H/N)[k0 Mx+m1 Mxg+(3/2)m2 Mg]`.

However **center Hessian alone is not enough**.  The divergence-free strict-lowpass shear

`V=(0,a[sin(r x1)-sin(2r x1)/2],0)`, `r=N/8`,

has `V(0)=0`, `grad V(0)=0`, `Hess V(0)=0`, but `d1^3 V2(0)=3 a r^3 !=0`.  Thus a bridge claiming that every non-affine moving-role residual is controlled only by the center tensor `B(X)` is false even for a smooth `B_(N/4)` transporter.

The natural replacement is a coherent, affine-invariant deformation observable.  For standard Gaussian coherent coordinate `z`, put

`F(z)=L^-1 grad V(X+Lz)L`,  `Abar=E_gamma F`,
`K_coh^2=E_gamma ||F-Abar||^2`.

Common affine flow makes `K_coh=0`.  Gaussian Poincare, componentwise, gives

`K_coh^2 <= E_gamma ||L^-1 Hess V(X+Lz)[L,L]||^2`.

So the unresolved outer-role defect is not an arbitrary packet error: it is **resolved deformation variance across the physical coherent eddy**.  This catches spatial non-affinity missed by point sampling while preserving affine covariance.

Stress: `50000` exact-commutator/curvature/Hermite/covariance states
- worst exact commutator relative residual: `5.848e-16`
- minimum Taylor-bound margin: `9.715e-02`
- minimum Gaussian-Poincare margin: `0.000e+00`
- worst affine-invariant gradient residual: `8.905e-15`
- countermodel normalized third derivative: `5.859375e-03`

The next theorem must either route this coherent/cellwise curvature variance into the existing H1/H3/source/transfer currencies, or prove that on the certified low-strain/near-extremal branch its contribution is uniformly absorbable by flat erosion.  No such closure is asserted here.
