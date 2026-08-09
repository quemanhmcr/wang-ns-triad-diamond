# Exact outer moving-role extraction

Status: **EXACT_OUTER_MOVING_ROLE_IDENTITY_AND_PERSISTENT_LOWLOW_EXCLUSION__INTERFACE_PROVENANCE_SUPPLIED**.

Let `V=S_(N/4)u`, `h=u-V`, and `L_V f=P div(V tensor f+f tensor V)`.  For a time-dependent scalar divergence-free Fourier role `w=Q(t,D)u`, direct differentiation of Leray Navier--Stokes gives exactly

`(dt+L_V-nu Delta)w = Q B(V,V) - Q B(h,h) + (dtQ+[L_V,Q])u`.

There is no pressure source, and scalar Fourier `Q` commutes with `Delta`.  Anchor `Q` at the already selected frozen physical transfer cell and transport its scalar symbol by the coherent affine dual flow.  On the low-strain branch `K<=1/30`, an anchored lower edge `3N/5` stays above `(3/5)e^(-1/30)N>N/2`, whereas `V tensor V` is supported in `B_(N/2)`.  Hence `Q B(V,V)=0` **throughout the whole role interval**, not merely at the anchor slice.

The selected outer-role equation is therefore

`(dt+L_V-nu Delta)w = -Q P div(h tensor h) + R_Q`,
`R_Q=(dtQ+[L_V,Q])u`.

The first term is the unique genuine quadratic high--high source; its positive coherent Hahn disintegration is already the physical binary causal work law.  For the second term, affine dual transport cancels the common affine advection Heisenberg term exactly, and constant affine stretching commutes with scalar `Q`.  The exact Egorov identity shows that `R_Q` is purely **nonaffine resolved low--high role-interface work** and vanishes for a genuinely affine transporter.

Stress: `50000` exact algebra/support/Egorov states
- worst outer-role identity residual: `5.245e-15`
- worst scalar-Q/viscosity commutator: `0.000e+00`
- worst affine Heisenberg residual: `0.000e+00`
- worst affine-subtracted Egorov residual: `5.995e-16`
- minimum sampled persistent low-low gap: `8.608317e-02`

This changes the continuum frontier: constructing the moving outer role itself is no longer the missing PDE step.  The only unclosed outer-role term is the **work-level routing of the nonaffine Heisenberg interface**.  It must be shown to enter exactly once as either coherent deformation/critical `D_V` or physical role-relink/transfer loss; it must not be promoted to a new currency or silently absorbed as representation `Xi`.  No Navier--Stokes global-regularity conclusion is asserted.
