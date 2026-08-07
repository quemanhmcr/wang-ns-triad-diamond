# Affine curvature connection and sideband source identity

For any smooth transporter `V`, along `dot X=V(X)`, `dot L=A L`, let `A=grad V`, `H=Hess V`, `F=D_t^V V`, and `B=L^-1 H[L,L]`.  Exact differentiation gives

`dot B+2 A_aff B = L^-1 Hess(F)[L,L]`,  `A_aff=L^-1 A L`.

Thus common affine deformation is a connection on curvature rather than a source.  For the resolved Navier--Stokes transporter
`F=-grad P-div R+nu Delta V`,

`Hess(F)=-nabla^3 P-nabla^2 div R+nu nabla^2 Delta V`.

Hence curvature-sideband dephasing is sourced by pressure third derivatives, differentiated SGS stress, or viscous fourth velocity derivatives after the affine connection is removed.  The pressure kernel is homogeneous of degree `-3`; three derivatives have degree `-6`, and 3D packet packing leaves the summable far exponent `6-3=3`.

This is a source/locality theorem, not yet a daughter-grain/coherence cost.

Stress checks: `50000`
- worst normalized connection residual: `1.361e-13`
- worst resolved-source split residual: `0.000e+00`
- worst `6-3=3` homogeneity residual: `0.000e+00`
