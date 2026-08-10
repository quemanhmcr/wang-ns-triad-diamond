# Hermite-helicity forcing ledger: field residual is not the same as spinor forcing

The forced relative-polarization identity was written for two-component helical
spinors.  Once an actual affine Gaussian packet is introduced, a PDE residual is
a **field** in `L^2`, not automatically a two-component spinor.  It must first be
projected onto the Hermite/helicity packet basis.

Let `g(z)` be the normalized affine Gaussian envelope, so `|g|^2 dz` is standard
Gaussian measure, and let `h_+,h_-` be the two constant triad-normal helical
polarizations at the packet carrier.  The degree-zero packet space is

\[
\mathcal H_0=\operatorname{span}\{g h_+,g h_-\}.
\]

Only

\[
\boxed{F_i^{(0)}=P_{\mathcal H_0}f_i}
\]

is the `F_i` appearing in the exact forced symplectic spinor identity.

## 1. Third-Hermite curvature is an orthogonal daughter mode

For a fully symmetric three-tensor `T`, define the probabilists' third Hermite
tensor

\[
H_{abc}^{(3)}(z)
=z_az_bz_c-(\delta_{ab}z_c+\delta_{ac}z_b+\delta_{bc}z_a).
\]

Gaussian chaos orthogonality gives

\[
\boxed{
\mathbb E[(T:H_3)P_2(Z)]=0
}
\]

for every scalar polynomial `P_2` of degree at most two, and the exact isometry

\[
\boxed{
\mathbb E[(T:H_3)^2]=3!\,\|T\|_F^2.
}
\]

The affine nonconformal envelope residual derived previously is

\[
F_{env}=-\frac14(\operatorname{Sym}B):H_3\,g,
\]

hence

\[
\boxed{
\frac{\|F_{env}\|_2^2}{\|g\|_2^2}
=\frac38\|\operatorname{Sym}B\|_F^2.
}
\]

It has **zero projection** onto the base spinor and onto all scalar Gaussian
parameter tangent modes.  Therefore it is not a direct phase/polarization force.
It is a coherent-leakage/daughter-packet channel.

## 2. Spatial polarization curvature is a first-Hermite sideband

If the local transfer-distinguishable polarization generator varies as

\[
D(z)-D(0)=\sum_c C_cz_c,
\]

then

\[
\boxed{
\mathbb E\left|\sum_cC_cZ_c\right|^2=\sum_c\|C_c\|_F^2,
\qquad
\mathbb E\left[\sum_cC_cZ_c\right]=0.
}
\]

Thus this curvature creates an `H_1` **polarization sideband**.  Its mean
projection into the base Gaussian spinor is zero, while the RMS theorem from
`affine_polarization_curvature` proves that the transfer-relevant sideband sees
at least one half of the physical strain-variation energy.

## 3. Correct forcing hierarchy

The packet PDE residual should therefore be split as

\[
\boxed{
f_i
=f_i^{(0)}
+f_i^{H_1,pol}
+f_i^{tan,\le2}
+f_i^{H_3,env}
+f_i^{higher}.
}
\]

- `f_i^(0)`: true base-spinor forcing; this alone enters `R_F` in the forced
  symplectic identity;
- `H_1,pol`: spatial polarization sideband; route to coherence/branching;
- scalar tangent `<=2`: center/carrier/covariance/chirp and bulk viscosity;
  quotient it;
- `H_3,env`: non-affine envelope daughter mode;
- higher modes: unresolved packet forcing/cross interactions.

This corrects the tempting but nonphysical shortcut of bounding a whole field
residual by `||f_i||_2` and inserting it directly as a spinor force.  Such a
bound is safe as an upper estimate but destroys the orthogonal structure and can
double-count curvature as both polarization forcing and packet leakage.

The remaining PDE task is now sharper: estimate the **degree-zero projection** of
role-dependent high-high/cross-cell forcing, while treating the orthogonal
Hermite sidebands through the fresh/reuse/coherence ledger.

## Current downstream status

The degree-zero/source-routing questions raised in this precursor were absorbed into the later covariant source, coherent averaged-strain, source compiler, and hard-event role architecture.  Higher Hermite towers are not a canonical closure variable in the present continuum spine.  This note is retained for the forcing/helicity provenance that motivated those later reductions.
