# Crossing-to-common-moat extraction

This note removes one previously conditional piece of the PDE-facing packet
bridge.  It does **not** construct a wave-packet decomposition of an arbitrary
Navier--Stokes solution.  It proves that once positive near-extremal transfer
has been resolved into crossing triadic edges at one physical SGS cut, a common
log-frequency moat follows by a four-bin argument; no Gaussian classification
is needed for this scale-localization step.

## 1. Signed-good crossing geometry

Work on the signed-good core used by the physical transfer/Hodge theorem,

\[
\eta_0=10^{-4}.
\]

The certified single-edge theorem gives

\[
u\le 50\eta_0=\frac1{200},\qquad
|g-\gamma_*|\le \sqrt{\eta_0}+25\eta_0=\frac1{80},
\]

where `u` is the parent log-scale imbalance and

\[
g=\ell_c-\ell_{p,top}.
\]

Fix a physical reference cut `tau_0` and retain positive child-transfer edges
which genuinely cross it,

\[
\ell_{p,top}\le \tau_0\le \ell_c.
\]

For the edge midgap

\[
m_e=\frac{\ell_{p,top}+\ell_c}{2}
\]

crossing gives the exact confinement

\[
\boxed{
|m_e-\tau_0|\le \frac{\gamma_*+1/80}{2}.
}
\]

Thus all such midgaps lie in one interval of length `gamma_*+1/80` even if the
absolute Fourier scale of the whole cascade is arbitrarily large.

## 2. Four bins are enough

Split that interval into four equal bins.  For any nonnegative transfer weights
`T_e`, one bin carries at least one quarter of the transfer.  Likewise, for any
nonnegative Hodge residual `D_e`, one bin carries at least one quarter of the
unnormalized Hodge numerator `sum T_e D_e`.

Let `c` be the center of one bin.  Use ideal shell centers

\[
L_p=c-\frac{\gamma_*}{2},\qquad
L_c=c+\frac{\gamma_*}{2}.
\]

The bin halfwidth is `(gamma_*+1/80)/8`; gap error costs another `1/160`, and
the lower parent costs `1/200`.  Therefore both parents and the child lie in
common shells of halfwidth

\[
\boxed{
\sigma_{cross}
=\frac{\gamma_*}{8}+\frac{5}{8}\frac1{80}+\frac1{200}.
}
\]

The 160-bit Arb certificate proves

\[
\boxed{
\sigma_{cross}<\frac2{25}.
}
\]

Hence the previously certified smooth-filter transition `delta=1/20` fits
between the two shells.  In fact

\[
\boxed{
\frac{\gamma_*}{2}-2\sigma_{cross}-\frac1{20}
>\frac1{25}.
}
\]

So a common smooth spectral moat is not an additional Gaussian-profile
hypothesis for crossing edges: it is a deterministic consequence of the
single-edge rigidity plus one physical reference cut.

## 3. Consequence for the physical Hodge branch

The physical change-of-measure theorem gives the good-core coefficient
`25/106`.  If one insists on realizing the Hodge numerator through one common
smooth-cut subcore, the four-bin selection loses at most a factor four.  A
conservative coefficient is therefore

\[
\boxed{
\frac{25}{424}.
}
\]

Thus a crossing block with physical Hodge threshold `h_H` may use the explicit
packet-level cost

\[
\boxed{
c_{0,H}^{cross}
=\min\left\{\frac1{20000},\frac{25}{424}h_H\right\}>0.
}
\]

This does not yet manufacture the triadic edge measure from a general PDE
solution.  It removes the *common log-shell/moat* part of that manufacture once
the physical crossing measure exists.
