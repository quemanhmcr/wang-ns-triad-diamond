# Sharp-cutoff log-scale energy-flux bridge

This note gives an exact finite-triad bridge from the progress factor in
`mathcal J` to a genuine Navier--Stokes energy-flux quantity.  It does **not**
yet solve the continuum localization/atomic-extraction problem.

## 1. Sharp spectral flux of one conservative triad

Let an ordered Fourier triad have magnitudes

\[
0<k\le p\le q
\]

and nonlinear modal energy rates

\[
\dot E_k+\dot E_p+\dot E_q=0.
\]

For the sharp spectral low-pass energy, define outward nonlinear flux by

\[
\Pi_K^{(e)}
=-\frac{d}{dt}_{\!NL}
\sum_{|\xi|\le K}E_\xi
\]

restricted to this triad.  Purely by which modes lie below the cutoff,

\[
\Pi_K^{(e)}=
\begin{cases}
0,&K<k,\\
-\dot E_k,&k\le K<p,\\
\dot E_q,&p\le K<q,\\
0,&K\ge q.
\end{cases}
\]

The third line uses energy conservation:
`-(dot E_k+dot E_p)=dot E_q`.  Therefore the Mellin/log-scale moment is

\[
\boxed{
\int_0^\infty \Pi_K^{(e)}\,\frac{dK}{K}
=-\dot E_k\log\frac pk
+\dot E_q\log\frac qp.
}
\]

This is exact for a finite Fourier/Galerkin triad.  No turbulence closure or
statistical assumption is used.

## 2. Why the progress factor is natural, and the correction we must keep

A tempting but false identity would keep only the second term.  The interval
`k<=K<p` contributes the first term whenever the parent magnitudes are unequal.
Thus the correct physics ledger contains two scale segments.

At the single-edge extremizer `k=p`, however,

\[
\log(p/k)=0,
\]

so the lower segment disappears exactly and

\[
\boxed{
\int_0^\infty \Pi_K^{(e)}\frac{dK}{K}
=\dot E_q\log\frac qp.
}
\]

Hence the logarithmic progress in `mathcal J` is not an ad hoc reward for
jumping scales: it is the exact logarithmic measure of the cutoff interval
across which the child receives energy.  The cusp direction measuring parent
imbalance is simultaneously the direction that turns on the omitted lower-cut
segment.

## 3. Helical maximizing orbit and flux retention

Normalize `q=1`, order `0<x<=y<1`, and take opposite-helicity parents
`(s_x,s_y)=(+,-)`.  For the child sign `s_q=-`, which realizes the sign envelope
when `x<y`, the common triad phase/amplitude factor `R` gives

\[
\dot E_x=(1-y)R,
\qquad
\dot E_y=-(1+x)R,
\qquad
\dot E_q=(x+y)R.
\]

Choosing the phase so the child gains energy means `R>0`.  The lower cutoff
segment is then a small backscatter term, while the upper segment is forward:

\[
\frac1R\int_0^\infty\Pi_K^{(e)}\frac{dK}{K}
=(x+y)\log(1/y)
-(1-y)\log(y/x).
\]

The first term is precisely the progress times the child-energy coefficient
that enters the sign-maximized single-edge functional; the geometric helical
factor `|g|` and the packet amplitude/phase factor are what become the physical
transfer weight in a block ledger.

In the certified Hodge coordinates

\[
u=\log(y/x),
\qquad
v=-\frac12(\log x+\log y)-\gamma_*,
\]

with `0<=u<=2/25`, `|v|<=2/25`, the adverse lower/upper ratio obeys

\[
\eta_{low}
=\frac{(1-y)u}{(x+y)\log(1/y)}<\frac1{10}.
\]

The Arb certificate checks this using only the certified bracket for `r_*` and
monotone endpoint bounds.  Consequently

\[
\boxed{
\int_0^\infty\Pi_K^{(e)}\frac{dK}{K}
\ge\frac9{10}\,
\dot E_q\log\frac1y
}
\]

for the adverse maximizing child sign.  For the opposite child sign the lower
segment is forward and adds to the upper contribution, so this is the worst
case among the two child helicities near the equal-parent extremizer.

This is the first direct PDE-facing reason that the same anisotropic stability
certificate is useful twice: `u` both costs multiplier efficiency and controls
the difference between the upper-edge progress ledger and the full signed
sharp-cutoff Mellin flux.

## 4. What this does and does not prove

The identity upgrades the status of the logarithmic factor: it is forced by
sharp spectral energy balance.  It also gives an exact finite-triad way to
factor a near-extremal block into

\[
\text{physical transfer weight}
\times
\text{dimensionless helical progress coefficient}.
\]

A later calculation showed that the **full all-scale Mellin moment is not the
right replacement for `mathcal J`**: the rational triad
`x=13/40`, `y=17/20`, helicities `(-,+,-)` has a full Mellin coefficient
Arb-certified above `3J*/2`, because lower-scale redistribution dominates.
Thus the all-scale moment is retained only as a conservation diagnostic.

The correct PDE-facing observable is developed in
`docs/smooth_log_flux_cocycle.md`. An even graded log filter and a
transfer-weighted common midgap give the exact identity

\[
2\int_\tau^\infty\Pi^\delta_{\rm core}(t)dt
=\sum_eT_e\log(q_e/p_e),
\]

with **no sharp-to-smooth comparison loss** on a common spectral moat. The
same module identifies the graded spectral flux with the space-average physical
SGS transfer and converts near-saturation into positive physical transfer
weights comparable to the Bellman/Hodge capacity measure.

What remains is therefore genuinely packet/PDE-level: construct the transfer
core and its common shell/moat from an arbitrary near-extremal block, then
control spatial-window commutators, pressure boundary work, backscatter and time
synchronization with summable errors.
