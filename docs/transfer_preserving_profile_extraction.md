# One-shot transfer-preserving Gaussian profile extraction

The previous ledger asked for a full iterative profile decomposition of every
near-extremal Navier--Stokes block.  For the no-escape architecture that is
stronger than necessary.  A block which is not near extremal already pays a
positive cost.  On a block which *is* near extremal, the inverse theorem for
sharp Young convolution supplies one dominant Gaussian profile immediately.
This note records the deterministic consequences needed by the physical SGS
ledger.

The analytic external input is Michael Christ's theorem: for
`1<p,q,r<infinity` in the Young range, every sufficiently near extremizer on
`R^d` is close in norm to a Gaussian extremizer (with the usual affine,
translation and modulation symmetries).  No numerical value for Christ's
modulus is invented here.

## 1. No global Fourier `L^{3/2}` assumption is needed

The Navier--Stokes critical exponent appearing in the trilinear convolution is
`p=3/2`.  A frequency-localized finite-energy block automatically belongs to
that space.  If

\[
\operatorname{supp}f\subset\Omega_N,
\qquad |\Omega_N|\le C_\Omega N^3,
\]

then Holder and Plancherel give

\[
\boxed{
\|f\|_{3/2}
\le C_\Omega^{1/6}N^{1/2}\|f\|_2.
}
\]

Thus the inverse Young theorem can be applied block by block even when the full
solution has no known global `||uhat||_{3/2}` bound.

## 2. Weighted near equality reduces to ordinary Young near equality

For

\[
T_m(f,g,h)=\iint m(p,q)f(p)g(q)\overline{h(p+q)}\,dp\,dq,
\qquad |m|\le m_*,
\]

normalize the three `L^{3/2}` norms to one.  If

\[
|T_m(f,g,h)|\ge(1-\delta)m_*A_3,
\]

then simply

\[
|T_m|
\le m_*T_1(|f|,|g|,|h|)
\le m_*A_3
\]

implies

\[
T_1(|f|,|g|,|h|)\ge(1-\delta)A_3.
\]

Hence Christ applies to the magnitudes.  The polarization theorem already in
the repository independently controls the signed phase/orientation defect with
respect to the same transfer-capacity measure.

## 3. One Gaussian profile preserves the transfer

Suppose Christ supplies Gaussian extremizing profiles `F,G,H` with

\[
\|f-F\|_{3/2},\ \|g-G\|_{3/2},\ \|h-H\|_{3/2}\le\varepsilon_G.
\]

The weighted trilinear form has operator norm at most `m_* A_3`.  Telescoping
one input at a time gives the exact normalized remainder bound

\[
\begin{aligned}
\frac{|T_m(f,g,h)-T_m(F,G,H)|}{m_*A_3}
&\le
\varepsilon_G
+(1+\varepsilon_G)\varepsilon_G
+(1+\varepsilon_G)^2\varepsilon_G\\
&=
\boxed{3\varepsilon_G+3\varepsilon_G^2+\varepsilon_G^3}.
\end{aligned}
\]

Therefore

\[
\boxed{
\frac{|T_m(F,G,H)|}{m_*A_3}
\ge
1-\delta-
(3\varepsilon_G+3\varepsilon_G^2+\varepsilon_G^3).
}
\]

At one-percent profile distance the deterministic replacement loss is only

\[
0.030301.
\]

The point is structural, not the decimal: a **single** Gaussian triple captures
all but `o(1)` of the actual weighted transfer as the block approaches equality.
There is no need to extract an arbitrarily long list of profiles inside one
efficient block.

## 4. The Gaussian profile already carries critical `L^2` mass

Because `f` is supported in `Omega_N` and `||F-f||_{3/2}<=epsilon_G`,

\[
\|F1_{\Omega_N}\|_{3/2}\ge1-\varepsilon_G.
\]

Holder on the finite frequency block then gives

\[
\boxed{
N\|F\|_2^2
\ge
C_\Omega^{-1/3}(1-\varepsilon_G)^2.
}
\]

This is exactly the scale-critical local-energy currency of the packet ledger.
After inverse Fourier transform the profile is a spatial Gaussian (possibly an
affine ellipsoid).  A fixed covariance ellipsoid contains a fixed fraction of
its `L^2` mass.  If that ellipsoid is much larger than the natural `N^{-1}`
scale, covering it by natural cells creates the replication/fresh alternative;
if it is scale matched, it is already one critical packet.

For the certified radial log shell used by the crossing theorem,

\[
\left|\log\frac{|\xi|}{N}\right|\le\frac2{25},
\]

the full spherical shell has

\[
|\Omega_N|
\le
\frac{4\pi}{3}
\left(e^{6/25}-e^{-6/25}\right)N^3.
\]

At the concrete profile distance `epsilon_G=1/100`, the 160-bit Arb certificate
therefore gives the clean frozen-time mass bound

\[
\boxed{
N\|G\|_2^2>\frac34.
}
\]

This number is independent of the absolute frequency `N`.  The only
non-explicit input is how small the weighted Young deficit must be for Christ's
theorem to guarantee one-percent Gaussian proximity.

Thus the uncertainty principle gives the desired physical dichotomy:

\[
\boxed{
\text{near-maximal transfer}
\to
\text{one Gaussian critical-mass grain}
\quad\text{or its spatial replication cost}.
}
\]

## 5. Interaction with the new crossing-moat theorem

The Gaussian profile is no longer responsible for producing a common
log-frequency shell.  That task is now cheaper: the crossing-to-common-moat
theorem selects a one-quarter transfer/Hodge subcore purely from the physical
reference cut and the certified single-edge gap rigidity.

The roles are therefore cleanly separated:

1. **Young/Christ:** one coherent phase-space profile preserves almost all
   trilinear transfer;
2. **crossing geometry:** four bins produce a common parent/child scale moat;
3. **single-edge certificate:** the selected transfer measure pays Hodge cost;
4. **Bernstein/Plancherel:** the coherent profile carries critical `L^2` mass;
5. **pressure multipole theorem:** spatial cancellation outside the grain is
   either summable or collides with fresh critical mass.

## 6. Remaining PDE work

What is still not proved is a fully quantitative continuum implementation with
all cutoffs and time dependence synchronized.  In particular one must:

- formulate the three Young roles from the smooth SGS block with a summable
  Coifman--Meyer/wave-packet remainder;
- track the Gaussian profile through a packet lifetime rather than at one frozen
  time;
- make the affine Gaussian grain compatible with the nested ancestry windows;
- insert the resulting critical mass into the fresh/reuse ledger without losing
  the transfer normalization.

But the need for a global Fourier `L^{3/2}` hypothesis and the need for an
infinite within-block profile extraction are both removed.
