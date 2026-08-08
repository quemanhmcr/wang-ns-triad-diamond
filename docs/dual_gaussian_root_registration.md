# Dual-Gaussian root registration: from inverse Young to actual energy quanta

The eventwise Gaussian mark must eventually support the root-energy bound used by
causal binary/Rényi reuse.  A dangerous shortcut would be to say that
`L^(3/2)` proximity to a Gaussian implies `L2` proximity.  It does not.

The correct bridge is duality: a near-Gaussian Young role has a large coefficient
against an explicitly chosen dual Gaussian probe.  This coefficient belongs to
the **actual Navier--Stokes role**, so global energy and Bessel bounds may be used
without introducing Gaussian synthesis coefficients.

## 1. Exact dual Gaussian

Assume first that the magnitude inverse-Young mark and the separate phase/polarization control have been assembled into a **complex phase-aligned Gaussian mark** satisfying

\[
\|f\|_{3/2}=\|G\|_{3/2}=1,
\qquad
\|f-G\|_{3/2}\le\varepsilon.
\]

This complex proximity is an explicit hypothesis of the present theorem.  The repository currently obtains the magnitude Gaussian directly from Christ and controls phase/polarization in separate modules; their quantitative lift to this complex norm statement remains part of the final registration bridge.

For a complex Gaussian profile define

\[
\boxed{h_G=|G|^{-1/2}G.}
\]

Then

\[
\|h_G\|_3=1,
\qquad
\langle G,h_G\rangle=1.
\]

Therefore

\[
\boxed{|\langle f,h_G\rangle|\ge1-\varepsilon.}
\]

This is `L^(3/2)-L^3` duality, not an `L2` approximation theorem.

## 2. Covariance quantization is quantitatively harmless

To use one Bessel theorem for many roots, quantize the Gaussian covariance.
Write the normalized profile as

\[
G_C(\xi)=A_Ce^{-\xi^TC\xi/2}
\]

and the `L3`-normalized dual associated with a representative covariance `D` as

\[
h_D(\xi)=B_De^{-\xi^TD\xi/4}.
\]

If `r_i` are the eigenvalues of `C^{-1/2}DC^{-1/2}`, direct Gaussian integration
gives

\[
\boxed{
\langle G_C,h_D\rangle
=
\prod_{i=1}^3
\left[
\frac{(3/2)r_i^{1/3}}{1+r_i/2}
\right]^{1/2}.
}
\]

Hence if

\[
d_{\log}(C,D)\le\delta,
\]

then every `|log r_i|<=delta`, giving an explicit uniform lower bound.  At

\[
\delta=0.4
\]

the pairing remains close to one.

## 3. L2 normalization produces a critical root quantum

Let `r_{g,rep}` be the geometric physical radius corresponding to the
representative **profile** covariance.  The `L3`-normalized dual satisfies

\[
\boxed{
\|h_D\|_2^2
=
\frac{3\sqrt\pi}{2r_{g,rep}}.
}
\]

Normalize

\[
\phi_D=h_D/\|h_D\|_2.
\]

Then

\[
|\langle f,\phi_D\rangle|
\ge
\frac{P_\delta-\varepsilon}{\|h_D\|_2},
\]

where `P_delta` is the exact Gaussian pairing lower above.

The shell theorem gives

\[
Nr_g>2/3.
\]

If the log-covariance representative lies within `delta`, determinant control
gives

\[
r_{g,rep}
\ge r_g e^{-\sqrt3\delta/6}.
\]

Thus

\[
\boxed{
N|\langle f,\phi_D\rangle|^2
\ge
\frac{2}{3\sqrt\pi}
(P_\delta-\varepsilon)^2
\frac23e^{-\sqrt3\delta/6}.
}
\]

At

\[
\varepsilon=1/100,
\qquad
\delta=0.4,
\]

the right side is **strictly larger than `1/5`**.

This recovers exactly the clean critical root quantum needed by causal reuse, but
as an analysis coefficient of the actual role.

## 4. Divergence-free projection does not alter the coefficient

For a unit polarization vector `e`, use the vector Gaussian probe `phi_D e`.
If `u` is divergence free,

\[
\boxed{
\langle u,\mathbb P(\phi_De)\rangle
=
\langle u,\phi_De\rangle.
}
\]

Thus the root quantum survives the canonical Leray-projected coherent analysis.
The vector family has no worse Gaussian Gram majorant because
`|e_a^*e_b|<=1`.

## 5. Covariance bins are finite on the clean branch

On

\[
2/3<Nr_g\le4,
\qquad
\operatorname{cond}L\le567/500,
\]

the eigenvalues of the normalized physical covariance `N^2 Sigma` lie in the
fixed interval

\[
4/9
<\lambda
\le
(567/500)^{4/3}16.
\]

Thus the normalized SPD matrices lie in a fixed compact subset of the six-dimensional Euclidean space `Sym(3)`.  To control the **affine-invariant** log-SPD metric rather than a commuting log-Euclidean surrogate, let `m=4/9` be the eigenvalue lower bound.  If

\[
\|A-B\|_F\le\epsilon<m,
\]

then for `E=A^{-1/2}(B-A)A^{-1/2}`,

\[
\|E\|_F\le\epsilon/m,
\qquad
\|\log(I+E)\|_F
\le\frac{\epsilon/m}{1-\epsilon/m}.
\]

Therefore choosing

\[
\epsilon=\frac{m\delta}{1+\delta}
\]

guarantees affine-invariant log distance at most `delta`.  A volumetric Frobenius net at this `epsilon` has a finite, scale-independent number of bins.  The numerical cover is intentionally crude; only finiteness independent of `N` and causal depth matters.

Within one covariance bin every root uses the same representative dual covariance.

## 6. Canonical phase-space coloring gives an analysis Bessel budget

For each fixed covariance representative, use the intrinsic phase coordinate of the **dual probe family** and take unit cells in that coordinate.  Color their integer addresses by residues modulo `5` in each coordinate.  There are

\[
5^6
\]

colors, and two distinct cells of one color are at least `4` apart even after
choosing arbitrary points inside the cells.

For each fixed covariance bin and color, the exact affine coherent theorem gives

\[
\boxed{
\sum_a|\langle u,g_a
angle|^2
\le\frac{25}{4}\|u\|_2^2.
}
\]

Therefore all registered root probes together have one finite effective analysis
budget

\[
P_{eff}
\le
(\#\text{covariance bins})\,5^6\,\frac{25}{4}.
\]

The constant is large but scale independent.  In the causal theorem it appears
only inside

\[
\log(P_{eff}E_{global}N_{base}/\eta),
\]

so it changes the finite depth offset but not the positive linear reuse slope.

## 7. Why this is not Gaussian synthesis

No claim is made that

\[
u=\sum_ac_ag_a.
\]

The numbers

\[
\langle u,g_a\rangle
\]

are actual analysis coefficients of `u`, and their Bessel sum is bounded directly
by physical energy.  Christ's Gaussian only determines which dual probe to test.

This distinction is exactly what avoids the old arbitrary redundant-Gaussian
coefficient problem.

## 8. What remains: phase lift and transfer-cell alignment

The theorem gives an actual energy quantum once a **complex phase-aligned** Gaussian event mark is supplied.  The current inverse-Young theorem directly marks the magnitudes, so the first remaining step is to combine the existing phase/polarization rigidity with that magnitude mark quantitatively.  After that, a causal parent slot still comes from an **actual positive transfer-selected coherent work atom**.

The final registration theorem must therefore show one of:

1. a fixed transfer-weighted fraction of the parent slots lie in the coherent
   root cells carrying the dual-Gaussian quantum above; or
2. the misaligned transfer is already one of the named physical currencies:
   selected cross/relink work, backscatter/cancellation, source/service, or a
   multiplicative profile/transfer loss.

A measurable Christ selector by itself would not prove this.  The remaining issue
is physical **transfer-cell alignment**, not selection theory.

No global regularity claim is made.
