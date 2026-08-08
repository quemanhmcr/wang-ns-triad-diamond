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

This complex proximity is the input to the duality calculation.  The companion `complex_young_parent_marking.md` shows how the existing physical symbol-freezing `Xi` reduces one efficient frozen scalar/helical parent cell to ordinary **complex** Young near-extremality; Christ's complex-valued stability theorem then supplies exactly this input, with its external modulus left symbolic.

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

If the selected parent is the range of an exact self-adjoint outer Fourier/helical projector `Q`, use `Q phi_D` as the probe.  Then `<u,Q phi_D>=<Qu,phi_D>` and `||Q phi_D||_2<=1`, so after normalization the lower coefficient cannot decrease.  Thus neither unrelated Fourier components nor pressure can cancel the registered root quantum.
The vector family has no worse Gaussian Gram majorant because
`|e_a^*e_b|<=1`.

## 5. Scale colors and covariance bins give a uniform global budget

The root scales in one depth layer need not be identical, so it would be wrong
to apply an equal-physical-covariance Bessel theorem to all roots at once.
Instead use the **outer role projector** which already carries the physical
frequency shell.

Put `log N` into bins of width

\[
h=2/25.
\]

A role at scale `N` is supported in the shell

\[
[e^{-2/25}N,e^{2/25}N].
\]

The union over one scale bin `b` is therefore

\[
[e^{bh-2/25},e^{(b+1)h+2/25}].
\]

Color scale-bin indices modulo `4`.  Since

\[
(4-1)h>2(2/25),
\]

distinct bins of one color have **disjoint physical frequency supports**.  Their
exact outer Fourier/helical projectors are orthogonal.  Consequently a Bessel
sum over all bins of one color does not acquire a factor equal to the number of
causal scales.

Inside one scale bin take the reference

\[
N_b=e^{(b+1/2)h}.
\]

Then `N/N_b in [e^{-h/2},e^{h/2}]`.  On

\[
2/3<Nr_g\le4,
\qquad
\operatorname{cond}L\le567/500,
\]

the eigenvalues of the **bin-rescaled** physical covariance `N_b^2 Sigma` lie in
one fixed interval

\[
e^{-h}\frac49
<\lambda
\le
e^h(567/500)^{4/3}16.
\]

Thus these matrices lie in a fixed compact subset of six-dimensional `Sym(3)`.
To control the affine-invariant log-SPD metric, let `m` be the lower eigenvalue
bound.  If

\[
\|A-B\|_F\le\epsilon<m,
\]

then

\[
\|\log(A^{-1/2}BA^{-1/2})\|_F
\le
\frac{\epsilon/m}{1-\epsilon/m}.
\]

Choosing

\[
\epsilon=\frac{m\delta}{1+\delta}
\]

therefore gives affine-log radius at most `delta`.  A volumetric Frobenius net at
`delta=0.4` has a finite, scale-independent number of covariance representatives.

## 6. Phase-space coloring and Bessel

For each fixed scale bin and covariance representative, use the intrinsic phase
coordinate of the **dual probe family** and take unit cells.  Color their integer
addresses modulo `5` in each coordinate.  There are

\[
5^6
\]

colors, and one color is `4`-separated.

Within one scale bin/covariance/color the exact affine coherent theorem gives

\[
\sum_a|\langle v,g_a\rangle|^2
\le\frac{25}{4}\|v\|_2^2.
\]

Apply this with `v=Q_bu`, where `Q_b` is the exact outer role projector for the
scale bin (and helicity sector, if one is used).  Across scale bins of one color,
the `Q_b` ranges are orthogonal.  Summing first over those bins and then over the
four scale colors gives one **depth-independent** effective analysis budget

\[
P_{eff}
\le
4\,(\#\text{covariance bins})\,5^6\,\frac{25}{4}.
\]

The probe may carry a fixed unit polarization depending on the selected frozen
scalar/helical role; this does not worsen the scalar Gaussian Gram majorant.

The constant is huge but scale independent.  In the causal root estimate it
appears only inside

\[
\log(P_{eff}E_{global}N_{base}/\eta),
\]

so it changes the finite depth offset, not the positive linear reuse slope.

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

## 8. What remains: transfer-cell alignment

The theorem gives an actual energy quantum once a complex Gaussian parent mark is supplied.  The companion complex-Young reduction obtains that mark from a frozen efficient physical cell using the same symbol `Xi` and Christ's complex-valued theorem.  A causal parent slot, however, still comes from an **actual positive transfer-selected coherent work atom**.

The final registration theorem must therefore show one of:

1. a fixed transfer-weighted fraction of the parent slots lie in the coherent
   root cells carrying the dual-Gaussian quantum above; or
2. the misaligned transfer is already one of the named physical currencies:
   selected cross/relink work, backscatter/cancellation, source/service, or a
   multiplicative profile/transfer loss.

A measurable Christ selector by itself would not prove this.  The remaining issue
is physical **transfer-cell alignment**, not selection theory.

No global regularity claim is made.

## Superseding amplitude-entropy use

The normalized dual coefficient is **not** promoted to an absolute root mass.
The physical role carries the factor `||f||_(3/2)^2`, exactly as required by
homogeneity.  The preferred downstream theorem keeps that amplitude explicitly:
`alpha=sqrt(N)|<u,phi>|` enters the multiplicative parent-product law and then the
root entropy estimate.  Thus no separate parent-amplitude floor is required.
The remaining continuum issue is the outer selected-role construction, not a
new Gaussian-registration mass hypothesis.
