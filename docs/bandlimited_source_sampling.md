# Band-limited source sampling: source replication becomes additive charge

After the H1 source calculus, one remaining loophole was spatial replication: many separated high-frequency grains might sample the **same** filtered pressure, SGS stress, or viscous field and apparently obtain many source exits from one reservoir.  This note removes that loophole at one generation by combining band limitation with Plancherel--Polya sampling.

The theorem uses a standard analytic input: if `f` is Fourier supported in a fixed multiple of `B_N` and points `x_a` are separated by at least `delta/N`, then for `1<p<infinity`

\[
\boxed{
\sum_a |f(x_a)|^p
\le C_{PP}(p,\delta,\Lambda)N^3\|f\|_p^p.
}
\]

This follows from a Schwartz reproducing kernel and the bounded overlap of its translates.  The constants below keep `C_PP` explicit rather than pretending to optimize it numerically.

## 1. Differentiated SGS source samples

For one affine grain, the source normalization obeys

\[
\rho_R
=N^{-4}\big|L^{-1}(\nabla^2\nabla\cdot R)[L,L]\big|
\le \kappa^2s\,N^{-5}|\nabla^3R|,
\qquad s=Nr_g.
\]

On `kappa<=kappa_0`, `s<=s_0`, apply Plancherel--Polya at `p=3/2` to `nabla^3R`.  Since the strict transporter has `supp Rhat subset B_(N/2)`, Bernstein gives

\[
\|\nabla^3R\|_{3/2}
\le C_{D3}N^3\|R\|_{3/2}.
\]

The powers cancel exactly:

\[
-5\cdot{3\over2}+3+3\cdot{3\over2}=0.
\]

Hence

\[
\boxed{
\sum_a\rho_{R,a}^{3/2}
\le
C_{PP}(\kappa_0^2s_0C_{D3})^{3/2}
\|R\|_{3/2}^{3/2}.
}
\]

The Germano increment identity already gives

\[
\|R\|_{3/2}^{3/2}
\le
(1+g_1)^{3/2}g_1^{1/2}Q_{inc}.
\]

Therefore many separated differentiated-SGS source grains force **additive cubic increment charge**.  The Onsager collision then routes it to dominant grain mass, Bellman/cycles, or high-frequency dissipation.

## 2. Viscous-fourth source samples

Similarly

\[
\rho_\nu
\le \kappa_0^2s_0\nu N^{-5}|\nabla^4V|.
\]

Using Plancherel--Polya at `p=2` and the strict support of `V`,

\[
\|\nabla^4V\|_2
\le C_{41}N^3\|\nabla V\|_2.
\]

Since `d_V=N^-1||grad V||_2^2`, again the powers cancel:

\[
-10+3+6+1=0.
\]

Thus

\[
\boxed{
\sum_a\rho_{\nu,a}^2
\le
C_{PP,2}(\kappa_0^2s_0\nu C_{41})^2\,d_V.
}
\]

Spatial replication of a viscous H1 source is therefore paid by additive normalized dissipation.

## 3. Filtered pressure-third near field

For the same strict filter,

\[
P=S_{N/4}p,
\qquad
-\Delta P=\partial_i\partial_j(V_iV_j+R_{ij}).
\]

Hence `P` is itself band limited.  Calderon--Zygmund/Riesz boundedness gives

\[
\|P\|_{3/2}
\le C_R(\|V\|_3^2+\|R\|_{3/2}).
\]

For low-pass `V`, Bernstein gives

\[
\|V\|_3^3
\le C_B^3\mu_V^{3/2},
\qquad
\mu_V=N\|V\|_2^2.
\]

Using `(a+b)^(3/2)<=sqrt(2)(a^(3/2)+b^(3/2))`,

\[
\boxed{
\|P\|_{3/2}^{3/2}
\le
\sqrt2 C_R^{3/2}
\left[C_B^3\mu_V^{3/2}+\|R\|_{3/2}^{3/2}\right].
}
\]

Apply the same `p=3/2` sampling estimate to `nabla^3P`:

\[
\boxed{
\sum_a\rho_{P,a}^{3/2}
\lesssim
\mu_V^{3/2}+Q_{inc}.
}
\]

Thus the **near pressure-third source** is no longer an independent continuum gap for separated source grains: it is paid by resolved low-pass critical mass or by the same SGS increment charge.  The previously proved `6-3=3` multipole estimate continues to handle the far pressure tail.

## 4. Scope

This theorem is a same-generation sampling/packing result.  It does not identify which low-pass packet becomes the persistent ancestor; that is handled by `ancestor_reservoir_sync.md`.  It also assumes a transfer-adapted family of source centers with physical separation comparable to `N^-1`; failure of that separation is a clustering/overlap branch and must be merged before applying the sampling inequality.
