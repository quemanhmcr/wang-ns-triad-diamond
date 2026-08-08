# Dual-Gaussian root registration: Christ proximity gives actual analysis energy

Status: **EXACT_DUAL_GAUSSIAN_ANALYSIS_ROOT_QUANTUM__TRANSFER_CELL_ALIGNMENT_REMAINS**.

A phase-aligned complex Gaussian near-profile should not be inserted into the velocity as a fictitious synthesis component.  Use it only to design a dual Gaussian **analysis probe**.  Christ supplies the magnitude Gaussian; the present theorem is conditional on the separate phase/polarization control being strong enough to lift that mark to complex `L^(3/2)` proximity.

Let `||f||_(3/2)=||G||_(3/2)=1`, `||f-G||_(3/2)<=eps`.  The exact L3-dual Gaussian is `h_G=|G|^(-1/2)G`, with `||h_G||_3=1` and `<G,h_G>=1`.  Quantize the profile covariance to a representative within log-SPD radius `delta`, and use the corresponding L3-normalized dual `h_rep`.  If `delta<=0.4`, exact Gaussian integration gives

`<G,h_rep> >= 0.97266240464`.

Therefore `|<f,h_rep>| >= <G,h_rep>-eps`.  After L2-normalizing `phi=h_rep/||h_rep||_2`,

`||h_rep||_2^2 = 3 sqrt(pi)/(2 r_g,rep)`.

The shell lower axis and covariance quantization give the scale-independent actual coefficient bound

`N |<f,phi>|^2 >= 0.207034423405`

at `eps=1/100`, `delta=0.4`.  In particular it is **strictly larger than 1/5 for the normalized role**.  For the actual physical role `f=a f_hat`, the bound is multiplied by `a^2=||f||_(3/2)^2`.  Therefore this theorem alone does **not** provide the absolute root quantum required by causal/Renyi when parent amplitudes are allowed to shrink.

After restoring the physical amplitude this is a coefficient of the actual selected role.  If that role is obtained by an exact self-adjoint outer Fourier/helical projector `Q`, use the probe `Q phi`; then `<u,Q phi>=<Qu,phi>` exactly and `||Q phi||_2<=1`, so normalization cannot reduce the coefficient.  For divergence-free roles, Leray projection is likewise coefficient-preserving.  No L2 closeness of the Christ remainder is asserted or needed.

Variable root scale also does not destroy the energy count.  Put `log N` in bins of width `2/25` and color the bins modulo `4`.  Because the outer role shell halfwidth is `2/25`, two distinct bins of one color have disjoint physical Fourier support, so their exact outer role projectors are orthogonal.  Within one scale bin use its reference `N_b`; the rescaled covariances `N_b^2 Sigma` lie in a fixed compact subset of six-dimensional `Sym(3)`.  A crude Frobenius net, chosen fine enough to guarantee affine-invariant log-SPD radius `delta=0.4`, uses at most `49876112440099344` bins.  For each fixed scale/covariance representative, unit cells in that representative dual-probe phase coordinate may be colored with `5^6=15625` colors so cells of one color are 4-separated.  This auxiliary coloring is only an analysis-budget device; it is not identified with the canonical material label.  Inside one scale bin/covariance/color the exact affine coherent Bessel theorem gives analysis budget `25/4`; orthogonality across same-color scale bins prevents any factor growing with causal depth.  Thus all registered roots have one finite effective budget

`P_eff <= 1.94828564219e+22`.

This constant is intentionally huge but **scale independent**.  In the causal root estimate it enters only through `log P_eff`, so it changes the finite depth offset and not the linear reuse slope `log(48/25)`.

The theorem therefore closes a subtle **shape-to-analysis** gap once complex Gaussian proximity is available: it produces a scale-critical coefficient **relative to the parent `L^(3/2)` amplitude** in the actual `L2` energy analysis, via duality, without pretending `L^(3/2)` closeness implies `L2` closeness.

Stress: `50000` covariance/pairing/Leray states
- minimum Gaussian pairing margin: `7.007e-04`
- minimum root-mass margin above the default endpoint: `9.590e-05`
- minimum covariance-cover radius margin: `-1.110e-16`
- worst Leray coefficient residual: `1.387e-14`

Two registration issues remain.  First, the current inverse-Young ledger directly supplies Gaussian proximity for magnitudes; the separate phase/polarization theorems must be assembled into the complex phase-aligned proximity assumed here.  Second,  A causal parent slot already has a transfer-selected material coherent label.  The dual probe above is centered at the Christ Gaussian mark.  The final assembly must prove that a fixed transfer-weighted fraction of roots are aligned with that mark (or else the misaligned physical work is already cross/relink/backscatter/service currency).  Measurable selection of a Christ mark is not the hard part; this **transfer-cell alignment** is.  No global-regularity claim is made.
