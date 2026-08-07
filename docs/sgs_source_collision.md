# Filtered differentiated-SGS and viscous source collision

The H1 covariant source theorem routes dephasing into pressure third derivatives, differentiated SGS stress, or viscous fourth derivatives. This note closes the first collision step for the latter two channels without differentiating the raw turbulent velocity.

The physical principle is that the curvature transporter is already coarse-grained:

\[
V=S_{N/4}u,
\qquad \operatorname{supp}\widehat V\subset B_{N/4}.
\]

Let

\[
R=S_{N/4}(u\otimes u)-V\otimes V.
\]

Then

\[
\boxed{\operatorname{supp}\widehat R\subset B_{N/2}.}
\]

Thus derivatives of the SGS source are derivatives of a **filtered stress**, not derivatives of raw turbulence.

## 1. Differentiated SGS source is controlled by the critical stress norm

Under the unitary Fourier convention, the Hilbert-valued symbol

\[
R\mapsto \nabla^2\nabla\cdot R
\]

has operator norm at most `|xi|^3`. Hausdorff--Young from `L^(3/2)` to `L^3`, followed by Holder in the ball `|xi|<=N/2`, gives

\[
\|\nabla^2\nabla\cdot R\|_\infty
\le C_RN^5\|R\|_{3/2},
\]

where

\[
C_R=(2\pi)^{-2}
\left(\frac{8\pi}{15}\right)^{2/3}
\left(\frac12\right)^5
<\frac1{800}.
\]

For the affine normalized source

\[
S_R=L^{-1}(-\nabla^2\nabla\cdot R)[L,L],
\]

put

\[
\Lambda_L=N\|L^{-1}\|\|L\|^2,
\qquad s=Nr_g,
\qquad r_g=|\det L|^{1/3}.
\]

If `kappa=cond(L)`, singular-value algebra gives

\[
\Lambda_L\le \kappa^2s.
\]

On the physical H1 mild-aspect branch `kappa<=21/20`, Arb certifies `kappa^2<6/5`, hence

\[
\boxed{
N^{-4}\|S_R\|
\le \frac3{2000}s\|R\|_{3/2}.
}
\]

Therefore if `s<=s0` and `rho_R=N^-4||S_R||`,

\[
\boxed{
\|R\|_{3/2}
\ge \frac{2000}{3s_0}\rho_R.
}
\]

Using the exact Germano increment inequality from `affine_sgs_boundary_ledger`,

\[
|R|^{3/2}
\le C_{inc}(G)
\int |G_\ell(r)|\,|\delta_ru|^3dr,
\]

this forces an Onsager cubic increment charge at the actual transporter filter scale `ell~4/N`:

\[
\boxed{
\iint |G_\ell(r)|\,|\delta_ru|^3\,drdx
\ge
\frac1{C_{inc}(G)}
\left(\frac{2000\rho_R}{3s_0}\right)^{3/2}.
}
\]

If instead `s>s0`, this is not called an aspect defect. The existing affine critical-grain theorem gives

\[
\boxed{
N\int_{E_2}|u|^2\ge\frac3{10}s_0,
}
\]

so the block is handed to fresh/reuse radius-energy ancestry.

## 2. Viscous fourth derivatives are a dissipation-density event

Since `supp Vhat` lies in `B_(N/4)`, Cauchy--Schwarz in Fourier space yields

\[
\|\nabla^2\Delta V\|_\infty
\le C_\nu N^{9/2}\|\nabla V\|_2,
\]

with

\[
C_\nu=(2\pi)^{-3/2}
\sqrt{\frac{4\pi}{9}}
\left(\frac14\right)^{9/2}
<\frac1{6000}.
\]

Define normalized instantaneous resolved enstrophy

\[
\mathfrak d_V=N^{-1}\|\nabla V\|_2^2.
\]

The same affine factor gives

\[
\boxed{
N^{-4}\|S_\nu\|
\le \frac{\nu s}{5000}\sqrt{\mathfrak d_V}.
}
\]

Hence on `s<=s0`, a normalized viscous source `rho_nu` forces

\[
\boxed{
\mathfrak d_V
\ge
\left(\frac{5000\rho_\nu}{\nu s_0}\right)^2.
}
\]

Persistence of this event through a positive fraction of one `N^-2` lifetime immediately pays normalized energy dissipation. Temporal concentration into a much shorter set is a separate CKN/burst branch rather than a free escape.

## 3. Scope

This closes the **source-to-physical-observable** step. It does not yet assert that every cubic increment charge contains one fresh Gaussian packet. The next theorem is an Onsager/Littlewood--Paley collision: cubic increments must resolve into a low/base critical-mass band, high-frequency normalized enstrophy, or packet collision entropy.
