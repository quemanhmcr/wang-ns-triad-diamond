# H1 source episode: source-weighted grains, entropy, and dissipation

The H1 source calculus produces one fixed physical source channel with

\[
\int_{t_0}^{t_1}\|S_*\|dt\ge\frac{I_1}{132T},
\qquad T=cN^{-2}.
\]

With scaled time `tau=N^2(t-t0)` and normalized source density

\[
\rho_*(\tau)=N^{-4}\|S_*(t)\|,
\]

this is the source-weight budget

\[
\boxed{
\Sigma_*:=\int_0^c\rho_*d\tau
\ge\frac{I_1}{132c}.
}
\]

The crucial point is that the downstream collision estimates have the correct homogeneity to use **source weight**, not time measure. Temporal concentration is therefore not a free SGS/viscous escape.

## 1. Differentiated SGS source is linear after the Onsager collision

On the scale-matched branch `s=Nr_g<=s_0`, filtered-source collision gives

\[
Q_{inc}\ge C_0\rho_R^{3/2}.
\]

The Onsager increment theorem uses

\[
X=\left[Q_{inc}/(g_1(C_{LP}C_B)^3)\right]^{2/3}.
\]

Hence the powers cancel exactly:

\[
\boxed{
\mu_{band}\ge c_\mu\rho_R
\quad\text{or}\quad
\mathfrak d_{>N}\ge c_d\rho_R.
}
\]

The coefficients `c_mu,c_d` are explicit functions of the filter, LP/Bernstein constants and the chosen radius cap `s_0`.

This is stronger than a persistence theorem. Partition the **source integral** itself.

- If times with `s>s_0` carry at least `Sigma_*/2`, the episode is a large-radius affine ancestry event, with
  \[
  N\int_{E_2}|u|^2\ge\frac3{10}s_0.
  \]
- Otherwise scale-matched times carry at least `Sigma_*/2`.
- On those times, either the mass or enstrophy branch carries at least `Sigma_*/4` of source weight.

Therefore the high-frequency branch pays directly

\[
\boxed{
D_{>N}:=N\int\|\nabla P_{>N}u\|_2^2dt
\ge\frac{c_d}{4}\Sigma_*.
}
\]

No assumption on the temporal measure of the source set appears.

## 2. The mass branch has a clean packet/ancestry constant

If the mass branch wins, choose `theta=1/4`. At each winning time either

\[
\mu_{atom}\ge\frac14\mu_{band},
\]

or `H_at>=log4`. Applying the exact atomic/component chain rule with `alpha=1/2` yields

\[
\boxed{H_{anc}\ge\log2}
\]

or

\[
\boxed{Q_{anc}-Q_{at}\ge\frac14.}
\]

Thus a winning low/base reservoir becomes one of:

1. a dominant packet, which is classified as fresh or reused;
2. a uniform Bellman entropy event;
3. a uniform same-ancestry cycle event.

Repeated source-weight pigeonholes retain a fixed fraction of `Sigma_*`, so fragmentation cannot be made cheap by shortening its time support.

## 3. Viscous source gets more expensive under temporal concentration

On `s<=s_0`, the viscous source theorem says

\[
\mathfrak d_V
\ge b\rho_\nu^2,
\qquad
b=\left(\frac{5000}{\nu s_0}\right)^2.
\]

If the scale-matched branch carries at least `Sigma_*/2`, then Cauchy on the whole scaled interval `[0,c]` gives

\[
\int\rho_\nu^2d\tau
\ge\frac{\Sigma_*^2}{4c}.
\]

Therefore

\[
\boxed{
D_V:=N\int\|\nabla V\|_2^2dt
\ge
\left(\frac{5000}{\nu s_0}\right)^2
\frac{\Sigma_*^2}{4c}.
}
\]

A shorter, more intermittent viscous source has a larger `L^2_tau` cost, so there is no separate temporal-concentration escape.

## 4. New master-facing source alternative

The H1 dephasing branch is now reduced to

\[
\boxed{
\text{pressure-third}
\ \lor\
\text{large affine radius/reuse}
\ \lor\
\text{dominant low/base packet}
\ \lor\
\text{Bellman entropy/cycles}
\ \lor\
\text{positive normalized dissipation}.
}
\]

The pressure-third channel retains its multipole locality theorem. The main unresolved continuum task is now **synchronization**, not source discovery: identify the dominant low/base packet as the same or a new affine ancestry component, and prevent one old coarse reservoir from servicing infinitely many high-frequency generations without paying holonomy/cycle cost.
