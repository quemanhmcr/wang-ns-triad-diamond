# H1 source episode: from dephasing to grains, entropy, dissipation, or CKN concentration

The H1 no-escape/source stack gives, on the H1-dominant low-strain mild-aspect branch, one fixed physical source channel satisfying

\[
\int_{t_0}^{t_1}\|S_*\|dt\ge\frac{I_1}{132T},
\qquad T=cN^{-2}.
\]

Use scaled time `tau=N^2(t-t0)` and normalized source density

\[
\rho_*(\tau)=N^{-4}\|S_*(t)\|.
\]

Then

\[
\boxed{
\int_0^c\rho_*(\tau)d\tau
\ge\frac{I_1}{132c}.
}
\]

## 1. Persistence versus temporal concentration

Set

\[
\boxed{
\rho_0=\frac{I_1}{264c^2}.
}
\]

The complement of `E={rho_*>=rho_0}` can carry at most half the required source integral, so

\[
\boxed{
\int_E\rho_*d\tau\ge\frac{I_1}{264c}.
}
\]

Fix any desired scaled-time persistence threshold `m_0>0`.

- If `|E|<m_0`, at least half the H1 source is concentrated on a time set of size `<m_0`: this is a genuine temporal-intermittency / CKN-burst branch.
- If `|E|>=m_0`, source collision is persistent enough to pay a spacetime ledger.

No time-average is used to hide a short violent event.

## 2. Differentiated SGS source

On the persistent set, first split according to the affine scale radius `s=Nr_g`.

If `s>s_0` on at least half the set, the affine critical-grain theorem gives a radius-energy event

\[
N\int_{E_2}|u|^2\ge\frac3{10}s_0.
\]

Otherwise a scale-matched subset of measure at least `m_0/2` remains. There the filtered-SGS source theorem turns `rho_0` into a fixed cubic increment threshold. The Onsager collision theorem then says that at every such time either a low/base dyadic band has critical mass above its explicit threshold or high-frequency normalized enstrophy does.

A second pigeonhole gives a subset of measure at least `m_0/4` carrying one of those alternatives.

For the mass branch choose packet fraction `theta=1/4`. Then exactly

\[
\boxed{
\mu_{atom}\ge\frac14\mu_{band}
}
\]

or `H_at>=log4`. Applying the existing atomic/component theorem with `alpha=1/2` gives

\[
\boxed{H_{anc}\ge\log2}
\]

or

\[
\boxed{Q_{anc}-Q_{at}\ge\frac14.}
\]

Thus spatial fragmentation of the increment reservoir is a uniform Bellman/cycle event.

For the enstrophy branch, if its pointwise threshold is `d_0`, normalized high-frequency dissipation pays

\[
\boxed{
N\int\|\nabla P_{>N}u\|_2^2dt
\ge\frac{m_0}{4}d_0.
}
\]

## 3. Viscous source

The same radius split applies. On the scale-matched subset the viscous source collision forces resolved normalized enstrophy `d_V>=d_{nu,0}`. Hence persistence gives

\[
\boxed{
N\int\|\nabla V\|_2^2dt
\ge\frac{m_0}{2}d_{\nu,0}.
}
\]

The factor is `1/2`, not `1/4`, because there is no second mass/enstrophy split after the radius branch.

## 4. Master-facing meaning

After this theorem, the H1 dephasing source branch has the following physical exits:

\[
\boxed{
\begin{array}{c}
\text{pressure-third source},\\
\text{large affine radius / ancestry},\\
\text{dominant low/base packet},\\
\text{Bellman entropy or ancestry cycle},\\
\text{high-frequency / resolved dissipation},\\
\text{temporal CKN concentration}.
\end{array}}
\]

The pressure-third branch retains its separate near/far multipole treatment. The remaining continuum task is no longer to discover a source norm; it is to synchronize the dominant coarse/base packet with the selected spacetime ancestry and to convert the temporal-concentration branch into a CKN-compatible positive cost.
