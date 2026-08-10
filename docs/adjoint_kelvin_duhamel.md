# Adjoint Kelvin--Duhamel causal gate

A synchronized causal ancestry must come from the actual Navier--Stokes evolution, not from a static interaction graph.  The natural one-step object is a selected divergence-free packet coefficient in the interaction picture of the **same low-frequency affine/Kelvin transport** already used by the smooth-SGS packet equation.

## 1. Put common physics in the propagator

Along a Kelvin carrier, in the objective transverse frame, write

\[
\dot c=G(t)c+F_{HH}(t)+R_{class}(t),
\]

where

\[
G(t)=-\operatorname{sym}(E^TAE)-\nu|k|^2I.
\]

`G` contains common low--high affine transport and bulk viscosity.  `F_HH` is the actual high--high generation of the selected role.  `R_class` contains only terms which already have a ledger destination: cross low--high cells not absorbed into the common transporter, moving-projector / localization commutators, H1/H3 non-affine sidebands, profile/cross errors, and similar selected-interface terms.

Pressure is **not** added here.  The microscopic role equation is already Leray projected, and the divergence-free localized test packet cancels pressure exactly.

Let the backward dual solve

\[
\boxed{\dot\psi=-G^*\psi.}
\]

Then

\[
\boxed{
\frac d{dt}\langle\psi,c\rangle
=\langle\psi,F_{HH}+R_{class}\rangle.
}
\]

Thus common affine strain and viscosity disappear from the Duhamel source ledger without being declared small.

## 2. Inherit or generate

Over one child slab,

\[
\boxed{z_1=z_0+I_{HH}+I_R.}
\]

Choose the terminal dual along the terminal coefficient, so

\[
|z_1|=\|c(t_1)\|=:A.
\]

The triangle inequality gives the clean alternative

\[
\boxed{|z_0|\ge A/4}
\]

or

\[
\boxed{|I_R|\ge A/4}
\]

or

\[
\boxed{|I_{HH}|\ge A/2.}
\]

The first branch is material inheritance under the common Kelvin propagator.  The second is already-classified residual/interface/source currency.  Only the third branch is called new high--high generation.

## 3. Positive causal generation without persistence

Decompose the generated Duhamel impulse into quadratic atoms

\[
I_{HH}=\sum_\alpha z_\alpha,
\]

where an atom contains the two parent roles at one physical interaction time.  Set

\[
e^{i\vartheta}=I_{HH}/|I_{HH}|.
\]

Then

\[
\sum_\alpha\Re(e^{-i\vartheta}z_\alpha)=|I_{HH}|.
\]

Therefore the positive aligned weights

\[
\boxed{
w_\alpha=[\Re(e^{-i\vartheta}z_\alpha)]_+
}
\]

obey

\[
\boxed{\sum_\alpha w_\alpha\ge|I_{HH}|.}
\]

After normalization they define a positive causal-generation probability law on **same-time quadratic parent pairs**.  No pointwise source persistence is required.

This law measures amplitude generation.  It is not silently identified with physical positive child-energy transfer.  In fact the stronger equality is false even for the scalar flat model `G=R=0`, `c(0)=0`, `c_dot=1`: normalized `dGamma=dt`, while normalized child-energy work is `2t dt`.  The preferred master-facing continuation is therefore the physical-energy causal bridge: keep the Duhamel atoms for same-time parent-pair **support**, but weight the causal layer by actual positive work `2[Re<c,F_HH,alpha>]_+dt`.  The half-slab/parabolic synchronization geometry is unchanged because it depends only on positivity and event support.

## 4. One-step parabolic time geometry

On the signed-good core,

\[
\frac35<\frac{N_p}{N_c}<\frac58.
\]

For the natural parabolic lifetime `T(N)=cN^{-2}`,

\[
\boxed{
\frac{64}{25}<\frac{T_p}{T_c}<\frac{25}{9}.
}
\]

Divide a child slab into two equal temporal halves.  One half carries at least half of the positive aligned generation mass.  If all interactions in that half are assigned their natural parent backward windows, those windows have a common temporal overlap at least

\[
\boxed{
\left(\frac{64}{25}-\frac12\right)T_c
=\frac{103}{50}T_c>2T_c.
}
\]

This last statement is **only natural-window geometry**.  It does not assert coherent packet persistence.  The historical synchronization seam was subsequently removed rather than solved by a packet-persistence theorem: common-slice coefficient registration, asynchronous support extraction, the physical-energy causal bridge, and the moving-role/interface theorems provide the required causal continuation using actual work and smooth roles.  High-tail continuation now has its own scale/time route in `high_tail_binary_work_reentry.md`, `high_tail_ultraviolet_locality.md`, and `high_tail_natural_window_reentry.md`.

Thus this note remains the exact adjoint/Duhamel support gate; it is no longer the current continuum frontier, and normalized Duhamel mass is not promoted to the physical causal probability law.
