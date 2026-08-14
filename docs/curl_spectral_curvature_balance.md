# Curl-spectral curvature law: one intrinsic nonlinear grammar behind energy, helicity, critical growth, strain and radial flux

Status: **draft theorem programme**.  This note records exact identities already bound to the repository's closed-helical-triad and canonical cyclic-work laws.  It does **not** claim global Navier--Stokes regularity.

The purpose is to remove case language wherever the PDE itself supplies a smaller law.

## 1. Positive state on the signed curl spectrum

For a helical Fourier mode `(k,s)`, `s in {+1,-1}`, put

\[
\lambda=s|k|,
\qquad
E_{k,s}=|a_{k,s}|^2.
\]

The physical energy state is the positive measure

\[
\rho_t=\sum_{k,s}E_{k,s}(t)\,\delta_{s|k|}.
\]

Its two affine moments are

\[
\int d\rho=\|u\|_2^2,
\qquad
\int \lambda\,d\rho=H(u).
\]

The critical positive moment is

\[
K(t)=\int |\lambda|\,d\rho
=\|u\|_{\dot H^{1/2}}^2.
\]

No hard cell, coherent packet, material name or checkpoint is a state node in this definition.

## 2. A closed triad has one signed nonlinear degree of freedom

Take one actual closed helical triad and order its curl eigenvalues

\[
a<m<b.
\]

Let the three rooted signed energy works be `(T_a,T_m,T_b)`.  The already-certified closed-triad law gives simultaneously

\[
T_a+T_m+T_b=0,
\]

and

\[
aT_a+mT_m+bT_b=0.
\]

Thus the work vector is one-dimensional.  With

\[
q:=-T_m,
\]

we have exactly

\[
T_a=q\frac{b-m}{b-a},
\qquad
T_m=-q,
\qquad
T_b=q\frac{m-a}{b-a}.
\]

`q>0` is a mean-preserving spread of positive energy mass from the median curl eigenvalue to the two extremes.  `q<0` is its reverse contraction.  These are not analyst branches; they are the two orientations of the same one-dimensional physical law.

## 3. The barycentric tent and the curvature equation

Define the triangular tent

\[
\Theta_{a,m,b}(\lambda)=
\begin{cases}
\dfrac{(b-m)(\lambda-a)}{b-a},&a\le \lambda\le m,\\[5pt]
\dfrac{(m-a)(b-\lambda)}{b-a},&m\le \lambda\le b,\\[5pt]
0,&\text{otherwise}.
\end{cases}
\]

Distributionally,

\[
\Theta''
=
\frac{b-m}{b-a}\delta_a
-\delta_m
+\frac{m-a}{b-a}\delta_b.
\]

Therefore one closed triad contributes

\[
\boxed{
\partial_t\rho_\triangle^{NL}
=
\partial_\lambda^2(q\Theta_{a,m,b}).
}
\]

At continuum measure level the existing normalization is retained exactly:

\[
dW=C_F T\,d\Lambda_{\rm closed},
\qquad
C_F=(2\pi)^{-3/2}.
\]

Hence define the signed physical curvature potential

\[
\kappa_t(\lambda)
=
C_F\int q_\triangle(t)\Theta_\triangle(\lambda)\,d\Lambda_{\rm closed}.
\]

For finite/Galerkin NS, and locally on continuum curl intervals where the registered Radon moments exist,

\[
\boxed{
\partial_t\rho
=
\partial_\lambda^2\kappa
-2\nu\lambda^2\rho.
}
\]

This is a state-balance law.  `kappa` is signed and is reconstructed **before Hahn**.  Canonical `dW+` remains the causal positive energy-transfer law.

## 4. Why `kappa` is intrinsic rather than a triad-decomposition artifact

Let

\[
S_{NL}:=\partial_t\rho+2\nu\lambda^2\rho.
\]

Energy and helicity conservation say

\[
\int dS_{NL}=0,
\qquad
\int \lambda\,dS_{NL}=0.
\]

The decaying solution of

\[
\kappa''=S_{NL}
\]

is therefore uniquely reconstructed by

\[
\boxed{
\kappa(\lambda)
=
\frac12\int_{\mathbb R}|\lambda-\mu|\,dS_{NL}(\mu).
}
\]

Thus two valid representations of the same physical nonlinear spectral source cannot produce different curvature potentials.  The triad tents are a disintegration of an intrinsic state-current potential, not its definition by analyst choice.

## 5. One weak law generates the familiar observables

For any sufficiently admissible spectral test function `phi`, integration by parts gives

\[
\boxed{
\frac d{dt}\int\phi(\lambda)\,d\rho
=
\int\phi''(\lambda)\kappa(\lambda)\,d\lambda
-2\nu\int\lambda^2\phi(\lambda)\,d\rho.
}
\]

### Energy and helicity

For `phi=1` and `phi=lambda`, `phi''=0`.  Nonlinearity disappears exactly.  Energy and helicity are therefore the affine kernel of the same curvature operator.

### Enstrophy / vortex stretching

For `phi=lambda^2`, `phi''=2`:

\[
\boxed{
\dot Z^{NL}=2\int\kappa(\lambda)\,d\lambda.
}
\]

Thus the full nonlinear enstrophy/strain reading is the **area** under the same curvature potential.

### Critical `H^{1/2}` production

For `phi=|lambda|`, distributionally

\[
\phi''=2\delta_0.
\]

Hence

\[
\boxed{
\dot K^{NL}=2\kappa(0).
}
\]

Critical growth is the **height of the same curvature potential at the helicity-sign interface**.  Homochiral triads have support on one side of zero, where `|lambda|` is affine, and therefore contribute exactly zero.

## 6. Full-field helicity-pair law

Put

\[
K_+=\int \lambda_+\,d\rho,
\qquad
K_-=\int \lambda_-\,d\rho,
\]

where `lambda_+=max(lambda,0)` and `lambda_-=max(-lambda,0)`.

Since

\[
(\lambda_+)''=(\lambda_-)''=\delta_0,
\]

both positive helicity magnitudes have the same nonlinear source:

\[
\boxed{
\dot K_+^{NL}=\dot K_-^{NL}=\kappa(0).
}
\]

Including viscosity,

\[
\dot K_+ +2\nu D_+^{(3)}=\kappa(0),
\qquad
\dot K_- +2\nu D_-^{(3)}=\kappa(0).
\]

Thus nonlinear critical growth is opposite-helicity **pair creation** in the exact deterministic energy ledger.  It is not a new probability law.

## 7. Radial cascade is the slope of the same potential

Let

\[
j(\lambda):=-\partial_\lambda\kappa(\lambda).
\]

Then the nonlinear balance is an ordinary one-dimensional continuity equation

\[
\partial_t\rho+\partial_\lambda j=0
\]

before viscosity.

For radial exterior stock

\[
E_{>R}=\int_{|\lambda|>R}d\rho,
\]

we get

\[
\boxed{
\dot E_{>R}^{NL}
=
\kappa'(-R)-\kappa'(R).
}
\]

Integrating the radial tail balance over all `R>=0` gives the layer-cake identity

\[
\boxed{
\int_0^\infty
[\Phi_\uparrow(R)-\Phi_\downarrow(R)]\,dR
=2\kappa(0)
=\dot K^{NL}.
}
\]

So critical `H^{1/2}` production, opposite-helicity pair creation and net radial first-moment transport are three readings of one state-current law.

## 8. The Tanaka--scale cocycle: critical production already contains native scale displacement

Consider a critical-growing spread.  Let

\[
Q=-T_m>0
\]

be its actual total physical donor work, let

\[
N_d=|m|,
\qquad
N_c=\max(|a|,|b|).
\]

The two extreme recipients form the barycentric image of the median donor.  Therefore

\[
\frac{\mathcal P_{crit}}{Q}
=
\mathbb E_{rec}|\lambda|-|m|
\le N_c-N_d.
\]

Define

\[
\eta_{crit}:=
\frac{\mathcal P_{crit}}{N_cQ}.
\]

Then purely from positive energy mass plus barycenter conservation,

\[
\boxed{
0<\eta_{crit}
\le
1-\frac{N_d}{N_c}
\le
\log\frac{N_c}{N_d}.
}
\]

If `T_N=cN^{-2}` is the native parabolic lifetime,

\[
\boxed{
\eta_{crit}
\le
\frac12\log\frac{T_d}{T_c}.
}
\]

No Waleffe envelope, phase efficiency, Young constant, log-progress `J`, capacity currency or event count appears.  Critical pair creation and parabolic scale displacement are two observables of the same elementary barycentric move.

## 9. Strong-production sharpening

The exact strict-UV variational problem for native critical production is **not** the existing log-progress `J` problem.

For the majority-helicity sector, with child radius normalized to one, same-helicity median donor ratio `D` and opposite-helicity recipient ratio `S`,

\[
\mathscr C_A(D,S)=2S(1-D)|g_{(+,-,+)}(D,S,1)|.
\]

The unique global interior maximizer satisfies

\[
\boxed{
S_*=4D_*^2,
\qquad
12D_*^3-5D_*^2+2D_*-1=0,
}
\]

with

\[
D_*\approx0.4539303,
\qquad
S_*\approx0.8242110.
\]

At this physical geometry the existing log-progress efficiency is only about

\[
J/J_*\approx0.2715.
\]

Therefore `J`-bad does **not** mean critically weak.  `J` is a scale-progress proof observable; it is not the intrinsic danger functional for the critical state.

A separate exact threshold theorem proves

\[
\boxed{
\mathscr C_A>\frac1{8\sqrt2}
\Longrightarrow
D<\frac58
\Longrightarrow
\frac{T_d}{T_c}>\frac{64}{25}.
}
\]

This is a sharpening of the continuous Tanaka--scale cocycle in the strongly critical-efficient region.

## 10. J-free generated energy-donor telescope

Once actual physical-energy reentry has established positive HH child-energy work, suppose the same event has a proved physical energy donor with

\[
N_d/N_c<5/8.
\]

A heavy half of the child work lies in temporal width

\[
|H|\le T_c/2 < \frac{25}{128}T_d.
\]

With common donor registration surface

\[
s=a_H-\frac25T_d,
\]

successive donor lifetimes expand backward by more than `64/25`.  The upper-ratio-only asynchronous theorem gives

\[
\boxed{
s_j-s_{j+1}\ge\frac{6859}{16000}T_j,
}
\]

and hence

\[
\boxed{
s_0-s_L
\ge
\frac{6859}{24960}T_0
\left[\left(\frac{64}{25}\right)^L-1\right].
}
\]

This theorem needs no lower donor/child ratio and no `J` classification.  The registration surfaces are not events; the conclusion is finite physical backward depth once a genuine generated donor lineage has been identified.

The canonical cyclic donor kernel supplies a mass-preserving same-event pushforward from critical-efficient recipient work to its unique median energy donor.  The full donor-scale distribution is retained; no argmax scale selector is introduced.

## 11. Why this still does not prove regularity

The balance grammar by itself is not enough.  An abstract signed curvature current could move positive energy outward on arbitrarily short time intervals while preserving the first two moments and accumulating only finite viscous exposure.  Such a construction would satisfy the state-balance form but need not be realizable by Navier--Stokes.

Therefore the remaining primitive is **constitutive**, not another owner type:

\[
\boxed{
\text{How strongly can the actual Waleffe quadratic law drive }\kappa
\text{ relative to }2\nu\lambda^2\rho?
}
\]

The intended small grammar is now:

\[
\boxed{
\text{positive energy state}
+\text{barycentric signed curvature current}
+\text{Waleffe constitutive rate}
+\text{viscous killing}
+\text{observer quotient}.
}
\]

Pressure/SGS/strain/material/checkpoint objects may remain genuine faces of derived-state equations or representations.  They do not automatically become new full-field energy generators.  Any claim of full-field energy ancestry must return to the fundamental state/current balance above.

## 12. Scope guards

This draft explicitly does **not** infer:

- global convex monotonicity: contractions reverse the sign of curvature production;
- a finite gross-transfer budget;
- `J`-bad work is harmless;
- a fresh Hahn law from `kappa`;
- that critical production itself is causal probability;
- FIFO/LIFO temporal matching;
- global Navier--Stokes regularity.

The next load-bearing theorem is a constitutive speed/regularity law for the actual NS curvature current, not another case split.
