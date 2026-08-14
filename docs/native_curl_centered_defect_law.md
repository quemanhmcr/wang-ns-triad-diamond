# Native curl-centered defect law

Status: **draft structural theorem**, not a global-regularity claim.

The purpose of this note is to record a reduction that sits below the current
cascade/owner/action architecture.  It starts from the Navier--Stokes equation
itself and asks which part of vorticity can actually drive the nonlinearity.
No packet clock, temporal matching rule, owner ledger, or reset budget is used.

## 1. The same intrinsic operator drives nonlinearity and dissipation

On divergence-free fields put

\[
C=\operatorname{curl},\qquad C^2=-\Delta.
\]

The rotational Navier--Stokes equation is

\[
\boxed{
\partial_tu=\mathbb P(u\times Cu)-\nu C^2u.
}
\]

Thus the Euler nonlinearity uses one power of the curl operator while viscosity
uses its square.

For every spatially constant scalar \(\lambda\),

\[
u\times Cu=u\times(C-\lambda)u,
\]

because \(u\times(\lambda u)=0\).  The nonlinear force is therefore blind to
the component of curl parallel to the velocity.

The state itself selects the canonical center

\[
\boxed{
\lambda=\frac{\langle u,Cu\rangle}{\|u\|_2^2}=\frac HE,
}
\]

which is the unique minimizer of \(\|(C-\alpha)u\|_2^2\).  Define the
curl-centered (Beltrami) defect

\[
\boxed{r=(C-\lambda)u},\qquad \langle u,r\rangle=0,
\]

and

\[
\boxed{B=\|r\|_2^2}.
\]

This is an intrinsic state quantity; no frequency threshold or analysis cell
enters its definition.

## 2. One defect, several exact meanings

Let

\[
E=\|u\|_2^2,\qquad H=\langle u,Cu\rangle,\qquad Z=\|Cu\|_2^2.
\]

Orthogonality gives

\[
\boxed{
Z=\lambda^2E+B,
\qquad
B=Z-\frac{H^2}{E}.
}
\]

In the Fourier--helical basis, curl has eigenvalue
\(a=s|k|\).  If modal energy is normalized to a probability measure \(\mu\),
then

\[
\lambda=\mathbb E_\mu[a],
\qquad
\boxed{B=E\,\operatorname{Var}_\mu(a)}.
\]

Hence the same \(B\) is simultaneously

- signed-curl spectral variance;
- squared \(L^2\) distance to the Beltrami eigenspace selected by the state;
- enstrophy in excess of the minimum compatible with \((E,H)\);
- the defect that contains the entire nonlinear Lamb force.

Indeed,

\[
\boxed{
\mathbb P(u\times Cu)=\mathbb P(u\times r).
}
\]

If \(B=0\), then \(Cu=\lambda u\), so the Euler nonlinearity vanishes exactly.
Zero defect is therefore not an externally declared rigidity branch: it is an
exact Beltrami state of the PDE.

## 3. The nonlinear defect is already inside physical energy dissipation

The energy equality is

\[
\frac12E'=-\nu Z.
\]

Using the centered decomposition,

\[
\boxed{
-\frac12E'=\nu\lambda^2E+\nu B.
}
\]

Thus the same defect that is necessary for nonlinear motion contributes
\(\nu B\) directly to actual viscous energy loss.  In particular,

\[
\boxed{
\nu\int_0^T B(t)\,dt\le \frac12E(0)
}
\]

on every interval where the energy equality is available.  This is not a
manufactured reset currency; it is a literal component of Navier--Stokes
energy dissipation.

A basic rigidity estimate is

\[
\boxed{
\|u\times\omega\|_{L^1}
=\|u\times r\|_{L^1}
\le \|u\|_2\,\|r\|_2
=\sqrt{EB}.
}
\]

So vanishing defect turns off the nonlinear Lamb vector in a direct physical
norm.

## 4. Exact defect evolution

Since Euler conserves \(E\) and \(H\), its contribution to \(B'\) is the
nonlinear enstrophy production.  The viscous contribution simplifies because
\(r=(C-\lambda)u\):

\[
\boxed{
(B')_{\rm visc}=-2\nu\|Cr\|_2^2.
}
\]

The nonlinear enstrophy production can also be centered.  Since
\(u\times Cu=u\times r\) and the discarded \(\lambda Cu\) component is
orthogonal to \(u\times r\),

\[
\boxed{
Q=2\langle Cr,u\times r\rangle.
}
\]

Therefore

\[
\boxed{
B'=2\langle Cr,u\times r\rangle-2\nu\|Cr\|_2^2.
}
\]

The active departure from Beltrami rigidity is governed by a native
**defect-stretching minus defect-diffusion** law.

## 5. Closed triads are the spectral shadow: affine null laws and curvature current

For a closed helical triad let \(a_i=s_i|k_i|\) and let \(T_i\) be its three
modal energy works.  The scalar-current law gives

\[
\sum_iT_i=0,
\qquad
\sum_i a_iT_i=0.
\]

Thus every affine observable \(\Phi(a)=c_0+c_1a\) is annihilated by the
nonlinear triad work.

Define the quadratic curvature current

\[
\boxed{
Q_\Delta=\sum_i a_i^2T_i.
}
\]

For pairwise distinct signed frequencies, every scalar observable obeys the
exact interpolation identity

\[
\boxed{
\sum_i\Phi(a_i)T_i
=Q_\Delta\,\Phi[a_0,a_1,a_2],
}
\]

where \(\Phi[a_0,a_1,a_2]\) is the second divided difference.  Coincident
nodes are obtained by the confluent continuous limit.

The nonlinear triad therefore annihilates affine structure and sees only the
curvature of the observable.  If \(\Phi\) is convex, the sign of \(Q_\Delta\)
controls the sign of every such convex response simultaneously.

For the distinguished observable

\[
\Phi(a)=a^2,
\]

the second divided difference is exactly one, hence

\[
\boxed{
\sum_i a_i^2T_i=Q_\Delta.
}
\]

After summing over triads this is precisely nonlinear enstrophy production,
i.e. vortex stretching.  In this language,

\[
\boxed{
\text{vortex stretching}
=\text{signed-curl spectral variance production}.
}
\]

This gives a direct state-space meaning to the earlier same-triad convex
branching law.

## 6. Why the sequence 0,1,2 matters

The two inviscid invariants seen by a closed triad are the moments

\[
\int 1\,d\mu,
\qquad
\int a\,d\mu.
\]

The first strictly convex polynomial moment not fixed by those affine null
laws is

\[
\int a^2\,d\mu,
\]

and the Navier--Stokes Laplacian reads exactly this same symbol because
\(C^2=-\Delta\).  Schematically,

\[
\boxed{
\text{Euler preserves degrees }0,1;
\qquad
\text{viscosity acts through degree }2.
}
\]

The current draft treats this as a structural observation, not yet as a
regularity theorem.  Its significance is that several phenomena previously
handled separately now have one common source:

\[
\text{curl defect}
\longrightarrow
\text{nonlinear Lamb force}
\longrightarrow
\text{spectral variance / vortex stretching}
\longrightarrow
\text{quadratic viscous killing}.
\]

## 7. High-frequency self-control suggested by the identity

For a state whose characteristic signed curl frequency \(|\lambda|\) is large,
small defect means near-Beltrami alignment:

\[
B\ll \lambda^2E
\quad\Longrightarrow\quad
Cu\approx\lambda u,
\]

which suppresses the nonlinear force while the baseline viscous loss
\(\nu\lambda^2E\) is large.  Large defect instead supplies the explicit excess
viscous loss \(\nu B\).

This is not yet a dichotomy theorem for arbitrary localized high-frequency
blocks.  It is a candidate intrinsic self-control principle that must next be
connected to the canonical positive true-upward spectral crossing law without
introducing analyst clocks or temporal matching.

## 8. The next research question

The next target is deliberately not another owner/case closure:

\[
\boxed{
\text{What is the sharp intrinsic conversion law from }
B=\|(C-\lambda)u\|_2^2
\text{ into genuine positive upward spectral transfer?}
}
\]

The strongest possibility is that the previously certified near-extremal
scale geometry \(r_*\approx0.61090410159\) is not merely an extremizer of a
chosen Mellin functional but the genuine geometry that most efficiently
converts Beltrami defect into canonical true-upward Navier--Stokes transport.
That identification is **open** in this draft.

No FIFO/LIFO/proportional temporal matching is used here.  Same-time triad work
is not promoted to between-time stock.  No owner label, packet partition, or
analysis scale is treated as a physical clock.  No global Navier--Stokes
regularity claim is made.
