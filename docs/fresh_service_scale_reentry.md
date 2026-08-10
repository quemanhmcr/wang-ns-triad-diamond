# Fresh coherent service as a refinement-invariant scale law

The fresh/new-new branch of coherent increment service currently uses a quarter-dominant **coherent cell** to expose a critical-mass cluster, with cell entropy/cycle as the alternative.  That is a useful ancestry refinement, but it is not the natural state variable for **renewal entrance**: a coherent-cell partition may be refined without changing the physical service law.

The correct quotient is to classify material ownership pointwise on the positive Moyal edge measure and then push the fresh part only to the canonical frequency-band index.

## 1. Freshness is a measurable edge property, not a cell name

Let `O` be the transported measurable old material set in intrinsic `zeta` space.  A positive coherent increment edge has two actual endpoints `zeta_0,zeta_1`.  Define

\[
\chi_{NN}=1_{O^c}(\zeta_0)1_{O^c}(\zeta_1).
\]

For the fixed LP band `j`, let `nu_j` be the positive Moyal increment-service measure and put

\[
F_j=\int \chi_{NN}\,d\nu_j.
\]

This definition contains no coherent cell partition.  Any measurable partition used to represent the phase-space integral merely subdivides the same positive measure, so summing all subcells at fixed `j` gives the same `F_j`.  The sequence `(F_j)_{j\le0}` is therefore the **cell-refinement quotient** of fresh service.

Fresh material provenance is still only an edge mark.  Nothing here declares a whole Fourier shell new material.

## 2. Fix the missing upper support of the canonical LP frame

The existing coherent-service theorem already fixes one smooth square-normalized dyadic analysis--synthesis frame.  Make its annular support explicit once and for all.

Choose a real nonnegative smooth bump `q` supported in `(1/2,2)` and positive on `[3/4,3/2]`.  For `M_j=2^jN`, set

\[
q_j(\xi)=q(|\xi|/M_j),
\qquad
\phi_j(\xi)=\frac{q_j(\xi)}{(\sum_kq_k(\xi)^2)^{1/2}}.
\]

The dyadic positive cores cover every nonzero frequency.  Thus

\[
\sum_j\phi_j^2=1,
\qquad
|\phi_j|\le1,
\qquad
\operatorname{supp}\phi_j\subset\{M_j/2<|\xi|<2M_j\},
\]

and with `u_j=phi_j(D)u`,

\[
u=\sum_j\phi_j(D)u_j
\]

exactly in `L^2`.  This is a Calderón analysis--synthesis frame, not a packet decomposition.  It does not alter the finite `L^3` LP/Bernstein constants already used upstream.

## 3. One smooth band touches only two hard shells

Split the annular support into the exact hard regions

\[
A_0=\{M/2<|\xi|\le M\},
\qquad
A_1=\{M<|\xi|\le2M\}.
\]

Let

\[
\mu_0=M\|P_{A_0}u\|_2^2,
\qquad
\mu_1=2M\|P_{A_1}u\|_2^2.
\]

Since `|phi_j|<=1`,

\[
M\|u_j\|_2^2
\le
M\|P_{A_0}u\|_2^2+M\|P_{A_1}u\|_2^2
=
\mu_0+\frac12\mu_1
\le
\frac32\max(\mu_0,\mu_1).
\]

Therefore

\[
\boxed{
\max(\mu_0,\mu_1)
\ge
\frac23 M\|u_j\|_2^2.
}
\]

This is the exact smooth-band to hard-shell bridge that the cell-based route did not need to state explicitly.

## 4. Fresh band service forces a hard-shell seed

For arbitrary displacement `r`, Moyal and the elementary increment bound give

\[
F_j
\le
M_j\int\|\delta_r u_j\|_2^2d\tau
\le
4\int M_j\|u_j\|_2^2d\tau.
\]

Integrating in the parent scaled time `d tau_N` over an interval of length `c`, some physical time has

\[
M_j\|u_j\|_2^2\ge\frac{F_j}{4c}.
\]

By the two-hard-shell bridge, one actual hard shell of `u` then satisfies

\[
\boxed{
\mu_{hard}\ge\frac{F_j}{6c}.
}
\]

No coherent cell is selected.

## 5. The native state is scale concentration, not a quarter threshold

Assume the already-certified **integrated** fresh branch lower on a parent scaled interval of length `c`,

\[
F:=\sum_{j\le0}F_j\ge Y/4.
\]

Normalize

\[
p_j=F_j/F.
\]

A countable positive probability law has an attained maximal atom.  Put

\[
p_{max}=\max_jp_j,
\qquad
H_\infty^{scale}=-\log p_{max}.
\]

Choosing the maximal **band**, not a coherent cell, gives

\[
\boxed{
\mu_{hard}
\ge
\frac{p_{max}Y}{24c}.
}
\]

Equivalently,

\[
\boxed{
\mu_{hard}e^{H_\infty^{scale}}
\ge
\frac{Y}{24c}.
}
\]

If one wants a collision-entropy coordinate, define

\[
H_2^{scale}=-\log\sum_jp_j^2.
\]

Since `p_max>=sum p_j^2`,

\[
\boxed{
\mu_{hard}e^{H_2^{scale}}
\ge
\frac{Y}{24c}.
}
\]

The `H_2` statement is a corollary; the native theorem is the exact max-scale concentration relation.  Neither quantity is a child-energy causal probability or a new stop class.

## 6. Generic shell reentry, with materiality still deferred

The hard shell enters the existing generic critical-shell theorem.  Its entrance is deliberately material-free.  The fresh NN service which exposed the scale remains provenance only; it is not used to declare the whole `u` shell fresh.

Starting from that event, the existing generic shell theorem inspects its own natural interval `c M^{-2}` (with the same canonical dimensionless block parameter `c` and its existing observed-history guard).  If the shell survives that full no-hit natural corridor, the own-scale service lower is linear in the shell mass.  Hence the same `H_infinity` weighting survives composition into renewed bounded service and integrated service.  Only after that new positive service exists is material OO/ON/NN reread from its actual endpoints.

For `j<=0`, the two candidate hard frequencies are `M_j` and `2M_j<=2N`.  This theorem does not manufacture signed-good scale progress relative to the previous block; that geometry remains a separate supplier property.

## 7. Architectural consequence

The fresh coherent service route can now be separated into two layers:

\[
\text{physical fresh NN service}
\to
\text{canonical scale pushforward}
\to
\text{hard critical shell},
\]

while coherent-cell dominance/entropy/cycle remains optional fine ancestry accounting.  A cell partition is no longer required to create a renewal seed.

No additive reset, frozen packet, per-cell mass floor, or Navier--Stokes global-regularity conclusion is asserted.
