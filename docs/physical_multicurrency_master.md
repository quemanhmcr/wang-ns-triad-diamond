# Physical multi-currency master: conservation budgets only change the finite offset

The finite-dimensional master theorem treats every non-flat block as a multiplicative transfer cost. The PDE architecture contains another legitimate kind of event: a block may spend a **globally finite physical resource** such as fresh coherent energy, viscosity/high-frequency dissipation, or initial-boundary mass. Those events can reset the flat episode without themselves being a pointwise transfer deficit. This note inserts them without losing the master linear rate.

## 1. Single-charge physical resources

Let additive resource classes be indexed by \(r\). A reset assigned to class \(r\) consumes at least \(b_r>0\), while the whole Navier--Stokes trajectory has a total available budget

\[
\sum_{e\in r}B(e)\le B_r<\infty.
\]

Every reset chooses **one primary resource**. It is not charged to several budgets simultaneously. Therefore

\[
\boxed{
N_A\le \sum_r\frac{B_r}{b_r}.
}
\]

Examples of the intended additive side are the already derived fresh-energy/radius budgets, high-frequency or viscous dissipation budgets, and initial-boundary coherent/Fourier budgets. A Bellman/Hodge/Rényi event which already produces multiplicative transfer loss belongs to the transfer-cost class instead; it must not also be counted here.

## 2. Flat erosion with two kinds of reset

Let

- \(N_F\): low-cost flat blocks;
- \(N_T\): multiplicative transfer-cost blocks, each with \(C_j\ge c_0>0\);
- \(N_A\): additive physical-resource resets;
- \(P_{\max}\): the upper reset value of the barycentric potential;
- \(\kappa_0>0\): physical flat erosion rate;
- \(Z=\sum\zeta_j\): total flat-potential perturbation.

Both types of non-flat event may begin a new flat episode. Hence there are at most \(N_T+N_A+1\) flat episodes, and

\[
\boxed{
N_F\kappa_0
\le
(N_T+N_A+1)P_{\max}+Z.
}
\]

Since

\[
L=N_F+N_T+N_A,
\]

one obtains

\[
\boxed{
N_T
\ge
\frac{\kappa_0L-P_{\max}-Z}{\kappa_0+P_{\max}}
-N_A.
}
\]

Combining with the physical resource count,

\[
\boxed{
N_T
\ge
\frac{\kappa_0L-P_{\max}-Z}{\kappa_0+P_{\max}}
-
\sum_r\frac{B_r}{b_r}.
}
\]

Thus finite conservation/dissipation resources can create only a finite number of extra episode resets.

## 3. Multiplicative efficiency remains exponentially expensive

Let the transfer-adapted block estimate be

\[
R_j\le e^{-C_j}+\eta_j,
\qquad
\xi_j=\log(1+\eta_je^{C_j}),
\qquad
\Xi=\sum_j\xi_j.
\]

Each multiplicative costly block has \(C_j\ge c_0\). Therefore

\[
-\log\prod_{j<L}R_j
\ge c_0N_T-\Xi.
\]

Substitution gives

\[
\boxed{
-\log\prod_{j<L}R_j
\ge
c_{\rm eff}L
-
\frac{c_0(P_{\max}+Z)}{\kappa_0+P_{\max}}
-c_0\sum_r\frac{B_r}{b_r}
-\Xi,
}
\]

with exactly the same asymptotic rate as the original master theorem,

\[
\boxed{
c_{\rm eff}
=
\frac{c_0\kappa_0}{\kappa_0+P_{\max}}>0.
}
\]

This is the central point: **finite energy/dissipation resources alter only the finite offset; they do not alter the linear depth coefficient.**

## 4. Physics interpretation

The theorem separates three fundamentally different mechanisms rather than forcing them into one artificial norm:

1. irreversible or finite-capacity physics — energy, enstrophy/dissipation, initial data — gives additive resource budgets;
2. loss of near-extremal transfer coherence — Bellman/Hodge/Rényi/sideband deficits — gives multiplicative transfer cost;
3. an almost ideal Kelvin-flat block gives neither immediately, but consumes barycentric geometric potential at rate \(\kappa_0\).

This mirrors the actual Navier--Stokes structure: conservation supplies finite reservoirs, viscosity supplies a monotone dissipative reservoir, and coherent inviscid transport can only remain cheap while a separate geometric potential continues to erode.

## 5. Scope

This is an exact master algebra once every physical reset has been assigned to one of the two costly classes and the additive classes have proven global budgets. It does **not** silently assert that every PDE source currency already has a global budget. The remaining continuum audit is concrete:

- list each physical branch from service-or-flat rigidity;
- classify it either as a multiplicative transfer cost or as one already-proved globally bounded additive resource;
- ensure every event is primary-charged exactly once;
- insert the now-summable representation/interface \(\Xi\) schedule.

After that audit, the master telescope is no longer missing a different combinatorial mechanism.
