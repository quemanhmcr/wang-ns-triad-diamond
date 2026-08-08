# Pressure reservoir synchronization: low-low pairs cannot service infinitely many generations

Band-limited source sampling routes the resolved pressure-third near field to the low-pass quadratic pressure reservoir plus the SGS stress.  The SGS part is already controlled by cubic increments.  This note synchronizes the remaining low-low pressure reservoir through generations.

Let `V_a,V_b` be two low-frequency packet reservoirs with frequencies `M_a,M_b`.  Riesz boundedness and band Bernstein give schematically

\[
\|\nabla^3P_{ab}\|_{3/2}
\lesssim
M_{max}^3\|V_a\|_3\|V_b\|_3,
\qquad
M_{max}=\max(M_a,M_b).
\]

With critical masses `mu_i=M_i E_i`,

\[
\|V_i\|_3\lesssim\sqrt{\mu_i},
\]

so after normalizing by the high block frequency `N`, the pair service coefficient before the final `3/2` source power is

\[
\boxed{
\mathcal P_{ab}(N)
\sim
{M_{max}^3\sqrt{M_aM_b}\over N^3}
\sqrt{E_aE_b}.
}
\]

On a materially reused low-strain lineage, each reservoir frequency grows by less than `21/20`, while the signed-good block scale advances by more than `8/5`.  Hence the coefficient per unit `sqrt(E_aE_b)` contracts by

\[
\boxed{
\left({21\over20}\right)^4
\left({5\over8}\right)^3
={194481\over655360}
<\frac13.
}
\]

Thus even if both reservoirs are adversarially allowed the entire global energy cap at every future service time, one fixed pair has

\[
\boxed{
\sum_{q\ge0}\mathcal P_{ab,q}
<\frac32\mathcal P_{ab,0}^{cap}.
}
\]

A persistent pressure-third near-field branch must therefore do one of the following:

1. relink to new low-frequency reservoir packets;
2. distribute service among many reservoir pairs, producing an atomic/component collision-entropy problem;
3. leave the low-strain material-reuse branch;
4. or use the SGS part of the pressure source, already routed through increments/dissipation.

This is the quadratic pressure analogue of the one-reservoir spectral half-life in `ancestor_reservoir_sync.md`.
