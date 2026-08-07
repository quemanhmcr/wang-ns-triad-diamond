# Odd-sideband pair rescue to spacetime ancestry

Parity makes every nonzero rescue interaction containing odd Hermite modes use exactly **two odd daughter endpoints** and one even/base role.  Therefore pair rescue projects canonically to an ordinary weighted graph on odd daughter atoms.

Split rescue mass into `W_cross` (odd endpoints in different ancestry components) and `W_same`.  A cross fraction at least `1/2` is already an existing `Xi`/cross-component branch.

On same-ancestry edges define

`w_v=d_v/(2 W_same)`.

For each ancestry class `A`, `W_A=sum_(v in A) w_v` equals **exactly** the fraction of same-edge rescue mass lying in `A`.  Hence

`Q_anc-Q_at = sum_A sum_(u!=v in A) w_u w_v`.

With thresholds `h>0`, `alpha in (0,1)`, the same-edge graph has the exact trichotomy:

1. `H_at<h` -> some odd daughter carries endpoint weight `> exp(-h)`: dominant reused daughter;
2. `H_at>=h` and `H_anc>=alpha h` -> ancestry/component Bellman collision entropy;
3. `H_at>=h` and `H_anc<alpha h` -> distinct same-ancestry pair mass at least `exp(-alpha h)-exp(-h)`, hence repeated attachments which become cycle-rank gain after contraction by the existing atomic-to-ancestry theorem.

Together with the cross branch, pair-sideband rescue is routed into the existing currencies `Xi`, dominant reuse, Bellman entropy, or ancestry cycles.  No spatial disjointness of Hermite modes is assumed.

Stress: `50000`
- branch counts: `{'cross_Xi': 34486, 'ancestry_Bellman_entropy': 5679, 'same_ancestry_pair_cycle': 9802, 'dominant_reused_daughter': 33}`
- worst ancestry-class mass identity residual: `5.551e-16`
- worst collision-chain identity residual: `9.992e-16`
- minimum branch margin: `1.388e-05`
