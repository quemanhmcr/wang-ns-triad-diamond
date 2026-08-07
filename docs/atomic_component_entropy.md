# Atomic-to-component entropy transfer

This module isolates an exact finite-dimensional bridge needed by the grain cascade program.
It is not a Navier--Stokes theorem.

Let atomic companion weights be `w_i`, `sum_i w_i=1`, and let `A(i)` be the ancestry-component label. Put
\[
Q_{\rm at}=\sum_iw_i^2,\qquad Q_{\rm anc}=\sum_A W_A^2,\quad W_A=\sum_{i:A(i)=A}w_i.
\]
Then
\[
\boxed{Q_{\rm anc}-Q_{\rm at}=\sum_A\sum_{i\ne j\in A}w_iw_j.}
\]
Thus the difference is exactly the probability that two independent size-biased atoms are distinct but belong to the same ancestry component.

Writing `q_A=sum_{i in A}(w_i/W_A)^2`,
\[
Q_{\rm at}=\sum_AW_A^2q_A,
\]
so with pair-biased ancestry weights `alpha_A=W_A^2/Q_anc`,
\[
\boxed{H_{\rm at}-H_{\rm anc}=-\log\sum_A\alpha_Aq_A.}
\]
This is an exact collision chain rule. Atomic entropy lost when passing to components becomes internal same-ancestry pair mass; it does not disappear.

For any `0<alpha<1`, if `H_at>=h`, then either
\[
H_{\rm anc}\ge \alpha h
\]
or
\[
\boxed{Q_{\rm anc}-Q_{\rm at}\ge e^{-\alpha h}-e^{-h}.}
\]
The first branch is component/Bellman entropy. The second is pair-cycle mass.

## Ancestry-cycle topology

Take the current active 3-uniform triad hypergraph and its incidence graph. Contract every previously connected ancestry component to one supervertex; keep fresh packet vertices separate. If an old ancestry component contributes `k_A` distinct current packet atoms, contraction raises cycle rank by at least
\[
\boxed{\sum_A(k_A-1).}
\]
Indeed before contraction the current incidence graph has some nonnegative cycle rank; identifying the `k_A` vertices of ancestry class `A` increases cyclomatic number by `k_A-1` (for a connected active component). Thus multiple distinct reused atoms from one old ancestry component create actual cycles in the union of old and new interaction graphs.

## Size-biased effective core

If `Q=sum_iw_i^2` and `lambda>1`, then the mass on atoms satisfying `w_i<=lambda Q` is at least `1-1/lambda`, because for `I` sampled with probability `w_i`, `E[w_I]=Q`. Hence high collision entropy supplies a large-mass core of individually small atoms; if that core is hidden in a small number of ancestry components, it necessarily supplies many repeated ancestry attachments.

The missing analytic step is weighted cycle conversion: pair-cycle mass must be converted either to transfer-weighted Hodge curvature, or in the Hodge-flat case to the spherical/balanced erosion ledger. That is now a flow/congestion problem on the combined ancestry-current incidence graph, rather than an entropy problem.
\n\n## Pair-biased multiplicity theorem\n\nLet\n\[\nd=H_{\rm at}-H_{\rm anc}.\n\]\nUnder the pair-biased ancestry law\n\[\n\alpha_A=\frac{W_A^2}{Q_{\rm anc}},\n\]\nthe exact chain rule says\n\[\n\mathbb E_\alpha q_A=e^{-d}.\n\]\nTherefore for every `lambda>1`, Markov gives\n\[\n\alpha\{A:q_A\le \lambda e^{-d}\}\ge 1-\lambda^{-1}.\n\]\nIf ancestry class `A` contains `k_A` distinct atoms, Cauchy gives `q_A>=1/k_A`. Hence on this pair-biased set\n\[\n\boxed{k_A\ge e^d/\lambda.}\n\]\nThus hidden collision entropy does not merely imply that some reuse occurs. A transfer-relevant ancestry class typically contains exponentially many distinct reused attachments in the hidden entropy gap. Contracting that old ancestry component into one vertex creates at least `k_A-1` cycle-rank gain, before counting cycles already present in the current interaction graph.\n