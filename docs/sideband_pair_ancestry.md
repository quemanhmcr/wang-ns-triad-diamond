# Pair-sideband rescue is an ancestry graph, not a residual norm

The odd-Hermite selection rule says that a nonzero triad interaction has even total Hermite parity.  Therefore an interaction containing odd sidebands has either zero or **two** odd roles; three odd roles also vanish.  After suppressing the even/base role, every pair-rescue interaction is an ordinary edge between two odd daughter atoms.

Let `a_e>=0` be rescue capacities and give each odd atom an ancestry label.  Split

\[
W=W_{cross}+W_{same}
\]

according to whether the two odd endpoints have different or equal ancestry labels.  If a fixed fraction of `W` is cross-ancestry, that fraction is already a true cross-component interaction and enters the existing `Xi` ledger.

For the same-ancestry graph define weighted degree

\[
d_v=\sum_{e\ni v}a_e,
\qquad
w_v={d_v\over2W_{same}}.
\]

This is a probability law.  For ancestry class `A`,

\[
W_A=\sum_{v\in A}w_v.
\]

Because every retained edge has both endpoints in the same ancestry class, the factor two cancels exactly and

\[
\boxed{
W_A={1\over W_{same}}\sum_{e\subset A}a_e.
}
\]

Thus `W_A` is not an artificial vertex weighting: it is exactly the pair-rescue edge-mass distribution across ancestry classes.

Define

\[
Q_{at}=\sum_vw_v^2,
\qquad
Q_{anc}=\sum_AW_A^2.
\]

Then the existing atomic/component collision chain rule becomes

\[
\boxed{
Q_{anc}-Q_{at}
=\sum_A\sum_{u\ne v\in A}w_uw_v.
}
\]

For `H_at=-log Q_at` and `H_anc=-log Q_anc`, fix `h>0` and `0<alpha<1`.  If `H_at<h`, then `Q_at>e^{-h}` and since `max_v w_v>=Q_at`, one daughter carries

\[
\boxed{w_v>e^{-h}},
\]

a dominant-reuse branch.  If `H_at>=h`, either `H_anc>=alpha h`, which is component/Bellman collision entropy, or

\[
\boxed{
Q_{anc}-Q_{at}
\ge e^{-\alpha h}-e^{-h}.
}
\]

The latter is distinct same-ancestry endpoint-pair mass.  By the repository's existing contraction theorem, multiple distinct attachments from an old ancestry component create cycle-rank gain in the union of old and current interaction graphs.

Consequently pair-sideband rescue is routed without a disjoint-support assumption:

\[
\boxed{
\text{cross }\Xi
\quad\text{or}\quad
\text{dominant daughter reuse}
\quad\text{or}\quad
\text{Bellman component entropy}
\quad\text{or}\quad
\text{same-ancestry cycles}.
}
\]

The theorem does not by itself price the dominant-reuse or cycle branches; those are handed to the existing scale/spin holonomy and Hodge/erosion machinery.  Its role is to prevent quadratic Hermite rescue from surviving as a new unstructured error term.
