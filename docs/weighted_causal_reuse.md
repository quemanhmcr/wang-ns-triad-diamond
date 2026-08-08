# Transfer-weighted causal reuse information

The raw layered ancestry theorem forces parent-slot reuse by counting distinct causal ancestors.  To connect this to physical transfer, one must not make the opposite mistake of charging the mere existence of the two parent **roles**: an exact single triad always has two parents and can still be an extremal, cost-free interaction.

The correct weighted quantity subtracts that role baseline exactly.

## One causal layer

Let `w(c)` be the normalized physical causal-transfer mass carried by child packet `c`.  Introduce the two parent-role slots `s=1,2` with joint probability

\[
\boxed{p(c,s)=w(c)/2.}
\]

Let `P=f(c,s)` be the physical parent packet/coherent-cell label of that slot, and let `w_parent` be the pushforward law.  Define

\[
\boxed{
\mathcal R
=H(C,S\mid P)
=H(w_{child})+\log2-H(w_{parent})\ge0.
}
\]

If the slot-to-parent map is injective on transfer support, then `R=0`: the full `log 2` increase is just ordinary fresh two-parent causality.  If different causal slots merge into the same physical parent, information about the child/role source is lost after the parent is known, and `R>0`.  Thus `R` is a genuine **transfer-weighted reuse information**, with the role baseline quotiented.

## Layered telescope

For a synchronized depth-`L` ancestry ending in one terminal packet,

\[
\boxed{
\sum_{j=0}^{L-1}\mathcal R_j
=L\log2-H(w_0).
}
\]

The root entropy is at most the log of the number of distinct root cells.  The causal scale / Moyal energy theorem gives

\[
n_0\le{P E_{global}N_{base}\over\eta}
\left(\frac{25}{24}\right)^L,
\]

so

\[
\boxed{
\sum_j\mathcal R_j
\ge
L\log\frac{48}{25}
-\log\frac{P E_{global}N_{base}}{\eta}.
}
\]

Consequently

\[
\boxed{
L\log\frac{36}{25}>
\log\frac{P E_{global}N_{base}}{\eta}
\Longrightarrow
\max_j\mathcal R_j>\log\frac43.
}
\]

For the clean Moyal / affine-grain constants `P=1`, `eta=1/5`, the right-hand budget is `log(5 E_global N_base)`.

This is the transfer-weighted analogue of the raw `rho_j>=1/4` theorem.  It resolves the caveat that combinatorial cycles might live on negligible transfer edges: at sufficiently large causal depth, at least one layer loses more than `log(4/3)` nats of the **physical causal transfer law** through parent reuse.

The next local bridge is narrower: convert such conditional reuse information into the already-existing Rényi collision-pair / multicommodity Hodge-resistance / component-entropy currencies.  The baseline `log 2` from the two parent roles must never be charged as a defect.

## Homogeneity-corrected preferred root bound

The root-count estimate above remains correct when every root independently
carries a fixed `N E>=eta`.  The preferred continuum-facing route no longer
requires that hypothesis.  Young near-extremality fixes shape but is separately
homogeneous in the parent amplitudes, so an absolute root mass cannot be inferred
from shape alone.

Use instead the scale-critical selected coefficient

\[
\alpha=\sqrt N|\langle u,\phi\rangle|.
\]

The companion amplitude--entropy theorem gives on each continuing generated node

\[
\alpha_{p_1}\alpha_{p_2}\ge\Lambda\alpha_c,
\]

and therefore a transfer-weighted lower bound for `E log alpha` at the root.
Bargmann/Moyal gives `N E_anchor>=beta alpha^2`, so log-sum yields

\[
H(w_0)+2E_{w_0}\log\alpha\le\log\sum_r\alpha_r^2.
\]

Substitution into the exact identity

\[
\sum_j\mathcal R_j=L\log2-H(w_0)
\]

retains the same linear `L log(48/25)` growth, with amplitude contributing only
a finite logarithmic offset.  See `amplitude_entropy_causal_reuse.md` and
`common_slice_coefficient_registration.md`.  This is now the preferred physical
root closure.
