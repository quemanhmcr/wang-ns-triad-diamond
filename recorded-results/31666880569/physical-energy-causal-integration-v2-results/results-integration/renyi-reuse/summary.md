# Renyi causal reuse: weighted binary ancestry reaches existing Bellman/cycle currencies

Status: **EXACT_RENYI_CAUSAL_REUSE_TO_EXISTING_CURRENCIES**.

For a causal child law `w`, duplicating each event into its two parent-role slots gives collision probability `Q_slot=Q_child/2`.  Pushforward through the physical parent label map has

`Q_parent = Q_child/2 + R_hidden`,

where `R_hidden` is exactly the weighted hidden pair mass of distinct causal slots sharing a parent.  Define

`theta = 2 R_hidden / Q_child`.

Then

`Q_parent=(Q_child/2)(1+theta)`.

Across a one-terminal depth-`L` ancestry,

`sum_j log(1+theta_j)=L log2+log Q_root`.

Using `Q_root>=1/n_0` and the coherent root-energy / signed-good scale bound gives

`sum_j log(1+theta_j) >= L log(48/25)-log(P E_global N_base/eta)`.

Therefore the clean depth condition `L log(36/25)>log(P E_global N_base/eta)` forces some layer with

`theta_j>1/3`.

At that layer there is a uniform existing-currency route:

- if `H2_child<log2`, then `Q_child>1/2` and `R_hidden>1/12`: a physical transfer-weighted parent-slot reuse pair;
- if `H2_child>=log2`, apply the existing atomic-to-ancestry collision chain.  Either ancestry/component collision entropy is at least `(1/2)log2`, or child same-ancestry hidden pair mass exceeds `1/sqrt(2)-1/2>1/5`.

Thus the baseline two parent roles remain free, but any sufficiently deep causal ancestry must eventually pay a **transfer-weighted** Bellman/component entropy or pair/cycle currency.  The remaining PDE issue is constructing the synchronized causal layers and registering these parent labels with the same coherent/nested ancestry used by the master.

Stress: `50000`
- worst layer recurrence residual: `2.220e-16`
- worst action telescope residual: `3.553e-15`
- minimum rich-depth action margin: `3.051e-04`
- branches: `{'child_component_Bellman_entropy': 22118, 'weighted_parent_slot_reuse_pair': 8704, 'child_same_ancestry_pair_cycle': 16772, 'not_renyi_reuse_rich': 2406}`
