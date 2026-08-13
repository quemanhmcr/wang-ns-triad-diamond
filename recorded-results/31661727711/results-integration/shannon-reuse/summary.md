# Transfer-weighted causal reuse information

Status: **EXACT_TRANSFER_WEIGHTED_CAUSAL_REUSE_INFORMATION**.

Let `w_(j+1)` be the normalized physical causal-transfer law on child packets of one synchronized layer.  Give each of its two parent-role slots weight `w/2`, and push this joint `(child,role)` law through the physical parent-label map to obtain `w_j`.

Define

`R_j = H(child,role | parent) = H(w_(j+1)) + log 2 - H(w_j) >=0`.

The `log 2` is the exact cost-free two-parent-role baseline.  `R_j`, not `log 2`, measures transfer-weighted causal merging/reuse.

For a one-terminal ancestry the identities telescope exactly:

`sum_j R_j = L log 2 - H(w_0)`.

Since `H(w_0)<=log n_0` and the coherent root-energy / signed-good scale theorem gives `n_0<=(P E_global N_base/eta)(25/24)^L`,

`sum_j R_j >= L log(48/25)-log(P E_global N_base/eta)`.

Therefore the same clean depth condition

`L log(36/25)>log(P E_global N_base/eta)`

forces at least one physical-transfer-weighted layer with

`R_j > log(4/3)`.

This removes the main caveat of the raw cycle-count theorem: the forced reuse signal can be measured with the causal transfer law itself.  What remains is a local conversion of this conditional reuse information into the existing collision-pair / Hodge-resistance / component-entropy currencies; one must not charge the baseline two parent roles as branching loss.

Stress: `50000`
- worst information-telescope residual: `5.329e-15`
- minimum sampled conditional reuse information: `0.000e+00`
- minimum rich-depth margin: `3.555e-05`
- branches: `{'transfer_weighted_reuse_information_rich': 35957, 'causal_depth_not_yet_forced': 14043}`
