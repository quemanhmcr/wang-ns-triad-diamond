# A transfer-adapted profile decomposition

Let `p=3/2`. Suppose frequency or phase-space grains split into interaction components `C`, with no trilinear interaction between distinct components. Normalize

\[
\sum_C X_C=\sum_CY_C=\sum_CZ_C=1,
\]

where `X_C=||f_C||_p^p`, and similarly for `Y,Z`. Sharp Young on each component gives

\[
|T|\le M_*A_3\sum_C(X_CY_CZ_C)^{2/3}.
\]

Define

\[
S(X,Y,Z)=\sum_C(X_CY_CZ_C)^{2/3}.
\]

By Holder,

\[
S\le \left(\sum X_C^2\sum Y_C^2\sum Z_C^2\right)^{1/3}\le1.
\]

Thus replication has a strict convexity cost. Two equal disconnected copies produce `S=1/2`, not `1`.

## Quantitative dominant-component certificate

Assume `S=s`. Choose `C_*` maximizing `X_C`. Since

\[
s\le(\sum X_C^2)^{1/3},
\]

we have `X_{C_*} >= s^3`. The contribution away from `C_*` is at most

\[
\left(\sum_{C\ne C_*}X_C^2\right)^{1/3}
\le(1-s^3)^{1/3}.
\]

Consequently

\[
Y_{C_*},Z_{C_*}\ge
\left[s-(1-s^3)^{1/3}\right]_+^{3/2}.
\]

This is an exact finite-dimensional lemma. Near-maximal transfer selects one common interaction component before any Gaussian classification is invoked.

## Approximate components

If discarded cross-component interactions have total magnitude at most `eta M_* A_3`, replace `s` by `1-epsilon-eta`. Refining a phase-space partition then gives a nested sequence of dominant components unless a definite transfer deficit appears. This is the transfer-adapted analogue of profile extraction.

## Greedy finite profile extraction

At tolerance `eta`, repeatedly extract an interaction component carrying transfer at least `eta`. The dominant-component certificate makes every extracted profile consume a definite amount `c(eta)` of `L^{3/2}` p-mass. Therefore only finitely many profiles can be extracted. The remainder is small in the nonlinear transfer norm, which is the property needed by the cascade argument; smallness in a stronger ambient norm is unnecessary.

The continuum theorem still requires a wave-packet or wavelet discretization with summable cross-interaction errors. Existing Banach-space profile decompositions provide a compactness framework, but the interaction-component selection above is tailored to the Navier--Stokes trilinear form.
