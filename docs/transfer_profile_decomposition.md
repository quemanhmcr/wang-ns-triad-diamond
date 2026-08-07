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

Assume `S=s`. Let

\[
Q_X=\sum_CX_C^2,\qquad Q_Y=\sum_CY_C^2,\qquad Q_Z=\sum_CZ_C^2.
\]

Hölder gives

\[
s\le(Q_XQ_YQ_Z)^{1/3}.
\]

Since every `Q` is at most one, each one is at least `s^3`. If `C_*` maximizes `X_C`, then

\[
X_{C_*}\ge Q_X\ge s^3.
\]

The key stability identity is

\[
\frac{x^2+y^2+z^2}{3}-(xyz)^{2/3}
\ge
\frac{(x-y)^2+(y-z)^2+(z-x)^2}{6},
\]

because the difference between the two sides is `(xy+yz+zx-3(xyz)^{2/3})/3 >= 0`. Summing over components and using `S>=1-epsilon` yields

\[
\|X-Y\|_2^2+\|Y-Z\|_2^2+\|Z-X\|_2^2\le6\epsilon.
\]

Therefore the same component dominates all three sides:

\[
X_{C_*}\ge(1-\epsilon)^3,
\]

\[
Y_{C_*},Z_{C_*}
\ge
(1-\epsilon)^3-\sqrt{6\epsilon}.
\]

This is an exact finite-dimensional lemma. Near-maximal transfer selects one common interaction component before any Gaussian classification is invoked. The `O(sqrt(epsilon))` loss comes from the stability of AM--GM/Hölder, not from a compactness argument.

## Approximate components

If discarded cross-component interactions have total magnitude at most `eta M_* A_3`, replace `s` by `1-epsilon-eta`. Refining a phase-space partition then gives a nested sequence of dominant components unless a definite transfer deficit appears. This is the transfer-adapted analogue of profile extraction.

## Greedy finite profile extraction

At tolerance `eta`, repeatedly extract an interaction component carrying transfer at least `eta`. The dominant-component certificate makes every extracted profile consume a definite amount `c(eta)` of `L^{3/2}` p-mass. Therefore only finitely many profiles can be extracted. The remainder is small in the nonlinear transfer norm, which is the property needed by the cascade argument; smallness in a stronger ambient norm is unnecessary.

The continuum theorem still requires a wave-packet or wavelet discretization with summable cross-interaction errors. Existing Banach-space profile decompositions provide a compactness framework, but the interaction-component selection above is tailored to the Navier--Stokes trilinear form.
