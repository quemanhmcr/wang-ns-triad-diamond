# Coherent increment service: from Onsager charge to physical phase-space edges

The previous Onsager collision produced a large low/base critical mass or high-frequency enstrophy, but the low/base mass was still an aggregate band quantity.  This note assigns the **actual increment service** to coherent phase-space cells before choosing any discrete packet synthesis.

Let

\[
Q=\int |G_N(r)|\,\|\delta_r u\|_3^3\,dr,
\qquad g_1=\|G\|_1.
\]

There is a filter displacement `r` for which

\[
\|\delta_r u\|_3^2\ge (Q/g_1)^{2/3}.
\]

For a dyadic decomposition `u=sum_j u_j`, `M_j=2^jN`, standard `L^3` square-function and Bernstein give

\[
\|\delta_ru\|_3^2
\le (C_{LP}C_B)^2\sum_jM_j\|\delta_ru_j\|_2^2.
\]

Put

\[
Y=\frac{(Q/g_1)^{2/3}}{(C_{LP}C_B)^2}.
\]

For `j>=1`, the high-band contribution is at most `2 d_high`, with `d_high=sum_{j>=1}2^j mu_j`.  Hence

\[
\boxed{S_{low}(r):=\sum_{j\le0}M_j\|\delta_ru_j\|_2^2\ge Y-2d_{high}.}
\]

In particular, on `d_high<Y/4`,

\[
\boxed{S_{low}(r)\ge Y/2.}
\]

## Exact coherent service edges

For each low band choose any normalized coherent window and phase-space partition `C`.  Define

\[
\boxed{s_{j,C}(r)=M_j\int_C|\mathcal V_{g_j}\delta_ru_j|^2d\mu.}
\]

Moyal gives

\[
\boxed{\sum_{j\le0,C}s_{j,C}(r)=S_{low}(r).}
\]

Translation covariance is exact:

\[
\boxed{\mathcal V_g(\delta_ru)(X,k)=e^{-ik\cdot r}\mathcal V_gu(X-r,k)-\mathcal V_gu(X,k).}
\]

Therefore

\[
\boxed{s_{j,C}(r)\le2M_j\{E_j(C)+E_j(C-r)\}.}
\]

Thus each service atom is a physical edge between two coherent neighborhoods separated by the actual filter displacement `r`; it is not an arbitrary packet label.

## Old / interface / new no-escape

Classify service edges relative to the transported old reservoir pool.

- old--old edges are bounded by the whole-old-pool spectral service capacity;
- old--new edges are genuine selected-interface/relink edges and belong to `Xi`;
- new--new edges create new coherent ancestry.

Once old-pool erosion gives

\[
C_{old}\le Y/8,
\]

the low-service lower bound yields the clean alternatives

\[
\boxed{d_{high}\ge Y/4}
\]

or

\[
\boxed{\Xi_{cell}\ge Y/8}
\]

or at least `Y/4` new--new coherent service remains.

Normalize the new--new service-edge weights.  With `theta=1/4`, either one edge carries at least one quarter, in which case

\[
\boxed{M_j(E_j(C)+E_j(C-r))\ge Y/32,}
\]

or the service-edge collision entropy is at least `log 4`.  Applying the existing atomic-to-ancestry chain rule with `alpha=1/2` gives

\[
\boxed{H_{ancestry}\ge\log2}
\]

or

\[
\boxed{Q_{ancestry}-Q_{atomic}\ge1/4.}
\]

This is the desired bridge from a cubic SGS/Onsager source to a **positive coherent phase-space service graph**.  The same old-pool half-life that prevented unlimited reservoir reuse now acts on actual increment service; what survives must cross a selected interface or generate new coherent ancestry.
