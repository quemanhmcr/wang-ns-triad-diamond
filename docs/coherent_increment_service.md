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

For the master-facing construction, fix one smooth square-normalized Littlewood--Paley **analysis--synthesis frame** once and for all.  Put

\[
u_j=\phi_j(D)u,
\qquad
\sum_j|\phi_j(\xi)|^2=1,
\qquad
u=\sum_j\phi_j(D)u_j,
\]

away from zero, with the high analysis shell `j` supported above `|xi|=M_j/2`.  Such a frame is obtained directly from any smooth dyadic annular cover: choose smooth covering bumps `q_j`, then set `phi_j=q_j/(sum_k |q_k|^2)^(1/2)` on the covered frequency region.  This preserves the annular supports and gives the quadratic partition exactly.  The last identity is then the exact Calderón reconstruction and avoids the false requirement `u=sum_j u_j` for overlapping smooth analysis bands.  This does not remove the ordinary finite `L^3` square-function/Bernstein constants `C_LP,C_B` used below.  It does, however, identify the exact `L^2` comparison needed by the later physical tail-energy theorem:

\[
\boxed{
D_{tail}:=N\int\|\nabla P_{>N}u\|_2^2dt
\ge \frac14D_{high}^{LP}.
}
\]

Indeed `|xi|>=M_j/2` implies `M_j^2<=4|xi|^2`, and square normalization gives `sum_j||phi_j(D)u||_2^2` with no overlap loss in `L^2`.  Thus the LP observable and the orthogonal PDE tail remain distinct objects, but their comparison constant is fixed by the same canonical decomposition rather than by a second observer choice.

For these analysis pieces `u_j=phi_j(D)u`, `M_j=2^jN`, standard `L^3` square-function and Bernstein give

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
