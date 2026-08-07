# Localized relative-polarization packet bridge

The exact relative-polarization transport theorem identifies the correct
observable, but a PDE packet still has finite frequency width and finite spatial
radius.  This note shows that, on the low-strain branch, those effects enter the
**existing** frequency-cell and curvature ledgers rather than creating a new
uncontrolled polarization error.

For a unit Kelvin direction `a`,

\[
\dot a=-(I-aa^T)A^Ta.
\]

For two gradients `A,B` and unit directions `a,b`, direct differentiation gives

\[
\boxed{
\|f_A(a)-f_B(b)\|
\le 4L\|a-b\|+\|A-B\|_{op},
\qquad L=\max(\|A\|_{op},\|B\|_{op}).
}
\]

At frequency `N`, suppose on one moving window

\[
\|A\|_{op}\le \sigma_0N^2,
\qquad
\|A_x-A_0\|_{op}\le \kappa M N^2,
\]

and the initial normalized carrier cell has directional diameter at most `h`.
For a packet lifetime `T=cN^{-2}`, Gronwall yields

\[
\boxed{
\delta_{dir}(T)
\le e^{4c\sigma_0}(h+c\kappa M).
}
\]

On the signed-good core the parent angle obeys `sin(theta)>9/10`.  Comparing two
triad-normal real frames directly gives the conservative endpoint bound

\[
\|E-E'\|_F\le 8\delta_{dir}
\]

when each of the three carrier directions differs by at most `delta_dir`.
Since trace-free projection is an orthogonal projection in Frobenius norm,

\[
\|(E^TSE)^0-(E'^TS'E')^0\|_F
\le
\|S-S'\|_F+2\|S'\|_{op}\|E-E'\|_F.
\]

Therefore the transfer-relevant generator pair satisfies

\[
\boxed{
\sqrt{
\|\Delta(D_1-D_2)\|_F^2+\|\Delta D_3\|_F^2}
\le
\sqrt5\big(\Delta S_F+16\sigma N^2\delta_{dir}\big).
}
\]

If the spatial strain variation obeys

\[
\Delta S_F\le \kappa M N^2,
\]

the exact capacity-weighted polarization transport identity gives the normalized
additional forcing over one lifetime

\[
E_{pol}
\le
2c\sqrt5\left[
\kappa M+16\sigma_0e^{4c\sigma_0}(h+c\kappa M)
\right].
\]

On the existing low-strain lifetime branch

\[
c\sigma_0\le\frac1{30},
\]

Arb verifies

\[
e^{2/15}<\frac65,
\qquad
\sqrt5<\frac94.
\]

Hence the clean conservative form is

\[
\boxed{
E_{pol}\le 3h+\frac{15}{2}c\,\kappa M.
}
\]

This is the key structural point: helical-frame localization adds no new spatial
currency.  Adding it to the previous window/filter ledger gives

\[
\boxed{
E_{total}(M)
\le
\frac aM+\left(b+\frac{15}{2}c\right)\kappa M+3h.
}
\]

Thus

\[
\boxed{
M_*=\sqrt{\frac{a}{(b+15c/2)\kappa}},
}
\]

and

\[
\boxed{
E_{total,*}
\le
3h+2\sqrt{a(b+15c/2)\kappa}.
}
\]

The `h` term is handled by the already-summable frequency-cell schedule.  The
`kappa M` term is absorbed into the already-required curvature-balanced spatial
moat.  If `c sigma_0` is not small, the existing shape/relative-polarization
strain observability and objective-strain source ledger, rather than this
low-strain perturbative branch, must be used.
