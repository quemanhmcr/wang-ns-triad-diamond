# Affine coherent-state Bessel budget inside one covariance cell

The whole-old-pool erosion theorem needs a frame-energy budget

\[
\sum_aE_a\le P\|u\|_2^2.
\]

For arbitrary redundant packet coefficients this is not automatic.  But it is automatic for **analysis coefficients of separated affine Gaussian coherent states**, which are the natural uncertainty-limited probes of a Gaussian reservoir.

Use the L2-normalized physical packet

\[
g_{X,k,L}(x)=C_L
\exp\!\left[-\frac14|L^{-1}(x-X)|^2\right]e^{ik\cdot x}.
\]

For two packets with the same covariance factor `L`, direct Gaussian integration gives

\[
\boxed{
|\langle g_{X,k,L},g_{Y,\ell,L}\rangle|
=
\exp\!\left[
-\frac18|L^{-1}(X-Y)|^2
-\frac12|L^T(k-\ell)|^2
\right].
}
\]

Define the intrinsic affine phase-space coordinate

\[
\boxed{
\zeta(X,k;L)=
\left(\frac12L^{-1}X,\;L^Tk\right)\in\mathbb R^6.
}
\]

Then

\[
\boxed{
|\langle g_a,g_b\rangle|
=e^{-|\zeta_a-\zeta_b|^2/2}.
}
\]

This coordinate is exactly invariant under a common physical affine change

\[
X'=SX,\qquad k'=S^{-T}k,\qquad L'=SL.
\]

Thus the orthogonality metric is grain intrinsic, not an Euclidean-aspect penalty.

## Six-dimensional packing and Schur

Suppose the `zeta_a` are `delta`-separated.  Around every point put a disjoint ball of radius `delta/2`.  For a fixed center `a`, the number of other centers in the shell

\[
n\delta\le|\zeta_b-\zeta_a|<(n+1)\delta
\]

is at most

\[
(2n+3)^6.
\]

Hence every absolute Gram row satisfies

\[
\sum_b|G_{ab}|
\le
1+
\sum_{n\ge1}(2n+3)^6e^{-n^2\delta^2/2}.
\]

At `delta=4`, Arb certifies

\[
\boxed{
1+\sum_{n\ge1}(2n+3)^6e^{-8n^2}<\frac{25}{4}.
}
\]

Schur's test for the Gram matrix therefore yields the explicit coherent-state Bessel inequality

\[
\boxed{
\sum_a|\langle f,g_a\rangle|^2
\le\frac{25}{4}\|f\|_2^2.
}
\]

This is precisely a packet-frame energy budget with `P=25/4` for a `4`-separated equal-covariance analysis family.

## Scope

This theorem does **not** say arbitrary coefficients in a redundant Gaussian synthesis obey the same bound: large cancelling coefficients can be artificial.  It also does not yet compare transfer-extraction coefficients to coherent analysis coefficients, nor synchronize different covariance cells.  The correct next bridge is therefore narrower: construct the iterative transfer-selected reservoir atoms so that their physical mass is controlled by separated coherent-state analysis coefficients, with covariance-cell changes and rejected interfaces entering `Xi`/fresh/relink bookkeeping.
