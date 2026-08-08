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
\n## Five-separated families are Riesz sequences, not merely Bessel families\n\nAt separation `delta=5`, the same packing estimate is much stronger.  Arb certifies\n\n\[\n\boxed{\n\sum_{b\ne a}|G_{ab}|\n<\frac3{50}.\n}\n\]\n\nTherefore the off-diagonal Gram operator `K=G-I` has `||K||<=3/50`, so\n\n\[\n\boxed{\n\frac{47}{50}I\le G\le\frac{53}{50}I.\n}\n\]\n\nFor **arbitrary synthesis coefficients** in this separated equal-covariance family,\n\n\[\n\boxed{\n\frac{47}{50}\sum_a|c_a|^2\n\le\n\left\|\sum_ac_ag_a\right\|_2^2\n\le\n\frac{53}{50}\sum_a|c_a|^2.\n}\n\]\n\nIn particular,\n\n\[\n\boxed{\n\sum_a|c_a|^2\n\le\frac{50}{47}\left\|\sum_ac_ag_a\right\|_2^2.\n}\n\]\n\nThis is the coefficient-energy budget required by the old-reservoir-pool erosion theorem, with the explicit frame constant `P=50/47`, whenever the old reservoir pool is realized as a 5-separated coherent Gaussian synthesis inside one covariance cell.\n\n
## Scope

The Riesz statement does control arbitrary coefficients **after** the family has been reduced to a 5-separated equal-covariance coherent synthesis.  It does not yet prove that an arbitrary transfer-selected Gaussian decomposition admits such a reduction with summable loss, nor synchronize different covariance cells.  The remaining bridge is therefore narrower: construct the iterative transfer-selected reservoir atoms as separated coherent synthesis families inside covariance cells, while rejected close clusters and covariance-cell changes enter `Xi`/fresh/relink bookkeeping.
