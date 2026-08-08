# Coherent localization operators: canonical synthesis with no reconstruction error

The divergence-free coherent Parseval frame removes compact spatial-window interfaces, but one could still worry that using continuous Moyal cells requires a later arbitrary discrete Gaussian synthesis. It does not.

## 1. Positive coherent localization operators

Let

\[
\phi^\sigma_{a,z}=\mathbb P(g_ze_a)
\]

be the divergence-free coherent Parseval frame. For a measurable phase-space cell \(C\), define weakly

\[
\boxed{
A_C
=
\sum_{a=1}^3
\int_C
|\phi^\sigma_{a,z}\rangle
\langle\phi^\sigma_{a,z}|\,d\mu(z).
}
\]

Every \(A_C\) is positive. For a measurable partition \(\{C\}\), Parseval gives

\[
\boxed{
\sum_C A_C=I
}
\]

on the divergence-free subspace. Since the complementary operator is also positive,

\[
\boxed{0\le A_C\le I.}
\]

## 2. Cell energy controls the actual synthesized piece

Define the positive Moyal cell energy

\[
E_C(f)=\langle f,A_Cf\rangle.
\]

Because \(0\le A_C\le I\), functional calculus gives \(A_C^2\le A_C\). Therefore

\[
\boxed{
\|A_Cf\|_2^2
=\langle f,A_C^2f\rangle
\le E_C(f).
}
\]

Summing over a partition,

\[
\boxed{
\sum_C\|A_Cf\|_2^2
\le
\sum_CE_C(f)
=
\|f\|_2^2.
}
\]

Thus the **canonical coherent synthesis pieces already have Bessel budget \(P=1\)**. The earlier 5-separated Riesz theorem remains useful when one wants literal Gaussian representatives, but it is not required to control the canonical continuous-cell coefficient budget.

## 3. Exact trilinear expansion

For a finite coherent partition and any continuous trilinear form \(T\),

\[
f=\sum_CA_Cf,
\qquad
g=\sum_DA_Dg,
\qquad\ h=\sum_EA_Eh.
\]

Hence exactly

\[
\boxed{
T(f,g,h)
=
\sum_{C,D,E}
T(A_Cf,A_Dg,A_Eh).
}
\]

For countable nested partitions the same identity follows by strong convergence whenever the trilinear form is continuous in the role spaces used by the block theorem.

There is therefore **no continuous-to-discrete reconstruction error** to place in \(\Xi\). The coherent cells themselves give an exact synthesis.

## 4. What remains a genuine interface

Suppose the transfer-selected ancestry keeps only a subset of cell triples. The omitted terms

\[
T(A_Cf,A_Dg,A_Eh)
\]

are actual physical cross-cell interactions. Their sum is a genuine `Xi`/Bellman interface, not a frame error.

The remaining analytic theorem is consequently very specific: use the Gaussian coherent triad kernel and the signed-good resonance/moat geometry to show that these **actual omitted cross-cell interactions** admit the already desired summable defect-space schedule, or else create a Bellman/fresh/reuse event.

This is substantially narrower than constructing an arbitrary Gaussian packet synthesis with controlled coefficients.
