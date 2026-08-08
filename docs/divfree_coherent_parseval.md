# Divergence-free coherent Parseval frame: remove the compact-window interface from the canonical master

The remaining `Xi` audit contains a difficult spatial moving-window term. A compact material cutoff reduces phase-space overlap to a literal region, but it also creates filter commutators, curvature-balanced moving-boundary errors, and localized pressure work. Those terms are genuine for a compact local-energy argument. They are **not necessary for the canonical coherent ancestry ledger**.

The reason is the orthogonality structure of incompressibility.

## 1. Leray projection does not change coherent coefficients of a divergence-free field

Let \(\mathbb P\) be the Leray projector on vector fields. For a scalar normalized affine coherent state \(g_z\), let

\[
\phi_{a,z}=g_z e_a,
\qquad a=1,2,3.
\]

Define the divergence-free probe

\[
\phi^\sigma_{a,z}=\mathbb P\phi_{a,z}.
\]

If \(u\) is divergence-free, \(\mathbb Pu=u\), and \(\mathbb P\) is self-adjoint. Therefore

\[
\boxed{
\langle u,\phi^\sigma_{a,z}\rangle
=
\langle u,\phi_{a,z}\rangle.
}
\]

The analysis coefficient is exactly the ordinary Gaussian coherent coefficient even though the test function used in the PDE is divergence-free.

## 2. Moyal remains exact on the divergence-free subspace

Componentwise Moyal gives

\[
\sum_{a=1}^3
\int
|\langle u,g_z e_a\rangle|^2d\mu(z)
=
\|u\|_2^2.
\]

Using the coefficient identity above,

\[
\boxed{
\sum_{a=1}^3
\int
|\langle u,\mathbb P(g_ze_a)\rangle|^2d\mu(z)
=
\|u\|_2^2.
}
\]

Thus \(\{\mathbb P(g_ze_a)\}\) is a continuous Parseval frame for the divergence-free subspace. Polarization also gives

\[
\boxed{
\sum_a\int
\langle u,\mathbb P\phi_{a,z}\rangle
\overline{\langle v,\mathbb P\phi_{a,z}\rangle}
\,d\mu(z)
=
\langle u,v\rangle
}
\]

for divergence-free \(u,v\).

This is the correct positive coherent energy/work ledger for incompressible Navier--Stokes.

## 3. Pressure cancels exactly without a compact spatial cutoff

Every projected coherent probe is divergence-free, hence

\[
\boxed{
\langle\mathbb P(g_ze_a),\nabla p\rangle=0.
}
\]

This is global Leray orthogonality, not a local integration-by-parts approximation. Therefore the canonical coherent coefficient equation does **not** acquire pressure boundary work.

Likewise there is no moving compact-window boundary and hence no `a/M+b kappa M` spatial-cutoff error in this canonical analysis.

## 4. What Leray projection changes

The projected probe \(\mathbb P(g_ze_a)\) is not exactly a scalar Gaussian times a constant polarization. In Fourier space it contains the smooth symbol

\[
P(\xi)=I-\frac{\xi\otimes\xi}{|\xi|^2}.
\]

But on the signed-good narrow shell/cap away from \(\xi=0\), this is precisely a smooth degree-zero multiplier. The existing symbol-freezing and helical-frame theorems give a cell error `O(h)`, and the canonical material-label schedule

\[
h_j=h_0 2^{-j}
\]

makes the total representation error summable:

\[
\Xi_{sym}\lesssim\sum_jh_j<\infty.
\]

So the price of using divergence-free coherent probes is **already in the frequency-representation ledger**, not a new spatial interface.

## 5. Consequence for the single `Xi` ledger

For the canonical global ancestry/master route:

- Moyal cell refinement: zero cost;
- common affine/Kelvin motion: zero cost;
- time-slab synchronization: zero cost;
- pressure: zero by Leray orthogonality;
- compact moving-window commutator: absent, because no compact window is introduced;
- Leray/helical variation across one frequency cell: existing summable `Xi_sym`;
- actual material cell switching/cross interaction: physical `Xi`/relink/Bellman currency.

Compact affine windows remain valuable for genuinely local statements such as CKN/local energy diagnostics. This theorem does not invalidate those modules. It says they need not be inserted into the **global coherent ancestry master** merely to define a packet coefficient or cancel pressure.

## 6. Scope

The remaining `Xi` audit is therefore substantially smaller. It concerns actual selected cross interactions, profile/defect-space rejection, and frequency/covariance representation errors. The difficult curvature-balanced compact-window term is not part of the canonical global coherent frame.
