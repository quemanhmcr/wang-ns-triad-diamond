# Master finite-dimensional no-escape theorem

This theorem glues the previously proved finite-dimensional modules.  It is an
abstract closure theorem for the grain-cascade model, not a Navier--Stokes
regularity proof.

## 1. Cost-or-potential architecture

For each block `j=0,...,L-1` let `R_j` be its normalized transfer efficiency.
After the transfer-adapted grain extraction, suppose the ideal block estimate
has the form

\[
R_j\le e^{-C_j}+\eta_j,
\qquad C_j\ge0,
\]

where `eta_j` is the omitted cross-component transfer.  Define the exact
cross penalty

\[
\xi_j:=\log(1+\eta_j e^{C_j}).
\]

Then

\[
-\log R_j\ge C_j-\xi_j.
\]

The previous modules supply positive contributions to `C_j` whenever a block
falls into one of the following costly classes:

1. component branching: transfer-weighted Bellman collision entropy;
2. hidden atomic entropy: ancestry/component chain rule followed by Bellman
   entropy or reused-pair cycles;
3. high electrical resistance: Poisson resistance stopping, giving either a
   Bellman cut or conductance collision entropy;
4. curved reused cycles: Hodge energy, converted to transfer deficit by the
   single-edge stability estimate;
5. balanced spherical companions: atomic collision entropy;
6. concentrated fresh companions: a fresh ancestry split; if they are reused,
   they return to cases 2--4.

Fix thresholds so that every costly block obeys

\[
C_j\ge c_0>0.
\]

The only block allowed to have `C_j<c_0` is therefore a synchronized flat,
hemispherical, no-fresh block.

## 2. Spherical potential for the zero-cost branch

For a synchronized flat component contained in the open hemisphere with pole
`n`, define

\[
P_j:=\sup_{x\in S_j}\Phi_n(x),
\qquad
\Phi_n(x)=-\log(n\cdot x).
\]

For an exact optimal midpoint step with no fresh direction, the spherical
midpoint theorem gives

\[
P_{j+1}\le P_j-\kappa_*,
\qquad
\kappa_*=-\log\cos(\theta_*/2)>0.
\]

For a near-extremal block allow a nonnegative perturbation `zeta_j`:

\[
\boxed{P_{j+1}\le P_j-\kappa_0+\zeta_j}
\]

with a fixed `kappa_0>0`.  This is the only permitted low-cost transition.
If a fresh direction is needed to keep the component alive, it leaves this
branch and is charged by the fresh/ancestry ledger.

## 3. Abstract no-escape theorem

Let `F` be the set of low-cost synchronized-flat blocks and `K` its complement.
Assume

\[
P_j\ge0,
\qquad
C_j\ge c_0\quad(j\in K),
\]

and for `j in F`

\[
P_{j+1}\le P_j-\kappa_0+\zeta_j.
\]

Then

\[
|F|\le \frac{P_0+\sum_j\zeta_j}{\kappa_0}.
\]

Indeed, summing the potential decrease over the low-cost transitions and
using nonnegativity of the potential gives this immediately.  Consequently

\[
\sum_{j=0}^{L-1}C_j
\ge
c_0\left(
L-\frac{P_0+\sum_j\zeta_j}{\kappa_0}
\right).
\]

Combining with the exact logarithmic cross-error penalty yields

\[
\boxed{
-\log\prod_{j=0}^{L-1}R_j
\ge
c_0L
-\frac{c_0}{\kappa_0}
\left(P_0+\sum_j\zeta_j\right)
-\sum_j\xi_j.
}
\]

Equivalently,

\[
\boxed{
\prod_{j=0}^{L-1}R_j
\le
\exp\!\left[
\frac{c_0}{\kappa_0}
\left(P_0+\sum_j\zeta_j\right)
+\sum_j\xi_j
\right]e^{-c_0L}.
}
\]

Thus if both perturbation series and the relative cross-error series are
summable, a depth-`L` cascade has exponentially decaying efficiency.

The theorem is deliberately independent of how a costly block pays its cost;
that is what makes the previous Bellman, Hodge, resistance and spherical
modules composable without a proliferating case tree.

## 4. Explicit spherical concentration sublemma

The balanced/concentrated spherical branch can itself be made quantitative.
Let

\[
\mu=\sum_iw_i\delta_{x_i},
\qquad b=\int x\,d\mu(x).
\]

If `|b|<=1-eta`, the exact barycenter--collision estimate gives

\[
H_2(\mu)\ge \log\frac{2}{2-\eta}.
\]

If `|b|>1-eta`, set `n=b/|b|`.  Since

\[
\int(1-n\cdot x)\,d\mu=1-|b|<\eta,
\]

Markov gives, for every angular radius `alpha in (0,pi)`,

\[
\boxed{
\mu\{x:d_{S^2}(x,n)\le\alpha\}
\ge
1-\frac{\eta}{1-\cos\alpha}.
}
\]

Choosing `eta=(1-cos alpha)/4` produces a cap carrying at least `3/4` of the
companion mass.  Hence the spherical exception is quantitatively either
collision-entropy rich or concentrated in a narrow cap.  The latter cap is
then sent to the nested-grain fresh/reuse classification.

## 5. What remains before a PDE theorem

All logical escape variables of the finite-dimensional atomic model now have a
ledger.  The remaining gap is analytic rather than combinatorial: prove that
an arbitrary near-extremal Navier--Stokes transfer block admits the Gaussian
atomic extraction with uniform synthesis constants and summable perturbation
errors, and certify the single-edge stability constant needed to convert Hodge
energy into actual transfer deficit.
