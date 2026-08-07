# Base-spinor cross forcing is already the trilinear cross-error ledger

The Hermite/helicity decomposition distinguishes the degree-zero base spinor from
orthogonal sidebands.  The remaining question is whether the degree-zero forcing
caused by the non-Gaussian remainder of the **other two roles** requires a new
PDE norm.  At the level relevant to transfer, it does not.

Let `T` be any trilinear form of operator norm at most one and normalize the
actual roles by

\[
\|f\|=\|g\|=\|h\|=1.
\]

Suppose the one-shot inverse theorem supplies Gaussian roles `F,G,H` with

\[
\|f-F\|,\|g-G\|,\|h-H\|\le\varepsilon.
\]

Then `||F||,||G||,||H||<=1+epsilon`.  Project the nonlinear parent forcing onto
the **base child profile `H`**.  Exact telescoping gives

\[
\begin{aligned}
T(f,g,H)-T(F,G,H)
&=T(f-F,g,H)+T(F,g-G,H),
\end{aligned}
\]

hence

\[
\boxed{
|T(f,g,H)-T(F,G,H)|
\le 2\varepsilon+3\varepsilon^2+\varepsilon^3.
}
\]

This is precisely the work/amplitude of parent-remainder interactions projected
onto the degree-zero child role.  They are **cross trilinear interactions**.
Therefore, once the transfer-preserving component extraction makes omitted
cross interactions `eta_j`, this degree-zero forcing belongs to the existing
logarithmic penalty

\[
\xi_j=\log(1+\eta_j e^{C_j}),
\qquad \Xi=\sum_j\xi_j,
\]

rather than requiring an independent `||F_i||_2` forcing currency.

The child representation mismatch is even simpler:

\[
\boxed{
|T(f,g,h)-T(f,g,H)|\le\varepsilon.
}
\]

Adding the two bounds gives exactly

\[
(2\varepsilon+3\varepsilon^2+\varepsilon^3)+\varepsilon
=3\varepsilon+3\varepsilon^2+\varepsilon^3,
\]

the replacement loss already certified in the one-shot profile module.
At one-percent profile distance the split is

\[
0.020301+0.010000=0.030301.
\]

## Interpretation

This closes an **algebraic identification**, not the full spacetime estimate.
It says:

- degree-zero cross forcing -> existing trilinear cross-error `Xi`;
- `H_1` polarization and `H_3` envelope components -> orthogonal
  sideband/coherence/fresh-grain ledger;
- Gaussian tangent modes -> quotient;
- a large `L^2` norm of the entire nonlinear residual is neither required nor
  implied by this work-level estimate.

The PDE bridge still has to show that the transfer-adapted packet extraction
keeps the relative cross interactions summable through a moving packet lifetime.
