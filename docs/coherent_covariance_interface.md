# Coherent covariance interfaces: window changes are a quantified `Xi`, not a new reservoir

The coherent Moyal ledger is independent of the Gaussian window globally, but an iterative packet construction may use different affine covariance representatives at different generations.  This note quantifies the representation change itself so that a small covariance update cannot masquerade as fresh reservoir mass.

For centered normalized Gaussian windows with physical `L^2` covariances `Sigma` and `Theta`, direct Gaussian integration gives

\[
\boxed{
\langle g_\Sigma,g_\Theta\rangle
=2^{3/2}\frac{(\det\Sigma\det\Theta)^{1/4}}{\det(\Sigma+\Theta)^{1/2}}.
}
\]

Let `exp(a_i)` be the eigenvalues of `Sigma^{-1/2} Theta Sigma^{-1/2}`.  Then

\[
\boxed{\langle g_\Sigma,g_\Theta\rangle=\prod_i\cosh(a_i/2)^{-1/2}.}
\]

Because `log cosh x <= x^2/2`,

\[
\boxed{
\|g_\Sigma-g_\Theta\|_2
\le\frac{1}{2\sqrt2}
\left(\sum_i a_i^2\right)^{1/2}.
}
\]

The right side is the affine-invariant SPD log-covariance distance up to the clean constant.

## Exact Moyal stability in the window slot

The coherent transform is linear in the signal and conjugate-linear in the window.  Applying Moyal to the window difference gives exactly

\[
\boxed{
\|\mathcal V_{g_\Sigma}f-\mathcal V_{g_\Theta}f\|_{L^2(d\mu)}
=\|f\|_2\,\|g_\Sigma-g_\Theta\|_2.
}
\]

Therefore

\[
\int\left||\mathcal V_{g_\Sigma}f|^2-|\mathcal V_{g_\Theta}f|^2\right|d\mu
\le2\|f\|_2^2\|g_\Sigma-g_\Theta\|_2,
\]

and hence

\[
\boxed{
\int\left||\mathcal V_{g_\Sigma}f|^2-|\mathcal V_{g_\Theta}f|^2\right|d\mu
\le\frac1{\sqrt2}
 d_{\log}(\Sigma,\Theta)\,\|f\|_2^2.
}
\]

For any common phase-space partition, the sum of absolute changes of the cell energies is bounded by the same quantity.

Thus covariance-cell synchronization has a clean rule:

- common affine transport of center/carrier/cell geometry is exact gauge;
- changing the Gaussian covariance representative by log-distance `delta` costs at most `delta ||f||_2^2/sqrt(2)` in the selected-interface ledger `Xi_cov`;
- a schedule with summable covariance-representative changes has summable representation error;
- a non-small covariance jump is not hidden here and must enter the existing covariance/strain/source/fresh-relink branch.

This separates two effects which must not be conflated: physical covariance dynamics and the harmless choice of nearby coherent analysis window.
