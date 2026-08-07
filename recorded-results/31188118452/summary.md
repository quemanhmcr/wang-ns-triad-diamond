# Divergence-free localized affine packet test equation

A localized role is tested with `Psi=M_N(chi phi)`, where `M_N` is a smooth shell-localized Leray/helical multiplier.  Since `div Psi=0`, pressure vanishes exactly in the weak full-NS packet coefficient equation:

`d<u,Psi>/dt=<u,partial_t Psi>+<u tensor u,grad Psi>-nu<grad u,grad Psi>`.

Because the shell multiplier is smooth away from zero, its kernel has finite first moment and
`||[chi,M_N]f||_2 <= m1(K) N^-1 ||grad chi||_inf ||f||_2`.
The affine shell bound gives the clean `O(1/M)` estimate
`<= (3/2)m1(K)Cchi M^-1 ||f||_2`.
Thus enforcing divergence-free/helical localization does not create a separate pressure force or aspect penalty; it renormalizes the existing moat commutator coefficient.

Stress checks: `50000`
- worst Leray pressure residual: `1.791e-15`
- worst discrete commutator/bound ratio: `0.996773131`
- worst weak-identity residual: `8.882e-16`
