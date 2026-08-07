# Curvature-balanced spatial moat

- exact localization law: `E(M)=a/M+b*kappa*M`
- optimizer: `M*=sqrt(a/(b kappa))`
- optimum: `E*=2 sqrt(a b kappa)`
- algebra stress checks: `100000`
- worst optimality residual: `2.994e-16`
- scalar unitary-Fourier Hessian Bernstein constant at lambda=1: `0.085071895494`
- old quadratic schedule, commutator partial sum at 100k: `0.394924`
- old quadratic schedule, curvature partial sum at 100k: `10.590166` (harmonic divergence)
- balanced schedule, each partial sum at 100k: `1.252497`
- balanced schedule tail upper after 100k: `0.006324`

This module records a correction to the previous localization heuristic: an
expanding moat cannot be chosen from commutator considerations alone.  The moat
width must balance filter nonlocality against the curvature of the transported
velocity field.
