# Continuum helical edge measure on actual Galerkin Navier--Stokes

A single real divergence-free smooth Fourier polynomial is evolved by the
unforced three-dimensional incompressible Navier--Stokes Galerkin system with
Leray projection, viscosity, 2/3 dealiasing and RK4.  For child mode
`(5, 1, 0)`, every retained ordered convolution parent is read from the
same evolved state, quotiented into unordered parent orbits, and then resolved
into all eight helicity edges.

The torus Fourier-series convolution has coefficient one.  To test the theorem's
unitary-R3 measure normalization without changing the physical Galerkin source,
each discrete unordered orbit is assigned quotient mass `1/C_F`, so the theorem
factor `C_F` cancels exactly against the discrete counting-measure embedding.

| n | cutoff | unordered pairs | helical edges | source residual | work residual | progress residual | NS energy balance |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 5 | 329 | 2632 | 3.601e-16 | 3.429e-16 | 7.502e-16 | 1.305e-09 |
| 24 | 7 | 1049 | 8392 | 3.377e-16 | 2.915e-16 | 6.215e-16 | 1.365e-09 |
| 28 | 8 | 1631 | 13048 | 6.163e-16 | 4.371e-16 | 4.662e-16 | 1.368e-09 |

The native cutoffs above define different Galerkin PDEs. Their final-child-energy spread `6.634e-02` and integrated-child-work spread `7.695e-02` are **diagnostics only**; no unproved convergence threshold is imposed.

The following table embeds the **same** Galerkin cutoff `5` on every FFT grid:

| n | cutoff | unordered pairs | helical edges | source residual | work residual | progress residual | NS energy balance |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 5 | 329 | 2632 | 3.601e-16 | 3.429e-16 | 7.502e-16 | 1.305e-09 |
| 24 | 5 | 329 | 2632 | 3.240e-16 | 2.836e-16 | 6.215e-16 | 1.305e-09 |
| 28 | 5 | 329 | 2632 | 3.340e-16 | 4.571e-16 | 7.216e-16 | 1.305e-09 |

Same-system final-child-energy representation spread: `1.196e-15`.
Same-system integrated-child-work representation spread: `1.487e-15`.

This is direct falsification evidence on finite Fourier--Galerkin Navier--Stokes,
not a continuum PDE proof.  It checks that the proposed signed edge measure is
actually the same nonlinear convolution/work law on evolved PDE states before
Hahn splitting.
