# Full-natural corridor physical PDE probe

Status: **DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__Q2_CARRIER_ENERGY_AND_ENDPOINT_HARD_SHELL_COVER__SAME_INTERVAL_BOUNDED_INCREMENT_SERVICE**.

The trajectory is the unforced 3D incompressible Navier--Stokes Fourier-Galerkin system with Leray projection, viscosity, 2/3 dealiasing, and RK4 time integration.  It is not a proxy evolution.  The experiment is numerical falsification evidence, not a continuum proof.

Physical corridor: `T=0.0025=c A^-2`, `c=0.01`, `A=2`, `nu=0.05`.

Intrinsic heat-law thresholds: full annular fraction `q=0.309953702966`; radius-3/A retained fraction `q_b=0.154976851483`.

| N | steps | Q2 identity | carrier balance | global balance | shell-cover margin | full heat/carrier | bounded heat lower/carrier | nonlinear work |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 80 | 8.895e-16 | 2.174e-13 | 5.726e-14 | 6.327e-01 | 1.035e+00 | 9.177e-01 | 5.006e-02 |
| 16 | 80 | 5.337e-16 | 2.504e-13 | 6.547e-14 | 6.327e-01 | 1.035e+00 | 9.177e-01 | 5.011e-02 |
| 20 | 80 | 8.894e-16 | 2.492e-13 | 6.654e-14 | 6.327e-01 | 1.035e+00 | 9.177e-01 | 5.011e-02 |

Final carrier-energy resolution spread: `1.449e-07`.

The `Q^2` carrier balance, exact two-hard-shell cover, positive bounded increment service, and physical time integration are all read from the same evolved PDE corridor; none is introduced as an extra event or synthetic resource.
