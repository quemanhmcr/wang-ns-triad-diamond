# Same-carrier checkpoint segmentation: physical PDE probe

Status: **DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__ONE_EVOLVED_TRAJECTORY_FIXED_Q_AND_TERMINAL_DUAL__ACTUAL_COMPLEX_HH_RESIDUAL_IMPULSES__CHECKPOINT_PARTITION_INVARIANCE**.

The trajectory is the unforced 3D incompressible Navier--Stokes Fourier-Galerkin system with Leray projection, viscosity, 2/3 dealiasing and RK4.  One fixed Q and one terminal complex dual are read on the same evolved trajectory.  This is numerical falsification evidence, not a continuum proof.

Physical interval: `T=0.015625`, `A=4`, `nu=0.05`, fixed natural windows `4`.

| N | steps | div | global balance | Q2 balance | nonlinear split | Duhamel | low-low moat | HH impulse | residual impulse | segmentation |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 80 | 3.919e-17 | 1.600e-11 | 1.963e-07 | 3.071e-16 | 2.664e-09 | 3.559e-18 | 7.936e-04 | 4.106e-04 | 0.000e+00 |
| 24 | 80 | 4.746e-17 | 1.600e-11 | 1.970e-07 | 3.277e-16 | 2.815e-09 | 2.878e-18 | 7.935e-04 | 4.106e-04 | 0.000e+00 |
| 28 | 80 | 5.579e-17 | 1.600e-11 | 1.970e-07 | 3.445e-16 | 2.815e-09 | 1.822e-18 | 7.935e-04 | 4.106e-04 | 0.000e+00 |

Terminal coefficient resolution spread: `1.928e-06`.

All cumulative complex impulses, resolved strain, carrier/global balances and checkpoint partitions are read from the same evolved PDE states. No proxy evolution, segment reset, diagnostic-scale duration or synthetic resource is introduced.
