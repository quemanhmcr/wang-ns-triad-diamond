# High-strain descending epoch: physical PDE probe

Status: **DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__ACTUAL_LOW_PASS_STRAIN_AND_NORMALIZED_DISSIPATION__PHYSICAL_DV_CRITICAL_SHELL_DISINTEGRATION__DESCENDING_HIGH_STRAIN_EPOCH_TELESCOPE**.

The probe evolves the unforced three-dimensional incompressible Fourier--Galerkin
Navier--Stokes system with Leray projection, viscosity, `2/3` dealiasing and RK4.
On the same evolved states it reads the strict low pass `S_(N/4)u`, integrates
the actual strain action `K_N`, normalized resolved dissipation `D_N`, and global
gradient reservoir `G_*`, and disintegrates the positive `D_V` density over the
physical dyadic ancestor shells.  No random operator or recurrence proxy supplies
the reported epoch.

The fixture uses `N=16`, `c=1`,
`nu=0.05`, amplitude `256` and the actual natural
duration `cN^-2`.  Its selected critical shell is `M=N/4` and renewal is
`A=3M/4`.  On this periodic falsification fixture the renewed cutoff `A/4<1`
lies below the first nonzero Fourier mode, so the next high-strain step is
physically absent on the evolved mean-zero Galerkin trajectory.

| resolution | steps | div/||u0|| | energy balance | K_N | D_N | collision margin | reservoir margin | retained D_V fraction | half-law margin | descendant/root grad | certified count upper |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 128 | 6.679e-17 | 4.148e-09 | 2.199e+00 | 3.534e+06 | 9.948e-01 | 3.027e-01 | 1.000000 | 5.000e-01 | 0.000e+00 | 9 |
| 16 | 128 | 6.933e-17 | 5.367e-09 | 2.230e+00 | 3.521e+06 | 9.946e-01 | 3.190e-01 | 1.000000 | 5.000e-01 | 0.000e+00 | 9 |
| 20 | 128 | 6.203e-17 | 6.305e-09 | 2.246e+00 | 3.521e+06 | 9.946e-01 | 3.217e-01 | 1.000000 | 5.000e-01 | 0.000e+00 | 9 |

Root-dissipation resolution spread: `3.777e-03`.

This is a numerical falsification test on an actual finite Galerkin NS system,
not a continuum proof and not evidence that every continuum high-strain epoch
terminates after one step.  The exact continuum conclusion still comes from
`D_*<=D_j<=N_jG_*` and the physical `N_(j+1)/N_j<=3/16` renewal law.
