# Smooth relink donor: physical PDE probe

Status: **DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__ACTUAL_RESOLVED_LINEARIZED_OPERATOR__SMOOTH_KPHYS_PAIR_FLUX_AND_DONOR_CLOSURE__NATIVE_WORK_MASTER_REPLAY**.

The trajectory is the unforced 3D incompressible Navier--Stokes Fourier-Galerkin system with Leray projection, viscosity, 2/3 dealiasing and RK4.  At every snapshot the probe applies the actual resolved linearized operator `L_V f=B(V,f)+B(f,V)` and reads its adjoint `K_phys/S` work split through physical pairings.  No proxy evolution or random matrix supplies the reported donor law.  This is numerical falsification evidence, not a continuum proof.

Physical interval: `T=0.01`, `A=4`, `nu=0.05`.

| N | steps | relink snapshots | mixed snapshots | div | global balance | K antisym | K rows | S rows | donor margin | max path | master failures |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 | 48 | 49/49 | 49/49 | 4.367e-17 | 1.008e-11 | 0.000e+00 | 6.628e-15 | 8.865e-15 | -6.196e-15 | 1 | 0 |
| 16 | 48 | 49/49 | 49/49 | 7.223e-17 | 1.173e-11 | 0.000e+00 | 4.427e-15 | 5.154e-15 | -2.315e-15 | 1 | 0 |
| 20 | 48 | 49/49 | 49/49 | 3.966e-17 | 1.173e-11 | 0.000e+00 | 2.626e-15 | 4.026e-15 | -2.196e-15 | 1 | 0 |

Final positive relink-work resolution spread: `1.234e-06`.

The smooth roles are one fixed positive quadratic Fourier partition `eta_a=A_a^2`, the resolved field and every work pairing come from the same evolved PDE state, and the master receives the replayed positive native work while any surviving strain bundle receives only its own positive component.
