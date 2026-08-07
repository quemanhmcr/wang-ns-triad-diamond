# Divergence-free localized affine packet tests: pressure stays out of the microscopic forcing

A strong spatial product `chi w` is generally not divergence free even when `w`
is.  Therefore the microscopic packet equation should not be derived by
multiplying a helical mode by a window and then silently discarding pressure.
Use a divergence-free **test packet** instead.

Let `M_{i,N}(D)` be a smooth matrix Fourier multiplier supported on the certified
role shell, containing the Leray/helical projection and the frozen frequency
cell.  The shell is bounded away from zero, so the complete symbol is smooth and
its physical kernel is rapidly decaying.  Define

\[
\boxed{
\Psi_i=M_{i,N}(D)(\chi_{L,M}\phi_i).
}
\]

Then `div Psi_i=0` exactly.  For the full Navier--Stokes equation,

\[
\partial_tu+\nabla\cdot(u\otimes u)+\nabla p=\nu\Delta u,
\]

the weak coefficient equation is therefore

\[
\boxed{
\frac d{dt}\langle u,\Psi_i\rangle
=\langle u,\partial_t\Psi_i\rangle
+\langle u\otimes u,\nabla\Psi_i\rangle
-\nu\langle\nabla u,\nabla\Psi_i\rangle.
}
\]

The pressure term is exactly zero.  This is the correct microscopic equation for
a spatially localized role.

## Localizing the Leray/helical multiplier costs only the existing moat error

Let the scale-covariant shell multiplier have kernel

\[
K_N(x)=N^3K(Nx),
\qquad m_1(K)=\int|x||K(x)|dx<\infty.
\]

For every Lipschitz window,

\[
[\chi,M_N]f(x)
=\int K_N(y)[\chi(x)-\chi(x-y)]f(x-y)dy,
\]

and Young's inequality gives

\[
\boxed{
\|[\chi,M_N]f\|_2
\le\frac{m_1(K)}N\|\nabla\chi\|_\infty\|f\|_2.
}
\]

For the certified affine ellipsoid window,

\[
N^{-1}\|\nabla\chi\|_\infty\le\frac{3C_\chi}{2M},
\]

so

\[
\boxed{
\|[\chi,M_N]f\|_2
\le\frac{3m_1(K)C_\chi}{2M}\|f\|_2.
}
\]

This applies not only to the coarse-graining filter but to any **shell-localized
Leray/helical/cell multiplier**.  Leray by itself has a singular long-range
kernel, but after multiplication by the compact smooth shell symbol away from
`xi=0`, the combined kernel has finite moments.

Thus spatially localizing the microscopic divergence-free role introduces no
new pressure source and no aspect-ratio defect.  It adds to the same `1/M`
commutator ledger already balanced against affine curvature.

## Two different appearances of pressure

This produces a clean two-level rule:

1. **microscopic full-velocity role coefficient:** pressure is annihilated by the
   divergence-free test packet;
2. **macroscopic resolved local energy:** pressure is spatial transport and is
   handled by combined work `G=Pi+div(PU)` plus the pressure-cancellation annular
   charge.

These are compatible descriptions of the same Navier--Stokes pressure and must
not be added as two independent packet forces.
