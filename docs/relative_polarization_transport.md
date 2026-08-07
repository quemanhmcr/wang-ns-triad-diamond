# Exact relative-polarization transport

At the equal-parent helical extremizer, factor out scalar amplification/damping
from each objective Kelvin spinor.  The remaining determinant-one polarization
variables satisfy

\[
\dot U=-D_1U,\qquad \dot V=-D_2V,\qquad \dot Z=-D_3Z,
\]

where each `D_i` is real symmetric trace free.  Let

\[
J=\begin{pmatrix}0&1\\-1&0\end{pmatrix},
\qquad \lambda=(1,-1)^T.
\]

The parent tensor is the symplectic wedge

\[
W=U^TJV.
\]

Every trace-free `2x2` matrix obeys

\[
D^TJ+JD=0.
\]

Therefore, without any commutativity or frozen-strain assumption,

\[
\boxed{\dot W=U^TJ(D_1-D_2)V.}
\]

A common arbitrary time-ordered `SL(2)` history cancels **pointwise**.  This is
stronger than a Magnus approximation and explains why a norm such as
`||M_i-I||` is not a physical polarization defect.

The signed extremal child-energy tensor also contains the child linear factor
`lambda^T Z`.  Thus

\[
\mathcal P=(U^TJV)(\lambda^TZ)
\]

satisfies the exact full time-ordered identity

\[
\boxed{
\dot{\mathcal P}
=\big[U^TJ(D_1-D_2)V\big](\lambda^TZ)
-(U^TJV)(\lambda^TD_3Z).
}
\]

Only the **relative parent generator** and the child generator enter.  By
Cauchy--Schwarz, `||J||=1`, and `||lambda||=sqrt(2)`, one obtains the clean
capacity-weighted bound

\[
\boxed{
|\dot{\mathcal P}|
\le
2\sqrt{\|D_1-D_2\|_F^2+\|D_3\|_F^2}
\,\|U\|\|V\|\|Z\|.
}
\]

Consequently

\[
|\mathcal P(T)-\mathcal P(0)|
\le
2\int_0^T
\sqrt{\|D_1-D_2\|_F^2+\|D_3\|_F^2}
\,\|U\|\|V\|\|Z\|\,dt.
\]

Scalar transverse strain and viscosity were factored into the amplitude/capacity
ledger before this statement.  The theorem concerns the determinant-one
polarization sector.

A useful countermodel is a large common hyperbolic deformation.  Taking
`D=diag(1,-1)` for time `8` makes the common propagator extremely far from the
identity and very ill-conditioned, while the parent wedge remains exactly
invariant.  Any future packet theorem must therefore estimate the wedge/relative
generator, not Euclidean propagator distance.
\n\n## Forced packet identity\n\nA localized Navier--Stokes packet has residual forcing.  If\n\n\[\n\dot U=-D_1U+F_1,\qquad\n\dot V=-D_2V+F_2,\qquad\n\dot Z=-D_3Z+F_3,\n\]\n\nthe polarization identity acquires exactly the additive term\n\n\[\n\boxed{\n\mathcal R_F\n=(F_1^TJV+U^TJF_2)(\lambda^TZ)\n+(U^TJV)(\lambda^TF_3).\n}\n\]\n\nIt obeys\n\n\[\n\boxed{\n|\mathcal R_F|\n\le\sqrt2\big(\n\|F_1\|\|V\|\|Z\|\n+\|U\|\|F_2\|\|Z\|\n+\|U\|\|V\|\|F_3\|\big).\n}\n\]\n\nThus nonlinear packet forcing is an additive cross-error/fresh-mass ledger term;\nit does not destroy the exact common-parent `SL(2)` cancellation.  The remaining\nPDE task is to estimate the forcing norms from SGS transport, curvature, pressure\nand packet overlap.\n

## Common material phase is not part of `F_i`

For phase variables, a parallel cancellation holds to all spatial orders.  If
all roles use the same resolved velocity `U`, then

\[
(\partial_t+U\cdot\nabla)(\phi_1+\phi_2-\phi_3)
=\rho_1+\rho_2-\rho_3.
\]

At exact carrier resonance, the common velocity-Hessian source cancels from the
chirp/Hessian lock as well.  If the three transport velocities differ, only the
differences `U-U_i` remain.  Consequently the forced polarization ledger should
not count common non-affine phase/chirp evolution inside `F_i`; it belongs to the
osculating packet gauge.  See `docs/material_phase_lock.md`.
\n\n## Field forcing versus spinor forcing\n\nThe symbols `F_i` in this note now mean the **degree-zero Gaussian/helicity\nprojection** of the full PDE field residual.  They must not be identified with\nthe whole residual norm.  `H_1` spatial-polarization curvature and `H_3`\nenvelope curvature have exactly zero base-Gaussian projection and are routed to\nthe sideband/coherence ledger; degree-`<=2` scalar modes are Gaussian tangent.\nAt the transfer/work level the degree-zero parent-profile remainder is already\nthe trilinear cross-error `Xi`; see `docs/hermite_helicity_forcing_ledger.md` and\n`docs/packet_cross_forcing.md`.\n