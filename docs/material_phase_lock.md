# Material triad phase lock: common non-affine advection is a gauge

Let three packet phases satisfy

\[
(\partial_t+U\cdot\nabla)\phi_i=\rho_i.
\]

The transfer phase is the signed combination

\[
\Phi=\phi_1+\phi_2-\phi_3.
\]

Linearity of the material derivative gives the exact identity

\[
\boxed{
(\partial_t+U\cdot\nabla)\Phi=\rho_1+\rho_2-\rho_3.
}
\]

Thus common advection is a phase gauge to **all spatial orders**.  If the three
roles are initially pointwise phase locked and have no differential phase
source, the lock remains exact even when `U` is nonlinear in space.

For the local carrier `k_i=grad phi_i`,

\[
D_tk_i=-(\nabla U)^Tk_i+r_i.
\]

Hence

\[
\boxed{
D_t(k_1+k_2-k_3)
=-(\nabla U)^T(k_1+k_2-k_3)+(r_1+r_2-r_3).
}
\]

At Hessian/chirp level, with `K_i=Hess phi_i` and
`H_{a jk}=partial_j partial_k U_a`,

\[
D_tK_i=-A^TK_i-K_iA-H[k_i]+G_i,
\]

where `H[k]_{jk}=k_aH_{a jk}`.  Therefore

\[
\boxed{
D_tK_{lock}
=-A^TK_{lock}-K_{lock}A-H[k_{lock}]+G_{lock}.
}
\]

Exact carrier resonance `k_lock=0` removes the common velocity-Hessian source.
This is precisely the cancellation behind the quadratic `q.B` chirp in
`docs/affine_gaussian_forcing.md`: each role acquires wavefront curvature, but
the common curvature is compatible with triad phase lock.

The PDE-facing phase residual must therefore be built from **differential**
resolved velocities, SGS forcing, packet partitions and other role-dependent
sources.  Common non-affine advection should not be charged twice.
