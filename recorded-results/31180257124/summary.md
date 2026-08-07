# Quadratic incompressible swirl kernel

For a quadratic normalized velocity `V_a=(1/2)B_abc z_b z_c`, the divergence-free
third-Hermite kernel has the exact representation

`B_abc=eps_abd M_dc+eps_acd M_db`, `M=M^T`, `tr M=0`,

hence `V(z)=z cross (M z)`.  It is tangent to Gaussian level spheres:
`V.z=0`; a carrier sees only the quadratic chirp `q.V`.

- random checks: `50000`
- worst full-symmetry residual: `8.665e-16`
- worst divergence residual: `0.000e+00`
- worst M reconstruction residual: `8.882e-16`
- worst velocity representation residual: `1.068e-14`
- worst radial Gaussian-advection residual: `7.105e-15`
- worst chirp reconstruction residual: `1.421e-14`

This five-dimensional mode is a real dynamical symmetry of the scalar Gaussian
envelope, not missing coercivity.  It must be routed to vector/polarization
transport rather than charged as scalar non-affine forcing.
