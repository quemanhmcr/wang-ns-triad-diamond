# Affine low-high Navier--Stokes packet equation

For the exact affine low-frequency jet `V=A(x-X)` in co-moving coordinates, the divergence-free linearized role obeys in Fourier space

`(partial_t-(A^T xi).grad_xi) w_hat = (-A+2 khat khat^T A) w_hat - nu|xi|^2 w_hat`.

Thus the characteristic laws are exactly
`dot k=-A^T k` and `dot a=(-A+2 khat khat^T A)a-nu|k|^2a`.
The pressure correction is exactly the `2 khat khat^T A` term; there is no additional packet pressure force after Leray projection.

In an objective transverse frame with `E^T dot E=-skew(E^T A E)`, the exact two-component coefficient equation is

`dot c=-sym(E^T A E)c-nu|k|^2c`,

the same generator used by the helical polarization ledger.  Hence the PDE low-high linearization, Kelvin carrier dynamics and objective helicity spinor are one identity, not separate modeling assumptions.

Stress checks: `50000`
- worst pressure/Kelvin residual: `8.702e-15`
- worst transversality-rate residual: `3.056e-14`
- worst frame orthogonality-rate residual: `2.062e-15`
- worst frame transverse-constraint residual: `1.790e-15`
- worst objective 2x2-coordinate residual: `2.920e-14`
