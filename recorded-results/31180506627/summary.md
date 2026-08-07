# Material triad phase-lock gauge

For three phases transported by the same resolved velocity, the signed triad
phase `Phi=phi1+phi2-phi3` is itself a materially transported scalar.  Exact
carrier resonance and quadratic chirp lock therefore persist under arbitrary
common non-affine advection.  Only differential packet/resolved-flow sources can
dephase the transfer.

- random algebraic checks: `50000`
- worst gradient-lock residual: `7.081e-15`
- worst Hessian/chirp-lock residual: `1.029e-14`
- worst common Hessian-source residual at exact resonance: `0.000e+00`
- minimum differential-source triangle margin: `1.232e-02`

This explains why the quadratic `q.B` term in the affine Gaussian forcing module
belongs to a common phase gauge rather than to the transfer-facing residual.
