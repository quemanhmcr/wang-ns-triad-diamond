# SGS / graded spectral flux regression

The equality is analytic; this pseudospectral run validates conventions.

- random divergence-free fields: `200` on `12^3`
- worst SGS-vs-resolved relative error: `4.903e-16`
- worst Leray-projection work relative error: `5.419e-16`
- worst relative divergence diagnostic: `7.735e-16`

The tested identity is `<Pi_SGS>=<ubar . div overline(u tensor u)>`; applying
the Leray projector to the nonlinear term leaves the global work unchanged.
