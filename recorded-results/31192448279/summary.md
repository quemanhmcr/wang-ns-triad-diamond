# Sideband coherence / first-Duhamel daughter theorem

Status: **CERTIFIED_ARITHMETIC_GIVEN_GAUSSIAN_HYPERCONTRACTIVITY** for the clean Young-footprint arithmetic constants, conditional only on the standard Gaussian hypercontractivity inequality used analytically.

For an interaction-frame sideband forcing `f(t)` in any Banach space,

`int ||f|| <= ||int f|| + (T/2) int ||f_dot||`.

Hence, writing `A=int||f||`, either the first Duhamel daughter impulse has norm at least `A/2`, or the forcing variation is at least `A/T`.

For the affine curvature connection `Bdot+2 A_aff B=S`, pulling back by the connection gives `dot B_tilde=P S`.  The H3 envelope forcing in that pulled-back grain has normalized L2 norm `sqrt(3/8)||Sym B_tilde||` and its dephasing derivative is sourced by `Sym(P S)`.

Critical Young footprint:
- H3 Gaussian sideband: `L^(3/2) relative norm > (1/160) * relative L2 norm`;
- two-component H1 Gaussian sideband: `L^(3/2) relative norm > (1/16) * relative L2/Frobenius norm`.

Thus a coherent H1/H3 first-Duhamel daughter cannot be invisible to the critical Young capacity.  If the first iterate is later cancelled, that cancellation must come from nonlinear sideband/cross interactions rather than from the common affine connection.

Stress checks: `20000`
- worst negative variation margin: `0.000e+00`
- worst exact curvature-pullback residual: `3.252e-15`
- worst H3 pulled-source derivative residual: `3.860e-15`
