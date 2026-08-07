# Affine critical-grain energy ledger

Status: **CERTIFIED_FROM_AFFINE_SHELL_INPUTS**.

Define `r_g=(det Sigma_x)^(1/6)` and the affine scale-critical local mass
`M_aff(E)=r_g^-1 integral_E |u|^2`.  The shell/aspect certificate gives on the
radius-two Gaussian covariance ellipsoid

`M_aff(E2) >= 3/10`.

For fresh grains with `M_aff>=eta` and overlap multiplicity `P`, physical energy
conservation gives the exact budget

`sum r_g <= P E_total / eta`.

The certified shell lower axis also implies, with `s=N r_g` and `A=N l_max`,

`A <= (9/4) s^3`.

Thus an affine grain with natural geometric scale has bounded aspect, while a
very elongated grain necessarily has a large physical geometric radius and is
more expensive in the fresh-energy ledger; no false Young/Bellman anisotropy
penalty is required.

- random checks: `50000`
- worst geometric-radius residual: `7.201e-11`
- local affine-mass coefficient: `0.311109826`
- minimum aspect-relation margin: `6.300e-01`
- minimum fresh-budget margin: `1.627e-05`
