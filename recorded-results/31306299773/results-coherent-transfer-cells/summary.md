# Coherent transfer cells: exact nonlinear work and relinking balance

Status: **EXACT_BY_POLARIZED_MOYAL_AND_POSITIVE_CELL_ENERGY**.

Polarized Moyal gives an exact cellwise work decomposition

`W_C = 2 Re int_C V_g f conjugate(V_g F) dmu`,  `sum_C W_C = 2 Re <f,F>`.

If the coherent cells are transported by the common affine phase map `L->ML, X->MX, k->M^-T k`, their intrinsic coordinate `zeta=(L^-1X/2,L^Tk)` is unchanged.  Thus common affine motion creates no coherent-cell interface forcing.

For a piecewise material selected family, switching from one selected cell set to another has jump bounded by the positive Moyal energy in the symmetric difference.  Integrating the exact cell work balance gives

`P_plus <= E_final + P_minus + R_switch`.

Hence positive nonlinear service cannot disappear: at least one of terminal coherent energy, backflow/cancellation, or relinking symmetric-difference energy is at least `P_plus/3`.

Stress: `50000`
- worst polarized-Moyal relative residual: `3.426e-16`
- worst cell-work relative residual: `6.209e-16`
- worst affine phase residual: `9.414e-15`
- worst switch jump / symmetric-difference ratio: `0.993763587`
- minimum one-third routing margin: `7.265e-03`
