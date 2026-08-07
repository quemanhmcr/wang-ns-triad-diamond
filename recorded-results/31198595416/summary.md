# Physical H1 covariant source calculus

Status: **CERTIFIED**.

Let `B_H=P_H B` be the five-dimensional hook curvature.  The full affine curvature equation is `Bdot=-2 A_aff B+S`, hence the fixed hook projector gives the exact identity `B_Hdot=P_H(-2 A_aff B+S)`.  With `G_c=L B_H,c L^-1`, `Ldot=A L`, the physical hook matrix obeys `Gdot=A G+L B_Hdot L^-1-G A`.

Using `cond(L)<=21/20`, the good-core triad-normal frame rate, the objective base generator bound and low-strain action `K<=1/30`, the three-role interaction-picture H1 forcing satisfies the clean pointwise/integrated estimate

`J1 <= 2 int||S|| + 54 int ||A|| ||B||`.

Therefore the H1 dephasing branch `J1>=I1/(11T)` forces

`int||S|| >= I1/(44T)`

or

`int||A||||B|| >= I1/(1188T)`.

On the H1-dominant full-curvature branch `I1>=I_B/2`, if `||A||_infty T<1/2376`, the strain-curvature alternative is impossible and the curvature source is mandatory.  Since `S=S_P+S_R+S_nu` with `S_P` from pressure third derivatives, `S_R` from differentiated SGS stress and `S_nu` from viscous fourth derivatives, one source channel has integrated norm at least `I1/(132T)`.

This source theorem does not charge base strain: if the `1/2376` action threshold is crossed, it is handed to the existing objective-strain/source branch.  Pressure-third far-field locality retains the previous `6-3=3` summable exponent.

Stress: `50000`
- worst exact hook-projector linearity residual: `3.472e-15`
- worst exact physical-hook product-rule residual: `2.278e-15`
- minimum clean-density margin: `1.481e-01`
- minimum routing margin: `1.327e-09`
