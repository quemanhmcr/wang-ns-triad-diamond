# Filtered SGS / viscous source collision

Status: **CERTIFIED**.

For the strict transporter `V=S_(N/4)u`, its SGS stress has Fourier support in `|xi|<=N/2`. Vector-valued Hausdorff--Young/Bernstein and the mild-aspect affine factor give

`N^-4 ||S_R|| <= (3/2000) s ||R||_(3/2)`,  `s=N r_g`,

while the viscous source obeys

`N^-4 ||S_nu|| <= (nu s/5000) sqrt(d_V)`,  `d_V=N^-1||grad V||_2^2`.

Thus, for `s<=s0`, a differentiated-SGS source level `rho_R` forces

`||R||_(3/2) >= 2000 rho_R/(3 s0)`,

and hence, by the exact Germano increment bound, cubic velocity-increment charge at the actual `N/4` filter scale. A viscous source level `rho_nu` forces

`d_V >= (5000 rho_nu/(nu s0))^2`.

If instead `s>s0`, the selected affine grain itself carries scale-critical physical mass `N int_E |u|^2 >= (3/10)s0`; this is a radius-energy/ancestry event, not an aspect defect.

Stress: `50000`
- worst affine-factor / `kappa^2 N r_g` ratio: `0.999999860`
- minimum SGS collision margin: `4.789e-08`
- minimum viscous collision margin: `2.377e-08`
