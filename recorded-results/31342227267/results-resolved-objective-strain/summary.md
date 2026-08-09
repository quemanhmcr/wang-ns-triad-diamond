# Resolved objective-strain source collision

Status: **ARB_CERTIFIED_RESOLVED_OBJECTIVE_STRAIN_COLLISION**.

The objective strain used by the affine/Kelvin packet is the strain of the strict transporter `V=S_(N/4)u`. Its filtered equation gives exactly

`S_circ = -S^2-Omega^2+[S,Omega] - Hess P - sym grad div R + nu Delta S`.

Here `supp Vhat subset B_(N/4)` and `supp Rhat,supp Phat subset B_(N/2)`.  The order-two Hilbert-valued Hausdorff--Young/Bernstein constant is `<1/380`; interpolation gives `||V||_3^2 < mu_V/15`; and `||D^3 V||_inf < N^(7/2)||grad V||_2/1500`. Therefore the four source channels have clean pointwise collisions

- quadratic stretching: `rho_Q <= d_V/(96 pi^2)`, hence `d_V>=96 pi^2 rho_Q`;
- resolved SGS strain source: `rho_R2 <= ||R||_(3/2)/380`, hence `||R||_(3/2)>=380 rho_R2`;
- filtered pressure Hessian: `rho_P <= mu_V/5700+||R||_(3/2)/380`, hence either `mu_V>=2850 rho_P` or `||R||_(3/2)>=190 rho_P`;
- viscosity: `rho_nu <= nu sqrt(d_V)/1500`, hence `d_V >= (1500 rho_nu/nu)^2`.

Thus the old unresolved **near pressure-Hessian coefficient disappears for the actual resolved Kelvin transporter**. Pressure strain-dephasing routes to coherent resolved mass or to the same SGS-increment service currency. The raw full-velocity strain identity remains mathematically valid, but it is not the correct source object for the affine transporter used in the service-or-flat gate.

If `A_obj=T int||S_circ||dt` and `T=cN^-2`, one of these four channels carries scaled source weight at least `A_obj/(4c)`. Quadratic/viscous channels route to critical normalized dissipation; SGS and the stress part of pressure route through Germano/Onsager to coherent service, ancestry or high-frequency enstrophy. The resolved-mass pressure branch enters the existing coherent reservoir/reuse mechanism.

Stress: `50000`
- worst corotational identity residual: `6.874e-16`
- minimum quadratic margin: `1.002e-10`
- minimum SGS margin: `1.002e-10`
- minimum pressure split margin: `5.008e-11`
- minimum viscous margin: `9.268e-13`
