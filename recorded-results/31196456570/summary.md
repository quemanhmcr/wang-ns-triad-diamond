# Physical H1/swirl bridge on mild-aspect affine grains

Status: **CERTIFIED**.

For the five-dimensional hook/swirl curvature sector, the physical polarization-only observable deliberately excludes the scalar shape term:

`Q_pol=sum_c (||D1(C_c)-D2(C_c)||_F^2+||D3(C_c)||_F^2)`.

At an isotropic grain, Arb certifies

`Q_pol >= (1/10)||B_hook||_F^2`.

For a general grain factor, after polar/global-scalar normalization its spectrum lies in `[cond(L)^(-1/2),cond(L)^(1/2)]`.  The physical slices obey `C_c=Sym(L B_c L^-1)` and

`||C(L)-C(I)|| <= (cond(L)-1)||B||`.

Since `sqrt(Q_pol(S))<=sqrt(5)||S||`, the observable triangle inequality gives

`sqrt(Q_pol(L)) >= [1/sqrt(10)-sqrt(5)(cond(L)-1)] ||B_hook||`.

Hence on the mild-aspect branch `cond(L)<=21/20`, Arb certifies the clean bound

`Q_pol(L) >= (1/25)||B_hook||^2`.

Every real symmetric trace-free 2x2 generator satisfies `D^2=(||D||_F^2/2)I`.  Thus the auxiliary relative-parent/child forcing coordinate has energy `Q_pol/2 >= (1/50)||B_hook||^2`.

For the **three physical role sidebands** needed by the odd-Hermite transfer theorem, parallelogram gives the robust lower bound

`sum_i ||F_i^H1||^2 >= Q_pol/4 >= (1/100)||B_hook||^2`.

The distinction matters because the relative-coordinate evolution is non-unitary and does not itself represent one physical Young role.  No `D_Pi` scalar-shape term is used.  The theorem is intentionally only mild-aspect; larger condition number remains in affine fresh/reuse ancestry.

Stress: `50000`
- worst isotropic `Q_pol/||B||^2`: `0.121618531`
- worst mild-aspect ratio: `0.104335578`
- minimum margin above condition-dependent perturbation bound: `6.255e-02`
- worst spinor-action identity residual: `7.105e-15`
- maximum condition number tested: `1.050000000`
