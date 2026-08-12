# Coherent increment service: SGS charge becomes phase-space edges

Status: **EXACT_MOYAL_CELL_ROUTING_GIVEN_STANDARD_LP_BERNSTEIN**.

Let `Q` be the cubic increment charge and `g1=||G||_1`.  Some filter displacement `r` satisfies

`||delta_r u||_3^2 >= (Q/g1)^(2/3)`.

Fix one smooth square-normalized dyadic annular LP analysis--synthesis frame `u_j=phi_j(D)u`, `sum|phi_j|^2=1`, `u=sum phi_j(D)u_j`, with `supp phi_j subset {M_j/2<|xi|<2M_j}` and `|phi_j|<=1`.  Its ordinary `L^3` square-function/Bernstein constants are the supplied finite `C_LP,C_B` for this fixed frame.  In `L^2`, the same frame gives the exact downstream comparison

`D_tail=N int||grad P_>N u||_2^2 dt >= D_high/4`.

With `Y=(Q/g1)^(2/3)/(C_LP C_B)^2`, standard LP/Bernstein gives

`S_low(r)=sum_(j<=0) M_j ||delta_r u_j||_2^2 >= Y-2 d_high`.

For any normalized coherent window and phase-space partition,

`s_(j,C)=M_j int_C |V_g delta_r u_j|^2 dmu`

is a positive **actual increment-service measure** and sums exactly to `S_low`. Translation covariance gives

`V_g(delta_r u)(X,k)=exp(-ik.r)V_g u(X-r,k)-V_g u(X,k)`,

hence

`s_(j,C) <= 2 M_j[E_j(C)+E_j(C-r)]`.

On `d_high<Y/4`, `S_low>=Y/2`.  Once whole-old-pool erosion gives `old_capacity<=Y/8`, either selected old/new interface service is at least `Y/8`, or new/new coherent service is at least `Y/4`.  With the clean quarter-dominance threshold, the latter yields

- a dominant new coherent cluster with critical Moyal mass at least `Y/32`; or
- ancestry Bellman entropy at least `log 2`; or
- same-ancestry hidden pair/cycle mass at least `1/4`.

Thus the old low-band reservoir selected by the Onsager collision is no longer merely an aggregate global mass: the **actual increment itself** generates nearby coherent-cell edges, and after old-pool erosion those edges must cross an interface or create new coherent ancestry.

Stress: `50000`
- worst STFT translation-covariance residual: `6.829e-14`
- worst local increment-capacity ratio: `0.855742519`
- minimum routing margin: `0.000e+00`
- branches: `{'high_frequency_dissipation': 12535, 'new_service_Bellman_entropy': 8327, 'dominant_new_coherent_cluster': 12093, 'selected_interface_Xi': 7056, 'old_pool_not_yet_eroded': 9483, 'new_service_same_ancestry_cycle': 506}`
