# Onsager cubic-increment to grain/enstrophy collision

Status: **EXACT_SEQUENCE_ALGEBRA_GIVEN_STANDARD_LP_BERNSTEIN**.

For a dyadic Littlewood--Paley decomposition at `N_j=2^j N`, standard `L^3` square-function and Bernstein estimates give

`Q_N <= g1(C_LP C_B)^3[(4/3)(beta R_G)^2 mu_low_max + 2 d_high]^(3/2)`,

where `mu_j=N_j||u_j||_2^2` and `d_high=sum_(j>=1)2^j mu_j = N^-1||grad P_>N u||_2^2` up to the fixed LP partition constants.

Thus if `X=[Q_N/(g1(C_LP C_B)^3)]^(2/3)`, every cubic increment event has the exact alternative

`mu_low_max >= 3X/[8(beta R_G)^2]`

or

`d_high >= X/4`.

A large aggregate band mass is not allowed to hide in many tiny packets.  For packet fractions `w_a`, either one atom has `w_a>=theta`, or `H_atomic>=-log(theta)`.  With ancestry labels the latter routes, by the exact existing collision chain rule, into component Bellman entropy or same-ancestry pair/cycle mass `>=theta^alpha-theta`.

If the enstrophy branch persists on a scaled-time set of measure `m`, normalized dissipation pays at least `m d_high`; failure of persistence is explicitly a temporal-concentration/CKN branch rather than a free escape.

Stress: `50000`
- worst exact dyadic square weight / coarse bound: `0.347252507`
- minimum mass/enstrophy routing margin: `1.135e+01`
- entropy branch counts: `{'ancestry_Bellman_entropy': 30745, 'dominant_packet': 15497, 'same_ancestry_pair_cycle': 3758}`
