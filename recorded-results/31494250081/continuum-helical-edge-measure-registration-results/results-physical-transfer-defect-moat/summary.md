# Physical transfer-weighted defect moat

Status: **EXACT_TRANSFER_WEIGHTED_DEFECT_MOAT_GIVEN_PHYSICAL_GOOD_CORE**.

On the `eta_0=10^-4` signed-good core, capacity and actual positive child-transfer laws have density ratio at most `53/50`.  If the block deficit is `epsilon<1/20000`, capacity good-core mass is at least `1/2`, while single-edge stability gives `sum w_e D_e<=2 epsilon`. Hence the normalized **physical** good-core law satisfies

`E_phys[D] <= (106/25) epsilon`.

At one recursive coherent-cell depth choose a scalar defect radius `R`, split `[R/2,R]` into `M` bins, and delete the bin with least physical transfer.  Connect packet/coherent vertices using all edges below the lower moat boundary. Every cross-component edge is then either in the chosen moat or has `D>=R/2`. Therefore

`eta_cross <= 1/M + 2 E_phys[D]/R`.

Choose `M_j=M0(j+2)^2` and `R_j=R0(j+2)^2`. Since `sum_(j>=0)(j+2)^-2<13/20`,

`sum_j eta_cross,j <= (13/20)[1/M0+2 Dbar/R0]`.

On the low-cost extraction branch `C_j<=Ccap`, `xi_j=log(1+eta_j exp(C_j))<=exp(Ccap)eta_j`, so this gives a summable and tunably small physical-transfer `Xi`.

No phase-space packet count or Gaussian synthesis tail is used. If the low-defect hypergraph remains connected instead of splitting, that is not an interface failure: each retained connected component obeys the exact incidence identity `(n-1)+beta=2m` and therefore routes to fresh-rich or cycle-rich ancestry.

Stress: `50000` random weighted triad graphs
- minimum cross-bound margin: `6.917e-04`
- minimum Markov-tail margin: `7.618e-05`
- minimum moat-pigeonhole margin: `3.769e-04`
- minimum square-schedule margin: `7.061e-05`
- incidence checks: `5817`
