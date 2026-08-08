# Physical multi-currency master telescope

Status: **EXACT_MULTICURRENCY_EPISODE_TELESCOPE**.

Separate non-flat resets into two kinds. A multiplicative transfer-cost block pays `C_j>=c0`. An additive physical-resource reset is assigned **one** primary currency `r`, consumes at least `b_r`, and that currency has total global budget `B_r`. Hence

`N_A <= sum_r B_r/b_r`.

If flat blocks erode barycentric potential by `kappa0` up to total perturbation `Z`, while every transfer/additive reset may restart the potential below `Pmax`, then

`N_F kappa0 <= (N_T+N_A+1)Pmax+Z`.

Solving with `L=N_F+N_T+N_A` gives

`N_T >= [kappa0 L-Pmax-Z]/[kappa0+Pmax]-N_A`.

Therefore, with one global cross/interface penalty `Xi`,

`-log prod_(j<L) R_j >= c0 N_T-Xi`

and the asymptotic depth rate remains

`c_eff=c0 kappa0/(kappa0+Pmax)>0`.

Finite fresh-energy, dissipation or initial-boundary resets change only the finite offset through `sum B_r/b_r`; they cannot destroy the linear-in-depth transfer-cost rate. Entropy/Hodge/resistance/Renyi events that already pay multiplicative Bellman cost remain in `N_T`, not in the additive resource count.

Stress: `50000` synthetic multi-ledger episodes
- minimum transfer-count margin: `0.000e+00`
- minimum log-efficiency margin: `0.000e+00`
- minimum resource-count margin: `1.414e-04`
- minimum sampled asymptotic rate: `2.862e-08`
