# Localized relative-polarization packet bridge

Status: **CERTIFIED** for the clean low-strain constants.

- Kelvin direction stability: `delta(T)<=exp(4 c sigma0)(h+c kappa M)`
- generator-freezing pair bound: `eps_pol<=sqrt(5)(deltaS_F+16 sigma delta_dir)`
- if `c sigma0<=1/30`, integrated additional polarization forcing:
  `E_pol <= 3 h +(15/2)c kappa M`
- combined with the old localization ledger:
  `E_total <= a/M +(b+15c/2) kappa M +3h`
- optimized spatial width: `M*=sqrt(a/((b+15c/2)kappa))`
- optimized error: `3h+2 sqrt(a(b+15c/2)kappa)`
- random checks: `50000`
- worst Kelvin Lipschitz ratio: `0.939777827`
- worst generator-bound ratio: `0.514395148`
- minimum simplified-bound margin: `4.552e-07`
- worst optimizer residual: `2.220e-16`

Thus localized helical polarization introduces no new spatial error currency on
the low-strain branch: frequency-cell variation contributes a summable `O(h)`
term and spatial frame variation is absorbed into the same `kappa M` curvature
term already balanced against the `a/M` commutator.
