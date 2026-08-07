# Objective helical polarization dynamics

- exact objective transverse generator: `-sym(A_perp)` after frame spin cancels `skew(A_perp)`
- circular/helical trace-free strain is off diagonal: first-order deformation is helicity conversion, not Berry phase
- exact commutator: `[D1,D2]=2(delta1 beta2-beta1 delta2) J`
- repository coherence thresholds `eps=1/20`, `dT=1/30` give second-Magnus bound
  `5.694444444e-05`
- random checks: `50000`
- worst objective-generator residual: `4.441e-16`
- worst circular-generator residual: `0.000e+00`
- worst strain-area commutator residual: `5.024e-15`
- worst frozen propagator residual: `4.700e-15`

The second-Magnus number is not advertised as a bound for the full time-ordered
propagator.  It identifies the correct local geometric-phase mechanism: rotation
of polarization appears from noncommuting strain orientations, while first-order
symmetric strain is a helicity mixer.  Failure of strain-orientation coherence is
already routed to the objective-strain/source ledger.
