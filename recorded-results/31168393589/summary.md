# Affine Gaussian grain dynamics

Status: **CERTIFIED** for the extremal strain-rigidity and frozen-strain cost.

- extremal planar strain coercivity: `> 43/100`
- certified coercivity enclosure: `[0.4359294243 +/- 9.23e-11]`
- frozen dimensionless strain-time: `d T <= 1/25`
- pointwise Hodge coefficient: `> 3/5`
- bracket enclosure: `[0.64139394722215711767964014324840433487233120981 +/- 5.72e-48]`
- time-averaged single-edge deficit: `>= 1/10 (dT)^2`
- adversarial affine/Gaussian checks: `50000`
- worst planar coercivity seen: `0.435929424`
- worst frozen Hodge ratio `H/(d t)^2`: `0.782764352`
- worst Gaussian dual-center residual: `3.010e-14`
- worst log-det residual: `2.653e-14`
- worst Kelvin transversality residual: `5.329e-15`
- worst Kelvin energy-work residual: `1.243e-14`
- worst isotropic-viscosity shape residual: `3.443e-15`

The theorem is a local affine/Kelvin-Gaussian dynamics statement, not a full
Navier--Stokes packet-lifetime theorem.  A common rotation and common planar
scalar strain are gauge directions; the charged quantity is the non-conformal
trace-free strain seen by the extremal triad plane.
