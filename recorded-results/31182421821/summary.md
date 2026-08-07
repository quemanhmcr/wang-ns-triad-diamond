# Affine ellipsoidal moving-window balance

Transport `Xdot=U(X)` and `Ldot=A(X)L`.  Then the normalized grain coordinate
obeys the exact non-affine remainder identity

`D_t z=L^-1[U(X+Lz)-U(X)-A(X)Lz]`.

For `chi(z/M)` this gives `O(kappa_aff M)` material leakage, while the certified
shell lower axis gives the complementary `O(1/M)` physical-gradient/commutator
scale.  Thus the affine window retains the same square-root balance
`a/M+b kappa_aff M` without an aspect-ratio penalty.

- random checks: `50000`
- worst material identity residual: `1.582e-11`
- worst Taylor leakage/bound ratio: `0.992529785`
- worst shell gradient/bound ratio: `0.999859538`
- worst affine-curvature invariance residual: `1.933e-08`
- worst optimizer residual: `4.441e-16`
- extreme transformed condition number: `1.444e+09`
- extreme leakage/bound ratio: `0.093410985`

This closes the geometry of an ellipsoidal moving moat.  It does not yet prove
the full SGS/filter commutator coefficient or pressure/window overlap for the
actual PDE packetization; those remain continuum terms.
