# Extremal helicity tensor and symplectic parent gauge

Status: **CERTIFIED** for the transfer-relevant strain observable.

- exact isosceles tensor: `T_(s1,s2,s3)=-i C(r) s3 epsilon_(s1,s2)`
- exact parent symmetry: `(M u)^T J (M v)=u^T J v` for every `M in SL(2)`
- exact relative formula: `(M1 u)^T J (M2 v)=u^T J (M1^-1 M2)v`
- transfer-relevant strain tomography:
  `||D_Pi||^2+||D1-D2||^2+||D_child||^2 >= 1/2 ||S||^2`
- exact positive remainder: `(3/2)d^2+(1+8C-8C^2)b^2+(7-8C)x^2+y^2`
- random checks: `50000`
- worst tensor-factorization residual: `5.551e-16`
- worst common-SL2 invariance residual: `7.553e-14`
- worst relative-matrix residual: `3.478e-14`
- worst observed transfer-relevant ratio: `0.510336893`

This corrects a tempting but false interpretation: absolute helicity conversion
is not by itself a transfer cost.  At equal-parent geometry a common determinant-
one deformation of both parent helicity spinors is an exact symmetry of the
nonlinear parent wedge.  The physical variables are relative parent polarization,
child polarization, and scalar triad shape.
