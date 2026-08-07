# Affine-covariant Gaussian forcing and Hermite projection

- intrinsic curvature tensor: `B=L^-1 H[L.,L.]`
- intrinsic carrier: `q=L^T k`
- exact unprojected quadratic-advection residual is the Wick formula recorded in the JSON
- after allowing Gaussian center/carrier/covariance/chirp to osculate the flow,
  the first transverse forcing is exactly third Hermite chaos:
  `||F_perp||/||psi|| = sqrt(6)/4 ||Sym B||_F <= sqrt(6)/4 ||B||_F`
- random affine checks: `50000`
- worst affine B residual: `4.049e-10`
- worst affine q residual: `1.419e-13`
- worst orthogonal grain-gauge residual: `3.789e-10`
- worst sampled transverse/bound ratio: `0.980183862`
- extreme transformed grain condition number: `1.091e+10`
- extreme affine-invariance residual: `5.608e-16`

A large Euclidean aspect ratio is therefore not itself a forcing cost.  The
physical curvature variable is the Hessian expressed in the grain's own affine
metric.  Quadratic wavefront curvature belongs to the Gaussian tangent manifold;
calling it residual forcing would double-count a packet degree of freedom.
