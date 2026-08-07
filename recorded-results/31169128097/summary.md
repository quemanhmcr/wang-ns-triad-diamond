# Objective-strain source collision

- exact source split: `S_circ = Q(A) - Hess p + nu Delta S`
- clean far pressure-Hessian kernel constant: `< 150`
- far shell exponent after 3D packing: `5-3=2`
- scalar gradient Bernstein constant (lambda=1): `0.100658424209`
- scalar third-derivative Bernstein constant (lambda=1): `0.075026359680`
- packet collision stress cases: `50000`
- worst quadratic threshold residual: `1.110e-16`
- worst viscous threshold residual: `1.110e-16`
- minimum pressure-Hessian collision margin: `1.095e-07`

A 5% strain-coherence failure over `T=cN^-2`, with initial non-conformal strain
`d=sigma N^2`, forces at least one normalized source channel at level
`sigma/(60c)`.  In the stated band-limited packet model the quadratic and
viscous channels force critical mass directly; the pressure-Hessian far field
has a stronger `2^-2n` summable packing gain, while its near field is passed to
the local critical-mass coefficient.
