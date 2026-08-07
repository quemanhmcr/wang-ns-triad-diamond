# Divergence-free curvature irreducible sideband split

For `B[a,b,c]=B[a,c,b]` with differentiated incompressibility `B[a,a,c]=0`, let `T=Sym B`, `t_c=T[a,a,c]`.  The exact divergence-free tensor carrying the envelope symmetrization is

`B_E = T - 1/2(delta_ab t_c+delta_ac t_b) + delta_bc t_a`.

Then `B_H=B-B_E` is orthogonal to `B_E`, has `Sym B_H=0` and remains divergence free.  Hence it is exactly the five-dimensional quadratic-swirl sector

`B_H[a,b,c]=eps[a,b,d] M[d,c]+eps[a,c,d] M[d,b]`, `M=M^T`, `tr M=0`.

Exact norm identities:
- `||B_E||^2=||T||^2+3||t||^2`;
- sharp symmetric trace bound `||t||^2<=(5/3)||T||^2`, hence `||B_E||^2<=6||T||^2`;
- `||B_H||^2=6||M||^2`;
- normalized hook H1 strain coefficient `C_hook=1/2(B_H+swap_output_input(B_H))` obeys `||C_hook||^2=(1/4)||B_H||^2`.

Therefore the two intrinsic sideband channels obey the aspect-free curvature observability

`||Sym B||^2 + ||C_hook||^2 >= (1/6)||B||^2`.

This is the representation-theoretic statement `15=(7+3)_envelope + 5_swirl`: normalized non-affine curvature cannot disappear simultaneously from the H3 envelope and H1/swirl sectors.  The `C_hook` norm is an intrinsic grain-coordinate sideband norm; comparison with physical Euclidean helicity curvature still uses the existing polarization/ancestry bridge and is not claimed aspect-uniformly here.

Stress: `50000`
- worst envelope-sym residual: `1.438e-15`
- worst hook-sym residual: `1.099e-15`
- worst hook-divergence residual: `1.246e-15`
- worst orthogonality residual: `7.088e-15`
- worst envelope-norm identity residual: `2.132e-14`
- worst hook-M reconstruction residual: `1.470e-15`
- worst hook-norm identity residual: `2.132e-14`
- worst H1-norm identity residual: `5.329e-15`
- minimum combined-observability margin: `4.789e-01`
- maximum sampled trace ratio: `1.435223791` (sharp bound `5/3`)
