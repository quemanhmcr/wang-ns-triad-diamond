# Physical `H1` covariant dephasing source calculus

The physical H1 no-escape theorem leaves one named escape: large variation of the interaction-picture H1 forcing.  This note routes that variation into actual Navier--Stokes sources or an already-existing strain-curvature branch.

Let `P_H` denote the fixed orthogonal hook projector on the differentiated-incompressible affine curvature space.  The full normalized curvature satisfies

\[
\dot B=-2A_{aff}B+S,
\qquad
S=L^{-1}\nabla^2(D_t^V V)[L,L].
\]

Therefore

\[
\boxed{
\dot B^H=P_H(-2A_{aff}B+S).
}
\]

No commutation of `P_H` with the affine driver is assumed.  For a physical hook slice

\[
G_c=L B_c^H L^{-1},
\]

`Ldot=A L` gives the exact identity

\[
\boxed{
\dot G_c=A G_c+L\dot B_c^H L^{-1}-G_cA.
}
\]

On `cond(L)<=21/20`, `||A_aff||<=kappa||A||`.  Since `P_H` is orthogonal on the curvature tangent space,

\[
\|\dot B^H\|
\le2\kappa\|A\|\|B\|+\|S\|.
\]

Consequently the physical matrix derivative is bounded by

\[
\|\dot G\|
\le\kappa\|S\|+(2\kappa+2\kappa^2)\|A\|\|B\|.
\]

The good-core triad-normal helical frame has complex derivative at most `(15/2)||A||` on the child and less on the parents.  Its underlying real transverse frame therefore has Frobenius rate at most

\[
15\|A\|/\sqrt2.
\]

Differentiating the trace-free transverse restriction adds twice this frame rate, while the interaction covariant derivative `dot D+[D0,D]` adds at most `2||A||||D||`.  Hence one physical role obeys

\[
\|\nabla_t D\|
\le
\kappa\|S\|+
(4\kappa+2\kappa^2+15\sqrt2\,\kappa)\|A\|\|B\|.
\]

There are three physical roles.  Pulling to their interaction pictures adds at most `e^(2K)` with the existing low-strain action `K<=1/30`.  Arb verifies at `kappa=21/20`

\[
\sqrt3e^{2K}\kappa<2,
\]

\[
\sqrt3e^{2K}
(4\kappa+2\kappa^2+15\sqrt2\kappa)<54.
\]

Thus the final clean H1 variation estimate is

\[
\boxed{
J_1
\le
2\int_0^T\|S\|dt
+54\int_0^T\|A\|\|B\|dt.
}
\]

If the H1 no-escape theorem is in its dephasing branch,

\[
J_1\ge I_1/(11T),
\]

then necessarily

\[
\boxed{
\int\|S\|dt\ge I_1/(44T)
}
\]

or

\[
\boxed{
\int\|A\|\|B\|dt\ge I_1/(1188T).
}
\]

On the full-curvature H1-dominant branch `I1>=I_B/2`, the second alternative implies

\[
\|A\|_{L^\infty_t}T\ge1/2376.
\]

Therefore below this explicit strain-action threshold the curvature source is mandatory.  For the resolved transporter

\[
S=S_P+S_R+S_\nu,
\]

with

\[
S_P\sim L^{-1}(-\nabla^3P)[L,L],\qquad
S_R\sim L^{-1}(-\nabla^2\nabla\cdot R)[L,L],\qquad
S_\nu\sim L^{-1}(\nu\nabla^2\Delta V)[L,L].
\]

By the triangle inequality, one channel has

\[
\boxed{
\int\|S_*\|dt\ge I_1/(132T).
}
\]

The pressure-third channel keeps the already-certified far-field exponent `6-3=3`.  The differentiated-SGS and viscous channels now have an explicit threshold but still require their own critical-mass/dissipation collision estimates.  If the strain-action threshold `1/2376` is crossed, that event is handed to the existing objective-strain/source ledger rather than called a sideband error.
