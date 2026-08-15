# Native curl--current law: one primitive operator grammar for critical Navier--Stokes

Status: **DRAFT STRUCTURAL THEOREM CANDIDATE — exact identities on smooth solutions; no global-regularity claim.**

This note records the smallest operator structure reached so far in the curl-centered programme.
Its purpose is not to accumulate representations.  It is to remove them.

The same Navier--Stokes state generates

- one self-adjoint differential operator, curl;
- one skew Cartan/Lamb action, projected cross product;
- one vorticity curvature/current complex;
- one symmetric strain commutator;
- one critical helicity-flip midpoint operator;
- one quadratic heat generator.

The old language of POWER, strain, triad class, spectral curvature, critical pair creation,
Poisson scale, Reynolds operator, Poynting residual and Gauss coercivity will be retained only
when it names a different **reading** of these primitive objects.  It will not be promoted to a
new physical mechanism.

Throughout, work on smooth mean-zero divergence-free fields on a periodic domain or on
`R^3` with sufficient decay.  Write

\[
C=\operatorname{curl},
\qquad C^*=C,
\qquad C^2=-\Delta,
\]

\[
J_uv:=\mathbb P(u\times v),
\qquad J_u^*=-J_u,
\qquad J_uu=0.
\]

The quadratic stocks used below are

\[
E=\|u\|_2^2,
\qquad
H=\langle u,Cu\rangle,
\qquad
Z=\|Cu\|_2^2,
\]

\[
\Lambda=|C|,
\qquad
\mathsf J=\operatorname{sgn}C,
\qquad
K=\langle u,\Lambda u\rangle,
\qquad
M_3=\|\Lambda^{3/2}u\|_2^2.
\]

The rotational Navier--Stokes equation is

\[
\boxed{
 u_t=J_uCu-\nu C^2u.
}
\]

If

\[
\lambda=H/E,
\qquad
r=(C-\lambda)u,
\]

then `J_u u=0` gives the centered form

\[
\boxed{
 u_t=J_ur-\nu C^2u.
}
\]

Pressure has already been quotiented by the Leray projection.  No shell, packet, owner,
stopping rule or analyst clock is used anywhere in the theorem statements below.

---

## 1. Curl defect is the native nonlinear opening

Put

\[
B=\|r\|_2^2
=Z-H^2/E,
\]

and define the critical Beltrami Gram determinant

\[
\boxed{
\Delta:=EB=EZ-H^2.
}
\]

If

\[
u=\sum_a u_a,
\qquad
Cu_a=a u_a,
\qquad
E_a=\|u_a\|_2^2,
\]

then

\[
\boxed{
\Delta
=\sum_{a<b}(a-b)^2E_aE_b.
}
\]

Thus the defect is exactly the total pairwise signed-curl shear capacity of the state.
The same signed-curl difference is forced by the Euler quadratic law.  If

\[
Cu=a u,
\qquad
Cv=b v,
\]

then

\[
\boxed{
\frac12\bigl(J_uCv+J_vCu\bigr)
=\frac{b-a}{2}\,\mathbb P(u\times v).
}
\]

Equal curl eigenvalues do not transfer energy through the symmetrized Euler interaction.
The nonlinear opening and the interaction coefficient are therefore the same native curl gap.

The raw Lamb field satisfies

\[
L:=u\times\omega,
\qquad
\omega=Cu,
\]

and because `u cross lambda u=0`,

\[
L=u\times r.
\]

Hence

\[
\boxed{
\|L\|_{L^1}^2\le E\|r\|_2^2=\Delta.
}
\]

This is only an `L^1` capacity statement.  It must not be upgraded to an `L^2` current bound;
that would erase the physical concentration problem which remains open.

### 1.1 Tangent gradient of the defect

Assume `B>0` and set

\[
J=\langle r,Cr\rangle,
\]

\[
\boxed{
h=Cr-\frac BEu-\frac JB r.
}
\]

Self-adjointness of curl gives

\[
\langle u,Cr\rangle=B,
\]

hence

\[
\boxed{
h\perp u,
\qquad
h\perp r.
}
\]

The common energy--helicity tangent space is

\[
T_u\mathcal M_{E,H}
=\{u,r\}^{\perp}.
\]

Since

\[
\frac12\nabla B=(C-\lambda)^2u,
\]

one has

\[
\boxed{
\nabla_T(B/2)=h=P_T C^2u.
}
\]

The Euler velocity

\[
F_E=J_ur=\mathbb P(u\times r)
\]

also lies in this tangent space, because Euler preserves `E` and `H`.  Therefore

\[
\boxed{
\frac12 B'_{NL}=\langle h,F_E\rangle.
}
\]

Defect growth is not a separate mechanism: it is the directional derivative of one
state function along the one Euler tangent velocity.

### 1.2 Curl--Euler commutator is physical strain

Because `C` is self-adjoint and `J_u` is skew-adjoint,

\[
[C,J_u]^*=[C,J_u].
\]

On divergence-free test fields direct vector calculus gives

\[
\boxed{
[C,J_u]=2\mathbb P S(u)=:2A_u.
}
\]

Consequently

\[
\boxed{
\langle r,[C,J_u]r\rangle
=2\langle Cr,J_ur\rangle
=2\int r\cdot S(u)r\,dx.
}
\]

The full defect law is

\[
\boxed{
B'
=\langle r,[C,J_u]r\rangle
-2\nu\|Cr\|_2^2.
}
\]

Thus POWER, vortex stretching, strain production and signed-curl spectral curvature are
quadratic readings of one self-adjoint commutator.

### 1.3 First curl--Krylov opening

Normalize

\[
q_0=u/\sqrt E,
\qquad
q_1=r/\sqrt B.
\]

With

\[
\beta_1=\sqrt{B/E},
\qquad
\alpha_1=\langle r,Cr\rangle/B,
\]

and, when `h\ne0`,

\[
\beta_2=\|h\|/\sqrt B,
\qquad
q_2=h/\|h\|,
\]

curl generates

\[
\boxed{Cq_0=\lambda q_0+\beta_1q_1,}
\]

\[
\boxed{Cq_1=\beta_1q_0+\alpha_1q_1+\beta_2q_2.}
\]

Since `F_E\perp q_0,q_1`,

\[
\boxed{
B'_{NL}
=2\beta_2B\sqrt E\,\Gamma,
\qquad
\Gamma=\langle q_2,\mathbb P(q_0\times q_1)\rangle.
}
\]

If `h=0`, then `B'_{NL}=0`; rank one is the Beltrami case, and rank three is the first
curl--Krylov level capable of first-order defect creation.

The intrinsic Gram hierarchy

\[
D_n=\det\bigl(\langle C^iu,C^ju\rangle\bigr)_{0\le i,j\le n}
\]

has the Vandermonde meaning

\[
\boxed{
D_n
=\frac1{(n+1)!}
\int\cdots\int
\prod_{i<j}(a_i-a_j)^2
\prod_i d\rho(a_i).
}
\]

In particular

\[
D_1=\Delta,
\qquad
\boxed{D_2=D_1\|h\|_2^2.}
\]

The hierarchy measures signed-curl complexity of the state; it is not a hierarchy of new
physical mechanisms.


---

## 2. One spectral three-current contains every curl observable

Let `P_a` be the spectral projector of curl.  Write

\[
u=\sum_a u_a,
\qquad
Cu_a=a u_a.
\]

On the mean-zero range of curl define the alternating Cartan form

\[
\Omega(x,y,z)
=\langle C^{-1}x,[y,z]\rangle
=-\int x\cdot(y\times z)\,dx.
\]

For ordered signed-curl nodes `a<b<c`, define

\[
\boxed{
\tau_{abc}
=-2\Omega(u_a,u_b,u_c)
=2\int u_a\cdot(u_b\times u_c)\,dx.
}
\]

The Euler energy rates of that one spectral triple are

\[
\boxed{
(S_a,S_b,S_c)
=\tau_{abc}(c-b,\ a-c,\ b-a).
}
\]

Therefore

\[
S_a+S_b+S_c=0,
\qquad
aS_a+bS_b+cS_c=0.
\]

Energy and helicity nullity are already built into the three-current.

For every real spectral observable `phi`,

\[
M_\phi=\langle u,\phi(C)u\rangle,
\]

and the Euler source is

\[
\boxed{
(M_\phi)'_{NL}
=\sum_{a<b<c}
\tau_{abc}
\det
\begin{pmatrix}
1&1&1\\
a&b&c\\
\phi(a)&\phi(b)&\phi(c)
\end{pmatrix}.
}
\]

Equivalently,

\[
\boxed{
(M_\phi)'_{NL}
=2\langle u,[\phi(C),J_u]Cu\rangle.
}
\]

For `a\ne b`,

\[
\boxed{
P_a[\phi(C),J_u]P_b
=\phi[a,b]\,P_a[C,J_u]P_b.
}
\]

Every non-affine curl observable is therefore functional calculus of the same primitive
commutator `[C,J_u]`.

The affine observables `1` and `a` are the universal kernel:

\[
\phi=1\Rightarrow E'_{NL}=0,
\qquad
\phi=a\Rightarrow H'_{NL}=0.
\]

The first nonlinear cases are

\[
\phi=a^2
\quad\Rightarrow\quad
\text{enstrophy/defect curvature},
\]

and

\[
\phi=|a|
\quad\Rightarrow\quad
\text{critical }H^{1/2}\text{ curvature}.
\]

### 2.1 Fisher velocity of the physical curl measure

Let

\[
E_a=\|u_a\|_2^2,
\qquad
S_a=2\operatorname{Re}\langle u_a,P_aF_E\rangle,
\]

and let

\[
\mathcal K_u=\overline{\{p(C)u\}},
\qquad
F_{spec}=P_{\mathcal K_u}F_E.
\]

Every first-order curl-spectral observable sees only `F_spec`.  Its exact speed is

\[
\boxed{
\mathcal A_{spec}:=\|F_{spec}\|_2^2
=\frac14\sum_{E_a>0}\frac{S_a^2}{E_a}.
}
\]

Define the actual Euler fitness

\[
\boxed{f(a)=\frac{S_a}{2E_a}.}
\]

Then

\[
\boxed{
\partial_t\rho=2(f(a)-\nu a^2)\rho,
}
\]

with

\[
\boxed{
\int f\,d\rho=0,
\qquad
\int af\,d\rho=0,
\qquad
\mathcal A_{spec}=\int f^2d\rho.
}
\]

This is the Fisher--Rao tangent speed of the physical energy measure, not an added probability
or ancestry variable.

### 2.2 Nijenhuis curvature is the same three-frequency law

The divergence-free Lie bracket obeys

\[
\boxed{
CJ_uv=-[u,v].
}
\]

For curl define its Nijenhuis torsion

\[
T_C(u,v)
=[Cu,Cv]-C([Cu,v]+[u,Cv])+C^2[u,v].
\]

If `Cu=a u` and `Cv=b v`, then

\[
\boxed{
T_C(u,v)=(C-a)(C-b)[u,v].
}
\]

Thus at a child curl level `c`,

\[
\boxed{
P_cT_C(u,v)
=(c-a)(c-b)P_c[u,v].
}
\]

The symmetric curvature

\[
\mathscr K_C(x,y)
=T_C(Cx,y)-T_C(x,Cy)
\]

has the full-field pairing

\[
\boxed{
Z'_{NL}=B'_{NL}
=\frac13\langle C^{-1}u,\mathscr K_C(u,u)\rangle
=2\int r\cdot Sr\,dx.
}
\]

So physical stretching, the curl commutator, the three-frequency Vandermonde law and the
failure of curl eigenspaces to close under the Lie bracket are one nonlinear curvature.

### 2.3 Jacobi completion: transfer and reorientation are one bracket

The Cartan form is a closed three-cocycle:

\[
\boxed{d\Omega=0.}
\]

Let `P=P_{\mathcal K_u}`, `Q=I-P`, and `R(x,y)=Q[x,y]`.  Projected Jacobi failure is exactly

\[
\boxed{
J_{\mathcal K}(x,y,z)
=-P\bigl([x,R(y,z)]+[y,R(z,x)]+[z,R(x,y)]\bigr).
}
\]

Thus horizontal spectral transfer and vertical reorientation are complementary pieces of one
Lie bracket.  Overlapping triad coefficients are constrained by Cartan/Jacobi compatibility;
they are not independent knobs.


### 2.4 Fixed-basis form: one Cartan tensor plus diagonal heat

Choose a fixed orthonormal curl eigenbasis `Ce_I=lambda_I e_I` and set

\[
\boxed{f_{IJK}=\Omega(e_I,e_J,e_K),}
\]

which is totally alternating and obeys the weighted Jacobi relations.  Writing
`u=sum_I z_I e_I`, the entire Euler field is

\[
\boxed{
\dot z_I^{NL}
=-\frac12\sum_{J,K}(\lambda_K-\lambda_J)f_{IJK}z_Jz_K.
}
\]

Navier--Stokes adds only

\[
\boxed{-\nu\lambda_I^2z_I.}
\]

Thus every triad phase, spectral current and reorientation is a contraction of one fixed Cartan
tensor, while viscosity is diagonal quadratic heat.  The network is not a changing rulebook.


---

## 3. One closed triad cannot be the ultraviolet escape engine

For a physical closed helical triad, the exact critical root geometry gives

\[
\boxed{
\frac{8|g|^2(a_j-a_k)^2}{r_i(r_j+r_k)}\le1,
\qquad r_i=|a_i|.
}
\]

Summing roots yields

\[
\boxed{
\|\Lambda^{-1/2}F_{spec,\triangle}\|_2^2
\le\frac12E_\triangle K_\triangle.
}
\]

Hence an isolated viscous triad satisfies

\[
\boxed{
K_\triangle(t)
\le K_\triangle(0)
\exp\!\left(\frac{E_\triangle(0)t}{4\nu}\right),
}
\]

with no dependence on absolute Fourier scale.  A single triad cannot create finite-time
critical escape by climbing to the ultraviolet.

This cannot be globalized by summing triads: full-PDE Cartan/Waleffe faces add before the
`H^{-1/2}` square is taken.  Their constructive cross terms are the coherent many-triad
loophole, constrained but not yet bounded by the same Cartan/Jacobi law.


---

## 4. Critical `H^{1/2}` is the positive modulus of helicity

The helicity operator has polar decomposition

\[
\boxed{C=\mathsf J\Lambda,
\qquad\Lambda=|C|,
\qquad\mathsf J=\operatorname{sgn}C.}
\]

Thus

\[
\boxed{K=\langle u,\Lambda u\rangle}
\]

is the canonical positive modulus of helicity.  In the critical coordinate

\[
\boxed{y=\Lambda^{1/2}u,}
\]

one has `K=||y||^2`, `H=<y,Jy>`.  The frozen Euler generator

\[
L_u=\Lambda^{1/2}J_u\mathsf J\Lambda^{1/2}
\]

obeys

\[
\boxed{L_u^*\mathsf J+\mathsf J L_u=0.}
\]

Its ordinary skew/symmetric split is

\[
A_u^{crit}=\tfrac12\Lambda^{1/2}\{J_u,\mathsf J\}\Lambda^{1/2},
\qquad
S_u^{crit}=\tfrac12\Lambda^{1/2}[J_u,\mathsf J]\Lambda^{1/2}.
\]

Here

\[
(A_u^{crit})^*=-A_u^{crit},
\quad[A_u^{crit},\mathsf J]=0,
\]

while

\[
(S_u^{crit})^*=S_u^{crit},
\quad\{S_u^{crit},\mathsf J\}=0.
\]

Therefore

\[
\boxed{y_t=A_u^{crit}y+S_u^{crit}y-\nu\Lambda^2y,}
\]

and only the helicity-flip boost changes the critical norm:

\[
\boxed{\kappa(0)=\langle y,S_u^{crit}y\rangle,
\qquad K'_{NL}=2\kappa(0).}
\]

The equal opposite-helicity production law is simply the block form of this one pseudo-unitary
boost, not a separate mechanism.


---

## 5. One Sobolev commutator law, with a unique critical midpoint

Set

\[
A_u=\mathbb P S(u)=\tfrac12[C,J_u],
\qquad D=\operatorname{ad}_{\log\Lambda}.
\]

For `G_s=\Lambda^sJ_uC\Lambda^{-s}`, split

\[
A_{even}=\tfrac12(A_u+\mathsf J A_u\mathsf J),
\qquad
A_{odd}=\tfrac12(A_u-\mathsf J A_u\mathsf J).
\]

Then

\[
\boxed{
\operatorname{Sym}G_s
=
\frac{\sinh((s-\tfrac12)D)}{\sinh(D/2)}A_{even}
-
\frac{\cosh((s-\tfrac12)D)}{\cosh(D/2)}A_{odd}.
}
\]

At the unique midpoint,

\[
\boxed{
\Sigma_c:=\operatorname{Sym}G_{1/2}
=-\operatorname{sech}(D/2)A_{odd}.
}
\]

Every helicity-preserving symmetric block vanishes exactly at `H^{1/2}`; the remaining
helicity-flip block is suppressed by log-scale separation.  The kernel is

\[
\boxed{
\operatorname{sech}\!\left(\frac12\log\frac rq\right)
=\frac{2\sqrt{rq}}{r+q}.
}
\]

Equivalently,

\[
\boxed{
\Sigma_c
=-2\Lambda^{1/2}\int_0^\infty e^{-t\Lambda}A_{odd}e^{-t\Lambda}dt\,\Lambda^{1/2}.
}
\]

The same map is a probability average under unitary log-scale conjugation, hence a contraction.
This is the intrinsic scale locality of critical growth; no shell cutoff is present.


---

## 6. The same `|curl|` calculus supplies the critical scale and the heat scale

Define

\[
w_t=e^{-t\Lambda}\Lambda u.
\]

Then

\[
\boxed{
d\mu_u(t)=\frac{2\|w_t\|_2^2}{K}\,dt
}
\]

is a probability measure, and if `m=K/E`,

\[
\boxed{
\mathbb E_{\mu_u}t=\frac{E}{2K}=\frac1{2m}.
}
\]

The critical boost rate is the state-generated Poisson-scale average

\[
\boxed{
\frac{\kappa(0)}K
=-\int_0^\infty
\frac{\langle w_t,A_{odd}w_t\rangle}{\|w_t\|_2^2}
\,d\mu_u(t).
}
\]

The Poisson semigroup itself is subordinate to the viscous heat semigroup:

\[
\boxed{
 e^{-t\Lambda}
=
\frac{t}{2\sqrt\pi}
\int_0^\infty
s^{-3/2}e^{-t^2/(4s)}e^{-sC^2}\,ds.
}
\]

Thus the scale filter intrinsic to critical Euler growth is generated by the square root of
the same `C^2` used by viscosity.  No natural-window ontology is needed.

### 6.1 Galilean null

For one Fourier triangle `q+\ell-k=0`, with critical input/output of opposite helicity,

\[
\boxed{
2\sqrt{|k||\ell|}\,|g(q,\ell,-k)|
\le
\frac{3\sqrt6}{16}|q|.
}
\]

A uniform advecting mode therefore creates no critical boost.  An almost-uniform mode loses
its helicity-flip strength linearly in its own wave number.

Dynamically, if `U` is constant,

\[
\boxed{
\Sigma_c(-(U\cdot\nabla)u)
=-i[U\cdot D,\Sigma_c(u)],
\qquad D=-i\nabla.
}
\]

Uniform sweeping only unitarily conjugates the dangerous operator.  It cannot regenerate its
norm or eigenvalues.

---

## 7. Productive Fisher action is the remaining scalar escape cost

Let

\[
p_a=E_a/E,
\qquad
r_a=|a|,
\qquad
m=\mathbb E_p r=K/E.
\]

The Euler fitness satisfies `E_p f=0`.  The critical curvature height is

\[
\boxed{
\kappa(0)=E\operatorname{Cov}_p(r,f).
}
\]

Since

\[
EZ-K^2
=E^2\operatorname{Var}_p(r),
\]

define the productive Fisher action

\[
\boxed{
\mathcal A_{prod}
:=\frac{\kappa(0)^2}{EZ-K^2}
=\frac{\operatorname{Cov}_p(r,f)^2}{\operatorname{Var}_p(r)}.
}
\]

It is exactly the squared regression component of the full spectral fitness which points in the
absolute-curl direction.  Cauchy gives

\[
\boxed{
0\le\mathcal A_{prod}\le\mathcal A_{spec}/E.
}
\]

The physical energy-loss clock is

\[
\boxed{
d\tau_E=-d\log E=2\nu N^2dt,
\qquad
N^2=Z/E.
}
\]

On every positive-energy interval,

\[
\boxed{
\int_{t_0}^{t_1}N^2dt
=\frac1{2\nu}
\log\frac{E(t_0)}{E(t_1)}.
}
\]

The critical mean-curl scale obeys

\[
\boxed{
\frac d{dt}\log\frac KE
\le
\frac{\mathcal A_{prod}}{2\nu N^2}.
}
\]

Therefore

\[
\boxed{
\log\frac{(K/E)(t)}{(K/E)(t_0)}
\le
\frac1{2\nu}
\int_{t_0}^t
\frac{\kappa(0,s)^2}
{N(s)^2[E(s)Z(s)-K(s)^2]}
\,ds.
}
\]

A critical mean-curl escape requires

\[
\boxed{
K/E\to\infty
\Longrightarrow
\int^T
\frac{\kappa(0,t)^2}
{N(t)^2[E(t)Z(t)-K(t)^2]}
\,dt
=\infty.
}
\]

This is a persistence requirement, not an event count.

### 7.1 Exact critical Hilbert square

Because `\Lambda u` lies in the curl-cyclic subspace,

\[
\langle\Lambda^{3/2}u,\Lambda^{-1/2}F_{spec}\rangle
=\kappa(0).
\]

Hence

\[
\boxed{
K'
=-2\nu\left\|
\Lambda^{3/2}u-rac1{2\nu}\Lambda^{-1/2}F_{spec}
\right\|_2^2
+
\frac1{2\nu}\|\Lambda^{-1/2}F_{spec}\|_2^2.
}
\]

Consequently

\[
\boxed{
K(t)\to\infty
\Longrightarrow
\int^T
\frac{\|\Lambda^{-1/2}F_{spec}(t)\|_2^2}{K(t)}dt
=\infty.
}
\]

If

\[
\pi(a)=\frac{|a|E_a}{K},
\]

then

\[
\boxed{
\frac1K\|\Lambda^{-1/2}F_{spec}\|_2^2
=
\mathbb E_\pi\left[
\left(\frac{f(a)}{|a|}\right)^2
\right].
}
\]

The local critical escape currency is therefore Euler fitness per unit curl frequency.

### 7.2 The same square extends to every Sobolev stock

For real `s`, with `K_s=\|\Lambda^su\|_2^2`, the shifted metric satisfies

\[
\boxed{
\nabla_{H^{s-1}}(K_s/2)=\Lambda^2u=C^2u.
}
\]

The corresponding exact square is

\[
\boxed{
K_s'
=-2\nu\left\|\Lambda^{s+1}u-rac1{2\nu}\Lambda^{s-1}F_{spec}\right\|_2^2
+\frac1{2\nu}\|\Lambda^{s-1}F_{spec}\|_2^2.
}
\]

Thus viscosity is the same gradient throughout the Sobolev hierarchy; only the metric reading
changes.  Section 10 rewrites this whole hierarchy as one operator Hilbert scale.


---

## 8. The midpoint operator itself obeys heat

The critical midpoint map `u\mapsto\Sigma_c(u)` is linear and translation-covariant.
Let

\[
D_j=-i\partial_j,
\qquad
\Delta_{op}B
=\sum_{j=1}^3[D_j,[D_j,B]].
\]

Because the kernel translation frequency is exactly the generating velocity wavevector,

\[
\boxed{
\Sigma_c(C^2u)
=\Delta_{op}\Sigma_c(u).
}
\]

Therefore full Navier--Stokes gives

\[
\boxed{
\partial_t\Sigma_c(u)
=\Sigma_c(F_E(u))
-\nu\Delta_{op}\Sigma_c(u).
}
\]

Pure operator heat acts by

\[
(e^{-\tau\Delta_{op}}B)_{k\ell}
=e^{-\tau|k-\ell|^2}B_{k\ell},
\]

or equivalently

\[
\boxed{
 e^{-\tau\Delta_{op}}B
=
\int_{\mathbb R^3}G_\tau(x)
T_xBT_x^{-1}\,dx.
}
\]

It is a Gaussian average of unitary translation conjugations and is norm-contractive.

The nonlinear source has no external object hidden inside it.  If

\[
\nabla u=S+\Omega,
\]

then the Euler strain of the acceleration is

\[
\boxed{
S(F_E)
=-(u\cdot\nabla)S-S^2-\Omega^2-\nabla^2p.
}
\]

Therefore `Sigma_c(F_E)` is the same fixed midpoint transform applied to transport, quadratic
deformation and the incompressibility pressure reaction of the existing velocity gradient.

Pressure and corotation can reorient strain without minting global strain amplitude:

\[
S:[S,\Omega]=0,
\qquad
\int S:\nabla^2p\,dx=0.
\]

Periodic incompressibility also gives the Betchov relation

\[
\boxed{
\int\operatorname{tr}S^3dx
=-\frac34\int\omega\cdot S\omega\,dx.
}
\]

Hence

\[
\boxed{
\frac12\frac d{dt}\|S\|_2^2
=-\frac23\int\operatorname{tr}S^3dx
-\nu\|\nabla S\|_2^2.
}
\]

The same gradient supplies amplitude, orientation and critical regeneration; the representation
changes, the physical source does not.

---

## 9. Critical growth is one Reynolds operator racing identity heat

Move to the critical dissipation coordinate

\[
\boxed{z=\Lambda^{3/2}u.}
\]

Define

\[
\boxed{
\mathscr R_c(u)
:=\nu^{-1}\Lambda^{-1}\Sigma_c(u)\Lambda^{-1}.
}
\]

Then

\[
\boxed{
\mathscr R_c^*=\mathscr R_c,
\qquad
\{\mathscr R_c,\mathsf J\}=0.
}
\]

The full critical balance is exactly

\[
\boxed{
K'
=2\nu\langle z,(\mathscr R_c-I)z\rangle.
}
\]

Thus heat is literally the identity in the dissipation coordinate, and all possible critical
Euler amplification is one self-adjoint Reynolds operator.

The scalar productive Reynolds number is only its Rayleigh quotient:

\[
\boxed{
\frac{\kappa(0)}{\nu M_3}
=
\frac{\langle z,\mathscr R_cz\rangle}{\|z\|_2^2}.
}
\]

Hence

\[
\boxed{
K'>0
\iff
\frac{\langle z,\mathscr R_cz\rangle}{\|z\|_2^2}>1.
}
\]

Because `mathscr R_c` anticommutes with helicity, every nonzero eigenvector is helicity-neutral:

\[
\boxed{
\mathscr R_cv=\lambda v,
\quad\lambda\ne0
\Longrightarrow
\langle v,\mathsf Jv\rangle=0.
}
\]

If

\[
H_3=\langle z,\mathsf Jz\rangle,
\]

then

\[
\boxed{
\left|
\frac{\langle z,\mathscr R_cz\rangle}{\|z\|^2}
\right|
\le
\|\mathscr R_c\|
\sqrt{1-\left(\frac{H_3}{M_3}\right)^2}.
}
\]

The operator heat equation becomes

\[
\boxed{
\partial_t\mathscr R_c
=
\nu^{-1}\Lambda^{-1}\Sigma_c(F_E)\Lambda^{-1}
-\nu\Delta_{op}\mathscr R_c.
}
\]

Pure heat cannot increase the top eigenvalue.  If `v` is a normalized top eigenvector of a
self-adjoint `mathscr R`,

\[
\boxed{
\langle v,\Delta_{op}\mathscr R\,v\rangle
=
2\sum_j
\langle D_jv,
(\lambda_{max}I-\mathscr R)D_jv\rangle
\ge0.
}
\]

Uniform sweeping contributes only a commutator and leaves every eigenvalue unchanged.
Therefore only genuine self-generated deformation can create or replenish an eigenvalue above
one.

---

## 10. On `R^3`, the critical midpoint transform is lossless

Remove viscosity from the Reynolds normalization and define

\[
\boxed{
\mathcal Q_c(u)
:=\Lambda^{-1}\Sigma_c(u)\Lambda^{-1}
=\frac12\Lambda^{-1/2}[J_u,\mathsf J]\Lambda^{-1/2}.
}
\]

Then

\[
\mathscr R_c=\nu^{-1}\mathcal Q_c.
\]

With the unitary Fourier convention

\[
\widehat f(k)
=(2\pi)^{-3/2}\int e^{-ix\cdot k}f(x)dx,
\]

the continuum map is an exact scaled Hilbert--Schmidt isometry:

\[
\boxed{
\|\mathcal Q_c(u)\|_{HS}^2
=\frac1{64}\|u\|_{\dot H^{1/2}}^2
=\frac K{64}.
}
\]

Polarization gives

\[
\boxed{
\langle\mathcal Q_c(u),\mathcal Q_c(v)\rangle_{HS}
=\frac1{64}\langle u,\Lambda v\rangle.
}
\]

The map intertwines physical heat and operator heat:

\[
\boxed{
\mathcal Q_c(C^2u)
=\Delta_{op}\mathcal Q_c(u).
}
\]

Therefore

\[
\boxed{
\sum_j\|[D_j,\mathcal Q_c(u)]\|_{HS}^2
=\frac{M_3}{64},
}
\]

and

\[
\boxed{
\langle\mathcal Q_c(u),\mathcal Q_c(F_E(u))\rangle_{HS}
=\frac{\kappa(0)}{64}.
}
\]

The critical NS balance is literally the Hilbert--Schmidt energy law

\[
\boxed{
\partial_t\mathcal Q_c
=\mathcal Q_c(F_E)-\nu\Delta_{op}\mathcal Q_c,
}
\]

\[
\boxed{
\frac12\frac d{dt}\|\mathcal Q_c\|_{HS}^2
=
\langle\mathcal Q_c,\mathcal Q_c(F_E)\rangle_{HS}
-
\nu\sum_j\|[D_j,\mathcal Q_c]\|_{HS}^2.
}
\]

Multiplying by `64` gives

\[
\frac12K'=\kappa(0)-\nu M_3.
\]

### 10.1 Spectral capacity

Since `Q_c` is self-adjoint and helicity-odd, its nonzero eigenvalues occur in pairs.
If the positive eigenvalues of `mathscr R_c` are `lambda_j>0`, then

\[
\boxed{
K=128\nu^2\sum_j\lambda_j^2.
}
\]

Hence

\[
\boxed{
\|\mathscr R_c\|_{op}
\le
\frac{\sqrt K}{8\sqrt2\,\nu},
}
\]

and

\[
\boxed{
N_{>1}^+
\le
\frac{K}{128\nu^2}.
}
\]

In particular,

\[
\boxed{
K<128\nu^2
\Longrightarrow
\|\mathscr R_c\|<1
\Longrightarrow
K'\le0
}
\]

on every smooth continuum interval where the identities apply.

This is a small-critical-data operator regime, not a large-data regularity theorem.

### 10.2 Lossless state representation, but no eigenvalue-only closure

For `T:u\mapsto\mathcal Q_c(u)`, polarization gives

\[
\boxed{T^*T=\frac1{64}\Lambda.}
\]

Hence on mean-zero critical states

\[
\boxed{u=64\Lambda^{-1}T^*\mathcal Q_c(u).}
\]

The transform is injective.  But every image operator is self-adjoint and helicity-odd, so

\[
\boxed{\operatorname{Tr}Q^{2m+1}=0}
\]

whenever defined.  Since `\kappa(0)` is generally nonzero, critical regeneration cannot be an
eigenvalue-only cubic Riccati law.

### 10.3 One operator Hilbert scale contains all Sobolev stocks

Functional calculus gives

\[
\boxed{T(\Lambda^\alpha u)=\Delta_{op}^{\alpha/2}T(u).}
\]

Therefore

\[
\boxed{
\|u\|_{\dot H^s}^2
=64\|\Delta_{op}^{s/2-1/4}\mathcal Q_c(u)\|_{HS}^2.
}
\]

The main stocks are simply

\[
\boxed{
(E,K,Z,M_3)
=64\bigl(
\|\Delta_{op}^{-1/4}Q_c\|_{HS}^2,
\|Q_c\|_{HS}^2,
\|\Delta_{op}^{1/4}Q_c\|_{HS}^2,
\|\Delta_{op}^{1/2}Q_c\|_{HS}^2
\bigr).
}
\]

Hence `H^{1/2}` is the unique physical Sobolev level which becomes plain operator `L^2`
energy.  Moment inequalities such as

\[
\boxed{K^3\le E^2M_3,}
\qquad
\boxed{Z^2\le KM_3}
\]

are ordinary interpolation in this one operator heat scale.


---

## 11. The smaller parent: a curved graded current already present in NS

The preceding curl/operator laws have a still smaller Cartan--Hodge parent.
Let

\[
\alpha=u^\flat,
\qquad
\beta=d\alpha,
\qquad
d\beta=0.
\]

For a real parameter `theta` define

\[
\boxed{
q_\theta=\nu\delta+\theta\,\iota_u,
\qquad
q_\theta^*=\nu d+\theta\,\alpha\wedge.
}
\]

These are not imported operators: `nu delta beta` and `i_u beta` are exactly the viscous and
Euler pieces of the physical electromotive current.

The exterior Leibniz identity is

\[
\{d,\alpha\wedge\}=d\alpha\wedge=\beta\wedge.
\]

Taking adjoints gives

\[
\boxed{
\{\delta,\iota_u\}=(\beta\wedge)^*.
}
\]

Since `delta^2=i_u^2=0`,

\[
\boxed{
q_\theta^2
=\theta\nu(\beta\wedge)^*.
}
\]

Thus vorticity is literally the failure of the NS current differential to be nilpotent.

### 11.1 Actual NS current

The physical member is

\[
\boxed{
q_1=\nu\delta+\iota_u,
\qquad
e=q_1\beta.
}
\]

The rotational momentum and vorticity equations are

\[
\boxed{
\alpha_t+dB=-q_1\beta,
}
\]

\[
\boxed{
\beta_t=-d q_1\beta.
}
\]

On the state itself,

\[
\boxed{
q_1^2\beta=\nu|\beta|^2.
}
\]

So the state current chain is

\[
\boxed{
\beta
\xrightarrow{q_1}
e
\xrightarrow{q_1}
\nu|\beta|^2.
}
\]

The same operator generates the vorticity PDE:

\[
\boxed{
\{d,q_1\}
=\mathcal L_u+\nu L.
}
\]

Associativity gives

\[
[\{d,q_1\},q_1]
=dq_1^2-q_1^2d.
\]

On `d beta=0`,

\[
\boxed{
[\mathcal L_u+\nu L,q_1]\beta
=\nu d|\beta|^2.
}
\]

Concentration therefore forces a noncommutation between the current and the generator made by
that same current.  This is a literal algebraic self-frustration law.

### 11.2 Actual acceleration already obeys a Gauss law

In vector notation

\[
e^\sharp
=-u\times\omega+\nu\operatorname{curl}\omega,
\]

and

\[
u_t+\nabla B=-e^\sharp.
\]

The self-return identity gives

\[
\boxed{
-\nabla\cdot e^\sharp
+\frac{u}{\nu}\cdot e^\sharp
=|\omega|^2.
}
\]

Equivalently,

\[
\boxed{
\left(\nabla\cdot-\frac{u}{\nu}\cdot\right)
(u_t+\nabla B)
=|\omega|^2.
}
\]

The physical NS acceleration itself is a positively sourced curvature current.

---

## 12. The canonical midpoint current contains Schrödinger, strain, Poynting and Gauss

The exact heat/Euler reflection square selects

\[
\boxed{
q=q_{1/2}
=\nu\delta+\frac12\iota_u,
\qquad
q^*=\nu d+\frac12u^\flat\wedge.
}
\]

Its curvature square is

\[
\boxed{
q^2
=\frac\nu2(\beta\wedge)^*.
}
\]

Put

\[
\boxed{
L_u^{Sch}
=-\nu^2\Delta+\frac{|u|^2}{4}.
}
\]

On scalars,

\[
\boxed{
qq^*\big|_{\Omega^0}=L_u^{Sch}.
}
\]

The drift cross term cancels exactly by incompressibility.

On one-forms,

\[
\boxed{
H^{(1)}:=q^*q+qq^*
=L_u^{Sch}I+\nu S.
}
\]

After Hodge identification of two-forms with vectors,

\[
\boxed{
H^{(2)}:=q^*q+qq^*
=L_u^{Sch}I-\nu S.
}
\]

Therefore

\[
\boxed{
2\nu S
=H^{(1)}-\star^{-1}H^{(2)}\star,
}
\]

and

\[
\boxed{
2L_u^{Sch}I
=H^{(1)}+\star^{-1}H^{(2)}\star.
}
\]

Both Hodge squares are positive, so

\[
\boxed{
-L_u^{Sch}I
\le\nu S
\le L_u^{Sch}I
}
\]

in quadratic-form order.

The symmetric velocity gradient is therefore the degree imbalance of one positive graded
current square.

### 12.1 Poynting and Gauss are the degree-two reading

On vorticity,

\[
\boxed{
(q\beta)^\sharp
=\nu\operatorname{curl}\omega-rac12u\times\omega
=-\frac12G,
}
\]

where

\[
\boxed{G=u\times\omega-2\nu\operatorname{curl}\omega.}
\]

Together with

\[
\boxed{q^2\beta=\frac\nu2|\beta|^2,}
\]

this gives

\[
\boxed{
\left(\nabla\cdot-\frac{u}{2\nu}\cdot\right)G=|\omega|^2.
}
\]

Thus the completed Poynting square and the Gauss source are degree-two readings of the same
midpoint current.

### 12.2 Mixed Gauss law acts on every strain direction

For any smooth vector field `b`,

\[
\boxed{\mathfrak G_\nu[u;b]=u\times b-2\nu\operatorname{curl}b}
\]

obeys

\[
\boxed{
\left(\nabla\cdot-\frac{u}{2\nu}\cdot\right)\mathfrak G_\nu[u;b]
=\omega\cdot b.
}
\]

For divergence-free `b`,

\[
\boxed{
\nu\langle b,Sb\rangle
=\nu^2\|Cb\|_2^2+rac14\|u\times b\|_2^2-rac14\|\mathfrak G_\nu[u;b]\|_2^2.
}
\]

Hence, in quadratic-form notation,

\[
\boxed{
4\nu A_u
=P(4\nu^2C^2+U_u^*U_u-\mathfrak G_u^*\mathfrak G_u)P.
}
\]

On the helicity-odd block the even heat term disappears.  The result is a difference of Gram
currents, not a negative square; odd projection does not preserve positivity.

### 12.3 Canonical curvature floor

Let

\[
D_u
=\nabla\cdot-\frac{u}{2\nu}\cdot,
\qquad
M_\omega b=\omega\cdot b.
\]

The mixed Gauss law is

\[
D_u\mathfrak G_u=M_\omega.
\]

Its scalar normal operator is

\[
\boxed{
D_uD_u^*
=-\Delta+\frac{|u|^2}{4\nu^2}.
}
\]

Hence, wherever the inverse is defined in the relevant energy space,

\[
\boxed{
\mathfrak G_u^*\mathfrak G_u
\ge
M_\omega^*
\left(
-\Delta+\frac{|u|^2}{4\nu^2}
\right)^{-1}
M_\omega.
}
\]

A dangerous strain direction with nontrivial `omega dot b` overlap therefore carries a
canonical curvature floor.

---

## 13. The midpoint current also gives the pressure-free critical tax

Let

\[
F_E=P(u\times\omega),
\qquad
G=u\times\omega-2\nu C\omega.
\]

The endpoint currents

\[
j_0=\nu C^2u,
\qquad
j_1=\nu C^2u-F_E=-u_t
\]

have midpoint

\[
\boxed{j_{1/2}=\nu C^2u-\tfrac12F_E=-\tfrac12PG.}
\]

For every `K_s=||Lambda^s u||^2`,

\[
\boxed{
K_s'
=\frac1{2\nu}
\left(
\|\Lambda^{s-1}F_E\|_2^2
-
\|\Lambda^{s-1}PG\|_2^2
\right).
}
\]

At criticality,

\[
\boxed{
K'
=\frac1{2\nu}
\left(
\|\Lambda^{-1/2}F_E\|_2^2
-
\|\Lambda^{-1/2}PG\|_2^2
\right).
}
\]

Euler energy nullity forces

\[
\boxed{
\|\Lambda^{-1/2}PG\|_2^2
\ge\frac{4\nu^2Z^2}{K},
}
\]

so

\[
\boxed{
K'
\le\frac1{2\nu}\|\Lambda^{-1/2}F_E\|_2^2
-2\nu\frac{Z^2}{K}.
}
\]

Helicity nullity supplies a second compulsory projection of the same residual onto the Gram
plane generated by the two invariant directions.  This is the two-null tax; it introduces no
new stock.

There is also a sharp no-go.  Since

\[
q_\theta^*\phi=\nu d\phi+\theta u^\flat\phi,
\]

requiring this scalar Gauss test to be divergence-free gives

\[
\nu\Delta\phi+\theta u\cdot\nabla\phi=0.
\]

On periodic/decaying domains incompressibility implies

\[
\boxed{\phi=\text{constant}.}
\]

Thus the constant test is the only scalar Gauss test living entirely in the Leray-horizontal
sector.  Static Schrödinger coercivity alone cannot close the projected critical gap; the missing
argument must control persistence dynamically.


---

## 14. Material coordinates show where the same strain writes heat memory

The critical operator does not live apart from the material vorticity law.
Let `Phi_t` be the incompressible flow map and define

\[
\widetilde\varpi=\Phi_t^*\varpi,
\qquad
 g_t=\Phi_t^*g_0.
\]

For the vorticity two-form, Navier--Stokes becomes

\[
\boxed{
\partial_t\widetilde\varpi
=\nu\Delta_{g_t}\widetilde\varpi.
}
\]

There is no explicit Euler stretching source left.  Euler moves the material metric:

\[
\boxed{
\partial_t g_t
=2\Phi_t^*S,
\qquad
\det g_t=1.
}
\]

Thus the same physical strain which forms the midpoint critical Reynolds operator is exactly the
metric velocity seen by material Hodge heat.

The feedback is one closed chain:

\[
\boxed{
\widetilde\varpi
\longrightarrow u
\longrightarrow S
\longrightarrow g_t
\longrightarrow\Delta_{g_t}
\longrightarrow\widetilde\varpi.
}
\]

The separate material-current note proves the transverse determinant and Minkowski memory laws.
The present curl note needs only their structural meaning: sustained critical reorientation is
not free of history, because the same strain continuously rewrites the metric in which heat acts.

---

## 15. Falsification guards: what the primitive law does **not** say

The reduction is useful only if its negative statements are kept explicit.

1. **No global-regularity claim.**  Concentrated many-mode deformation can still create a large,
   time-dependent Reynolds operator and reorient the state toward expanding directions.
2. **No spectral-only monotone beyond energy/helicity.**  Physical phase reversal flips every
   non-affine triad current while leaving the signed-curl energy measure fixed.
3. **No universal Cauchy, Schrödinger or pressure gap.**  The exact state-dependent bounds can
   approach saturation, and Leray projection has no fixed loss fraction.
4. **No low Krylov-rank closure.**  High curl--Krylov fitness components can reconfigure the state
   even when they do not instantaneously change enstrophy.
5. **No automatic viscous dephasing.**  Heat damps amplitudes and translation-frequency carriers;
   it does not rotate normalized curl-eigenspace directions by fiat.
6. **No unrestricted Nambu--Poisson theorem.**  The exact structure retained is Cartan/Jacobi.
7. **No eigenvalue-only nonlinear closure.**  The lossless midpoint operator retains orientation
   and translation structure which its eigenvalues forget.
8. **No negative-square theorem after helicity projection.**  The odd Gram difference has no fixed
   sign, so the graded-current identity does not imply `||mathscr R_c||<=1` at large data.
9. **No purely local finite-dimensional closure.**  Any regularity theorem must use the full
   spatial current, its transport compatibility and its heat coupling.


---

## 16. Minimal ontology and the remaining theorem

After all collapses, the primitive objects are

\[
\boxed{
C=\operatorname{curl},
\qquad
J_uv=\mathbb P(u\times v),
\qquad
q_1=\nu\delta+\iota_u.
}
\]

Their basic laws are

\[
\boxed{u_t=J_uCu-\nu C^2u,}
\qquad
\boxed{[C,J_u]=2\mathbb P S,}
\]

\[
\boxed{
\beta_t=-dq_1\beta,
\qquad
q_1^2=\nu(\beta\wedge)^*.
}
\]

The unique critical symmetric Euler operator is

\[
\boxed{
\Sigma_c
=-\operatorname{sech}\!\left(
\frac12\operatorname{ad}_{\log|C|}
\right)(\mathbb P S)_{odd},
}
\]

and

\[
\boxed{\mathscr R_c=\nu^{-1}\Lambda^{-1}\Sigma_c\Lambda^{-1}.}
\]

Critical growth is exactly

\[
\boxed{
K'
=2\nu\langle\Lambda^{3/2}u,
(\mathscr R_c-I)\Lambda^{3/2}u\rangle.
}
\]

On `R^3`, the same state is represented losslessly by

\[
\boxed{
\mathcal Q_c=\Lambda^{-1}\Sigma_c\Lambda^{-1},
\qquad
\|\mathcal Q_c\|_{HS}^2=K/64.
}
\]

The midpoint current simultaneously obeys

\[
\boxed{q_{1/2}^2=\frac\nu2(\beta\wedge)^*,}
\]

and

\[
\boxed{
\left(\nabla\cdot-\frac{u}{2\nu}\cdot\right)
\mathfrak G_\nu[u;b]=\omega\cdot b.
}
\]

Thus the apparently separate theories collapse to one chain:

\[
\boxed{
\text{state}
\to\text{current}
\to\text{curvature}
\to\text{strain}
\to\text{critical Reynolds boost}
\to\text{operator/material heat}.
}
\]

The unresolved theorem is now narrow:

> **Can the state sustain a helicity-odd Reynolds expansion above identity heat for the infinite
> productive action required by critical escape, while the same curvature/current law forces
> reflected-current cost, Cartan/Jacobi compatibility forces reorientation structure, and the
> same strain writes material heat memory?**

Equivalently one needs a dynamical implication of the form

\[
\boxed{
\text{persistent critical backward-heat alignment}
\Longrightarrow
\text{integrable current-turning / curvature-memory cost}.
}
\]

If it forces

\[
\int^T
\frac{\kappa(0,t)^2}
{N(t)^2[E(t)Z(t)-K(t)^2]}
\,dt<\infty,
\]

then the necessary escape condition of Section 7 is contradicted.  No such persistence theorem
is proved here.

The methodological rule is therefore strict:

\[
\boxed{
\text{do not add a mechanism unless Navier--Stokes itself generates the object.}
}
\]

Curl, Cartan contraction, codifferential, strain, helicity modulus, Poisson calculus, operator
heat and material Hodge heat survive this test.  A future proof should go deeper into this closed
current/curvature law rather than rebuild shells, owners or clocks around it.
