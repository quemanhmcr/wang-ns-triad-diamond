# Native curl--current law: one primitive operator grammar for critical Navier--Stokes

Status: **DRAFT STRUCTURAL THEOREM CANDIDATE — exact identities on smooth solutions; no global-regularity claim.**

This note records the smallest operator structure reached so far in the curl-centered programme.
Its purpose is not to accumulate representations.  It is to remove them.

The same Navier--Stokes state generates one actual graded current operator

\[
Q=\nu\delta+\iota_u,
\]

whose curvature, Hodge square and state velocity recover the vorticity/strain/Lamb geometry.
Curl, the Cartan three-current, the critical midpoint/Reynolds operator and quadratic heat are
retained below only as exact readings of this smaller endogenous law.

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

Define first the **lossless** critical Euler operator and then its dimensionless comparison
with molecular heat,

\[
\boxed{
\mathscr E_c(u):=\Lambda^{-1}\Sigma_c(u)\Lambda^{-1},
\qquad
\mathscr R_c(u):=\nu^{-1}\mathscr E_c(u).
}
\]

Both are self-adjoint and helicity-odd.  For `M_3>0` the state itself selects the signed coefficient

\[
\boxed{
\nu_E(u):=\frac{\kappa(0)}{M_3}
=\frac{\langle z,\mathscr E_c(u)z\rangle}{\|z\|_2^2}.
}
\]

This has the dimensions of viscosity but is not a constitutive parameter: it is the actual Euler
Rayleigh coefficient generated by the current state.  Its dimensionless version is exactly the old
productive Reynolds quotient,

\[
\boxed{
\frac{\nu_E}{\nu}
=\frac{\kappa(0)}{\nu M_3}
=\frac{\langle z,\mathscr R_cz\rangle}{\|z\|_2^2}.
}
\]

The full critical balance therefore collapses to

\[
\boxed{K'=2M_3(\nu_E-\nu).}
\]

Thus Euler supplies one signed two-way state coefficient `nu_E`, molecular heat supplies the fixed
one-way coefficient `nu`, and

\[
\boxed{K'>0\iff \nu_E>\nu.}
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

## 10. The helicity-odd midpoint transform is a compressed critical representation

The Reynolds/Krein operator remains exact, but it is no longer the lowest ontology.  Define

\[
\boxed{
\mathcal Q_c(u)
:=\Lambda^{-1}\Sigma_c(u)\Lambda^{-1}
=\frac12\Lambda^{-1/2}[J_u,\mathsf J]\Lambda^{-1/2},
\qquad
\mathscr R_c=\nu^{-1}\mathcal Q_c.
}
\]

On continuum `R^3`, with the unitary Fourier convention,

\[
\boxed{
\|\mathcal Q_c(u)\|_{HS}^2=\frac K{64},
\qquad
\langle\mathcal Q_c(u),\mathcal Q_c(v)\rangle_{HS}
=\frac1{64}\langle u,\Lambda v\rangle.
}
\]

The transform intertwines physical and operator heat,

\[
\boxed{
\mathcal Q_c(C^2u)=\Delta_{op}\mathcal Q_c(u),
\qquad
\sum_j\|[D_j,\mathcal Q_c(u)]\|_{HS}^2=\frac{M_3}{64}.
}
\]

Hence its Hilbert--Schmidt energy balance is exactly the critical Navier--Stokes balance:

\[
\boxed{
\partial_t\mathcal Q_c
=\mathcal Q_c(F_E)-\nu\Delta_{op}\mathcal Q_c,
\qquad
\frac12\frac d{dt}\|\mathcal Q_c\|_{HS}^2
=\frac{\kappa(0)}{64}-\nu\frac{M_3}{64}.
}
\]

If the positive eigenvalues of `mathscr R_c` are `lambda_j`, helicity oddness gives paired
spectrum and

\[
\boxed{K=128\nu^2\sum_j\lambda_j^2.}
\]

Therefore

\[
\boxed{
\|\mathscr R_c\|_{op}\le\frac{\sqrt K}{8\sqrt2\,\nu},
\qquad
N_{>1}^+\le\frac{K}{128\nu^2}.
}
\]

In particular `K<128 nu^2` implies `||mathscr R_c||<1` and `K'<=0` on every smooth continuum
interval where the identities apply.  This is a small-critical-data radius, not a large-data
regularity theorem.

For `T:u\mapsto\mathcal Q_c(u)`, polarization gives

\[
\boxed{T^*T=\frac1{64}\Lambda,\qquad
u=64\Lambda^{-1}T^*\mathcal Q_c(u).}
\]

Thus the transform is lossless on the mean-zero critical state, but its eigenvalues are not:
self-adjoint helicity-odd image operators have vanishing odd spectral traces, while `kappa(0)` is
generally nonzero.  No eigenvalue-only cubic Riccati closure exists.

Finally,

\[
\boxed{
T(\Lambda^\alpha u)=\Delta_{op}^{\alpha/2}T(u),
\qquad
\|u\|_{\dot H^s}^2
=64\|\Delta_{op}^{s/2-1/4}\mathcal Q_c(u)\|_{HS}^2.
}
\]

So the entire Sobolev hierarchy is one operator heat scale and `H^{1/2}` is its plain HS energy
level.  Section 11 now identifies a still lower full-graded current representation from which the
critical two-way channel itself is generated.

---

## 11. One scalar unit, one current creator, one curvature and one critical channel

Put

\[
\alpha=u^\flat,
\qquad
\beta=d\alpha,
\qquad
\boxed{Q=\nu\delta+\iota_u,\quad Q^*=\nu d+\alpha\wedge.}
\]

The scalar unit is a canonical basepoint of the graded algebra:

\[
\boxed{Q^*1=\alpha,\qquad (Q^*)^21=\nu\beta.}
\]

Thus state and curvature are successive applications of one operation.  Exterior Leibniz gives

\[
\boxed{
Q^2=\nu(\beta\wedge)^*,
\qquad
(Q^*)^2=\nu\beta\wedge.
}
\]

Curvature is therefore both

\[
\boxed{
\beta=\frac1\nu(Q^*)^21,
\qquad
Q^2=\nu(\beta\wedge)^*.
}
\]

### 11.1 Current and motion are the next step of the same chain

Let

\[
\mathbb H=QQ^*+Q^*Q.
\]

Since `Q1=0`,

\[
\boxed{
[\mathbb H,Q^*]1
=Q(Q^*)^21
=\nu Q\beta
=\nu e.
}
\]

The rotational momentum law `alpha_t+dB=-e` therefore becomes, after Leray projection,

\[
\boxed{
\partial_t(Q^*1)
=-\frac1\nu P[\mathbb H,Q^*]1,
\qquad
Q_t=\iota_{u_t}=-\iota_{P(Q\beta)^\sharp}.
}
\]

State, curvature, current and motion of the current operator are one chain; no second forcing
object is required.

Energy is the simplest one-way reading.  Because `alpha` is divergence-free, the pressure
projection disappears in the pairing and

\[
\boxed{
\frac12E'
=-\frac1\nu\langle Q^*1,[\mathbb H,Q^*]1\rangle
=-\frac1\nu\|(Q^*)^21\|_2^2
=-\nu\|\beta\|_2^2.
}
\]

The monotone direction is therefore adjointness inside the current chain, not an external damping
principle.

### 11.2 Master one-way / two-way metric law

Let `G=f(C)>=0` be a positive curl functional metric on divergence-free one-forms, extended to
vorticity two-forms by the same Hodge functional calculus, and put

\[
K_G=\langle\alpha,G\alpha\rangle.
\]

The exact NS identity is

\[
\boxed{
\frac12K_G'
=-\frac1\nu\|G^{1/2}(Q^*)^21\|_2^2
-\frac1\nu
\left\langle
[Q^*,G]Q^*1,(Q^*)^21
\right\rangle.
}
\]

The first term is the one-way curvature loss.  Every possibility of motion against that loss is
contained in the single metric/current commutator.  For `G=I` it vanishes.  At criticality
`G=Lambda=|C|`, the entire two-way channel is `[Q^*,Lambda]`.  The contrast with helicity is
primitive: Euler preserves the signed operator `C`, while the positive modulus `|C|` can fail to
commute with current creation.

### 11.3 The normalized critical channel creates its own curvature derivative

Extend

\[
\Lambda=|D|=(-\Delta)^{1/2}
\]

coefficientwise to the full exterior algebra; on divergence-free one-forms it equals `|C|`.
Define

\[
\boxed{
\mathcal A:=[\Lambda^{-1},Q^*]
=[\Lambda^{-1},\alpha\wedge],
\qquad
[Q^*,\Lambda]=\Lambda\mathcal A\Lambda.
}
\]

This is the normalized failure of the actual current creator to respect critical curl polarity.
Since `(Q*)^2=nu beta wedge`, put

\[
\boxed{
\mathcal B
:=\frac1\nu[\Lambda^{-1},(Q^*)^2]
=[\Lambda^{-1},\beta\wedge].
}
\]

The commutator Leibniz rule gives immediately

\[
\boxed{\{\mathcal A,Q^*\}=\nu\mathcal B.}
\]

Thus one more application of the same creator turns the two-way critical channel into a curvature
commutator.

On continuum `R^3`, use the full graded Hilbert--Schmidt norm on exterior forms.  The common raw
translation integral is

\[
\boxed{
\int_{\mathbb R^3}
\left(\frac1{|p+e|}-\frac1{|p|}\right)^2dp
=4\pi.
}
\]

For completeness, angular integration reduces `I(1)/(2pi)` to

\[
\int_0^1\!\left[r\log\frac{1+r}{1-r}-4r+2\right]dr
+\int_1^\infty\!\left[r\log\frac{r+1}{r-1}-2\right]dr=1+1.
\]

The second integral follows from the antiderivative
`((r^2-1)/2) log((r+1)/(r-1))-r`; the first is elementary by the same expansion/integration by
parts.  Thus the constant is analytic, not a fitted quadrature value.

With `C_F^2=(2pi)^-3`, exterior creation by a one-form has graded multiplicity `4`, while creation
by a two-form has multiplicity `2`.  Therefore

\[
\boxed{
\|\mathcal A\|_{HS,gr}^2=\frac{2}{\pi^2}K,
\qquad
\|\mathcal B\|_{HS,gr}^2=\frac1{\pi^2}M_3.
}
\]

These are continuum identities; a hard periodic/Galerkin cutoff carries boundary weights and is
not assigned these constants exactly.

Let `D_j=-i partial_j` and `epsilon_j=dx^j wedge`.  In the translation variable `q=k-l`,
`beta_q=i q wedge alpha_q`, hence

\[
\boxed{
\mathcal B
=i\sum_{j=1}^3\epsilon_j[D_j,\mathcal A].
}
\]

The exact constants sharpen this to the Dirichlet identity

\[
\boxed{
\sum_j\|[D_j,\mathcal A]\|_{HS,gr}^2
=2\|\mathcal B\|_{HS,gr}^2
=\frac{2}{\nu^2}\|\{\mathcal A,Q^*\}\|_{HS,gr}^2.
}
\]

So **heat erosion of the critical two-way channel is exactly the square of its incompatibility
with its own NS current creator.**  The escape channel carries the derivative through which heat
acts on it.

### 11.4 Pair-area closes the Euler loop and produces the exact critical residual

Write

\[
\mathsf E:=\alpha\wedge,
\qquad R:=\Lambda^{-1},
\qquad
\boxed{\mathcal V:=\mathsf E R\mathsf E.}
\]

Since `mathsf E^2=0`,

\[
\boxed{\{\mathcal A,\mathsf E\}=0,
\qquad
\mathsf E\mathcal A=\mathcal V,
\qquad
\mathcal A\mathsf E=-\mathcal V.}
\]

On continuum `R^3`, the Riesz kernel of `R` is `(2 pi^2 |x-y|^2)^-1`, so

\[
\boxed{
\|\mathcal V\|_{HS,gr}^2
=\frac1{2\pi^4}
\iint\frac{|u(x)\times u(y)|^2}{|x-y|^4}\,dx\,dy.}
\]

With

\[
\langle X,Y\rangle_{pair}
:=\frac1{2\pi^2}
\iint\frac{X\cdot Y}{|x-y|^4},
\quad
A_u=u(x)\times u(y),
\quad
D\omega=\omega(x)-\omega(y),
\]

fractional-Laplacian polarization gives

\[
\boxed{
\kappa(0)=\langle A_u,D\omega\rangle_{pair},
\qquad
M_3=\|D\omega\|_{pair}^2.}
\]

The kernels of `mathcal V` and `mathcal B` represent `A_u` and `-D omega`; hence

\[
\boxed{
\|\mathcal V\|_{HS,gr}^2=\frac1{\pi^2}\|A_u\|_{pair}^2,
\qquad
\kappa(0)=-\pi^2\Re\langle\mathcal B,\mathcal V\rangle_{HS,gr}.}
\]

Since `mathcal V^*=mathcal A^* i_u`, the second identity is also

\[
\boxed{
\kappa(0)=-\pi^2\Re\operatorname{Tr}
(\mathcal B\mathcal A^*\iota_u).}
\]

Thus the missing critical cubic is an orientation/current loop, not an eigenvalue-only trace.
The exact critical balance becomes

\[
\boxed{
K'=\frac{\pi^2}{2\nu}
\left(\|\mathcal V\|_{HS,gr}^2-\|\mathcal G_c\|_{HS,gr}^2\right),
\qquad
\mathcal G_c:=\mathcal V+2\nu\mathcal B,}
\]

or equivalently

\[
\boxed{
K'=\frac1{2\nu}
\left(\|A_u\|_{pair}^2-\|A_u-2\nu D\omega\|_{pair}^2\right).}
\]

For `M_3>0`, the same endogenous coefficient from Section 9 is the unique pair least-squares
projection coefficient,

\[
\boxed{
\nu_E
=\frac{\langle A_u,D\omega\rangle_{pair}}{\|D\omega\|_{pair}^2}
=\arg\min_{\lambda\in\mathbb R}\|A_u-\lambda D\omega\|_{pair}^2.}
\]

Write the forced orthogonal decomposition

\[
\boxed{
A_u=\nu_E D\omega+C_\perp,
\qquad
\langle C_\perp,D\omega\rangle_{pair}=0.}
\]

Then

\[
\|A_u\|_{pair}^2=\nu_E^2M_3+\|C_\perp\|_{pair}^2,
\qquad
\|A_u-2\nu D\omega\|_{pair}^2
=(\nu_E-2\nu)^2M_3+\|C_\perp\|_{pair}^2,
\]

so the whole reconfiguration component cancels **exactly** from the instantaneous critical balance:

\[
\boxed{
\|A_u\|_{pair}^2-\|A_u-2\nu D\omega\|_{pair}^2
=4\nu M_3(\nu_E-\nu).}
\]

This cancellation forbids a static closure based only on making `C_perp` large.  Its role can only be
historical: it records the 3D incompatibility which turns the state away from the productive road.

The equality case is rigid.  If `C_perp=0`, then `A_u=nu_E D omega` is a two-point coboundary, hence
for almost every triple

\[
0=A_u(x,y)+A_u(y,z)+A_u(z,x)
=[u(x)-u(z)]\times[u(y)-u(z)].
\]

Thus the essential velocity image lies on one affine line: `u=a+phi b`.  Incompressibility gives
`b dot grad phi=0`; after removing the uniform vector `a` by a Galilean frame this is a shear with
zero Euler self-interaction.  On finite-energy `R^3` a nonzero such shear is excluded by its invariant
direction.  Therefore exact perfect pair alignment is lower-dimensional/null geometry, but this is
only equality rigidity -- no uniform quantitative gap away from it is asserted.

So `K'>0` exactly when the state-selected Euler coefficient beats molecular viscosity, while the
orthogonal 3D component controls where that productive direction can persist next.

### 11.5 The residual deciding `K'` is itself a curved current

The residual law extends without adding a new state variable.  For every `lambda>0`, use the two
connections generated by the same `d` and `alpha`,

\[
D_R^{(\lambda)}:=\lambda d,
\qquad
D_L^{(\lambda)}:=\lambda d+\mathsf E,
\]

with induced Hom derivative `nabla_lambda`.  Put

\[
\boxed{\mathcal J_\lambda:=\mathcal V+\lambda\mathcal B.}
\]

Because `{d,mathcal A}=mathcal B`, `mathsf E mathcal A=mathcal V`, and
`(D_L^{(lambda)})^2=lambda mathsf F`, one has exactly

\[
\boxed{
\mathcal J_\lambda=\nabla_\lambda\mathcal A,
\qquad
\nabla_\lambda\mathcal J_\lambda
=\lambda\mathsf F\mathcal A.}
\]

Hence

\[
\boxed{\mathcal J_\lambda=0\Longrightarrow \mathsf F\mathcal A=0\qquad(\lambda>0).}
\]

Under the physical kernel identification, `mathcal J_lambda` is precisely
`A_u-lambda D omega`.  Thus in the productive regime `nu_E>0`, perfect productive alignment would
make `mathcal J_{nu_E}=0` and force the same dimensional-collapse source `mathsf F mathcal A=0`.

The physical residual is the member `lambda=2nu`.  Indeed if
`nabla=(1/2)nabla_{2nu}`, then

\[
\boxed{
\mathcal G_c=\mathcal J_{2\nu}=2\nabla\mathcal A,
\qquad
\nabla\mathcal G_c=\nu\mathsf F\mathcal A,}
\]

which is exactly the previous midpoint normalization.  Thus the family reveals rigidity of the
same current; it does not introduce an analyst-selected family of mechanisms.

This is the critical Gauss--Bianchi law on **the exact residual which decides `K'`**:

\[
\boxed{
\mathcal A\xrightarrow{\nabla}\mathcal G_c
\xrightarrow{\nabla}\nu(\beta\wedge)\mathcal A.}
\]

The kernels are explicit,

\[
\boxed{
\mathcal A(x,y)
=\frac{[u(y)-u(x)]^\flat\wedge}{2\pi^2|x-y|^2},}
\]

\[
\boxed{
(\mathsf F\mathcal A)(x,y)
=\frac{\omega(x)\cdot[u(y)-u(x)]}{2\pi^2|x-y|^2}\,dV,}
\]

and therefore

\[
\boxed{
\|\mathsf F\mathcal A\|_{HS}^2
=\frac1{4\pi^4}
\iint\frac{|\omega(x)\cdot[u(y)-u(x)]|^2}{|x-y|^4}\,dx\,dy.}
\]

The source is precisely local vorticity piercing the affine geometry of the velocity values.

### 11.6 Null geometry, six-dimensional coercivity and the sharpened gap

If `mathsf F mathcal A=0`, then

\[
\boxed{\omega(x)\cdot[u(y)-u(x)]=0\quad\forall x,y.}
\]

The resulting null geometry is classified more primitively by the positive critical
carre-du-champ tensor in Section 11.7 below.

For a top Hom test `Phi=phi(x,y)dV`,

\[
\langle\nu\mathsf F\mathcal A,\Phi\rangle
=\langle\mathcal G_c,\nabla^*\Phi\rangle.
\]

The two pieces of `nabla*` lie in orthogonal Hodge channels; incompressibility removes the
integrated drift cross term.  Thus, in quadratic-form sense,

\[
\boxed{
\mathcal L_u^{(6)}:=\nabla\nabla^*
=-\nu^2(\Delta_x+\Delta_y)+\frac{|u(x)|^2}{4}.}
\]

Whenever the inverse is defined,

\[
\boxed{
\|\mathcal G_c\|_{HS,gr}^2
\ge
\left\langle
\nu\mathsf F\mathcal A,
(\mathcal L_u^{(6)})^{-1}\nu\mathsf F\mathcal A
\right\rangle.}
\]

The six variables are only the two physical points already forced by `|D|`.  Combining this with
the critical square gives

\[
\boxed{
K'\le\frac{\pi^2}{2\nu}
\left[
\|\mathcal V\|_{HS,gr}^2
-
\langle\nu\mathsf F\mathcal A,
(\mathcal L_u^{(6)})^{-1}\nu\mathsf F\mathcal A\rangle
\right].}
\]

The earlier self-derivative law remains exact,

\[
\boxed{
\frac12\frac d{dt}\|\mathcal A\|_{HS,gr}^2
=\frac{2}{\pi^2}\kappa(0)
-\frac2\nu\|\{\mathcal A,Q^*\}\|_{HS,gr}^2,}
\]

with zero-set rigidity

\[
\boxed{\{\mathcal A,Q^*\}=0\Longrightarrow M_3=0\Longrightarrow u=0.}
\]

What is **not** automatic is quantitative stability:

\[
\|\mathsf F\mathcal A\|_{(\mathcal L_u^{(6)})^{-1}}\ll1
\Longrightarrow
u\ \text{near the intrinsic 2D/shear/null manifold}.
\]

That implication is a separate stability/compactness theorem.  The historical gap is then to show
that a genuinely 3D trajectory cannot stay close enough to this null geometry for the infinite
productive critical action required by escape while still regenerating `mathcal A` nonintegrably.

### 11.7 The critical carre-du-champ is the positive metric behind the Gauss source

For scalar functions let

\[
\Gamma(f,g):=f\Lambda g+g\Lambda f-\Lambda(fg),
\]

and for the velocity define the matrix `Gamma_u=(Gamma(u_i,u_j))`.  On continuum `R^3`,

\[
\boxed{
\Gamma_u(x)=\frac1{\pi^2}\int
\frac{[u(x)-u(y)]\otimes[u(x)-u(y)]}{|x-y|^4}\,dy\succeq0.}
\]

Thus criticality is already an intrinsic positive increment metric:

\[
\boxed{\int\operatorname{tr}\Gamma_u\,dx=2K.}
\]

The Gauss source from Section 11.5 is exactly vorticity measured in this same metric,

\[
\boxed{
\|\mathsf F\mathcal A\|_{HS}^2
=\frac1{4\pi^2}\int \omega^T\Gamma_u\omega\,dx.}
\]

It is also the finite-difference completion of vortex stretching.  For `y=x+h`,

\[
\boxed{
\omega(x)\cdot[u(x+h)-u(x)]
=(S\omega)(x)\cdot h+O(|h|^2),}
\]

because the antisymmetric part of `nabla u` annihilates `omega`.  The same source whose curvature
controls the critical residual therefore has tangent `S omega`, the 3D Euler stretching vector.

Let `V_u=span{u(y)-u(z)}`.  Positivity gives

\[
\boxed{\ker\Gamma_u(x)=V_u^\perp}
\]

for almost every `x` (with the usual essential-span interpretation).  Hence zero Gauss source has
exactly the previous lower-dimensional alternatives: full span forces `omega=0`; rank two is a
fixed embedded 2D flow; rank one is a shear `u=phi b`, `b dot grad phi=0`; rank zero is uniform.
On finite-energy continuum `R^3` every nonzero 2D/shear/uniform representative is excluded by its
translation-invariant direction, while curl-free divergence-free `L^2` fields vanish.  Therefore

\[
\boxed{
\int\omega^T\Gamma_u\omega\,dx=0\quad\Longrightarrow\quad u=0
\qquad (u\in L^2(\mathbb R^3)).}
\]

This is injectivity, not a uniform coercive gap.  On periodic geometry the nonzero null classes do
exist, but they are invariant under NS: embedded 2D remains in its fixed plane, while shear has
`F_E=0` and evolves only by heat.  A trajectory cannot use exact null geometry to rotate its plane
and restart a 3D critical channel.

The rank itself is intrinsic.  Writing `delta u_a=u(x)-u(y_a)`, Cauchy--Binet gives

\[
\boxed{
e_2(\Gamma_u)=\frac1{2\pi^4}\iint
\frac{|\delta u_1\times\delta u_2|^2}
{|x-y_1|^4|x-y_2|^4}\,dy_1dy_2,}
\]

\[
\boxed{
\det\Gamma_u=\frac1{6\pi^6}\iiint
\frac{\det(\delta u_1,\delta u_2,\delta u_3)^2}
{|x-y_1|^4|x-y_2|^4|x-y_3|^4}\,dy_1dy_2dy_3.}
\]

Thus `det Gamma_u=0` exactly when the essential increment span is planar, while `e_2(Gamma_u)=0`
exactly when it is one-dimensional.  No best plane or analyst-defined 3D coefficient is needed.

The tensor also carries the literal one-way/two-way split.  Put `L=Lambda^2=-Delta`.  Since the
Laplacian product defect is `-2 grad f dot grad g`, the full NS law is

\[
\boxed{
(\partial_t+\nu L)\Gamma_u
=\Gamma(F_E,u)+\Gamma(u,F_E)
-2\nu\sum_{j=1}^3\Gamma_{\partial_j u}.}
\]

Every `Gamma_{partial_j u}` is positive semidefinite.  Heat is therefore a one-way Loewner sink of
the critical metric, while Euler is the only two-way regeneration/reorientation term.  Integrating
the trace recovers only the scalar shadow

\[
\boxed{K'=2\kappa(0)-2\nu M_3.}
\]

On the exact critical null geometry `S omega=0`, so enstrophy also becomes one-way; on shear
`F_E=0`, every Sobolev quadratic is heat-only.  This gives the intrinsic dimensional ladder

\[
\boxed{
\text{3D two-way criticality}
\to \omega\in\ker\Gamma_u
\to \text{dimensional collapse}
\to \text{a stronger one-way law}.}
\]

The tangent identity does **not** supply the closure.  A physical-Fourier referee rejects the
proposed sharp factor `1/2` in
`||S omega||_{Hdot^-1/2}^2 <= C int omega^T Gamma_u omega`; a refined search has approached `C=1`
from below (about `0.99799`) without proving it.  A separate numerical referee also rejects the
tempting operator shortcut `A(F_E)=-2 d_op^* V` (order-one residual), while differentiating the
affine defect simply climbs the old `H^{3/2}->H^{5/2}->...` ladder.  None of these shortcuts is used
below; the exact pair dynamics makes the remaining obstruction a history problem instead.

### 11.8 The two-particle law lies below the critical metric

Put

\[
x=c+\frac r2,\qquad y=c-\frac r2,\qquad
U:=\frac{u(x)+u(y)}2,\qquad v:=u(x)-u(y),
\]

and `pbar=(p(x)+p(y))/2`, `delta p=p(x)-p(y)`.  These are the common and relative
motions of the two actual endpoints.  At fixed physical separation `r`, NS itself splits exactly as

\[
\boxed{
\partial_tU+(U\cdot\nabla_c)U+\frac14(v\cdot\nabla_c)v+\nabla_c\bar p
=\nu\Delta_cU,}
\]

\[
\boxed{
\partial_tv+(U\cdot\nabla_c)v+(v\cdot\nabla_c)U+\nabla_c\delta p
=\nu\Delta_cv.}
\]

The endpoint origin forces the compatibility identities

\[
\boxed{
\nabla_rv=\nabla_cU,\qquad
4\nabla_rU=\nabla_cv,\qquad
\Delta_cv=4\Delta_rv.}
\]

Thus center motion and separation motion do not carry independent heat operators.  If
`D_pair=partial_t+U dot nabla_c+v dot nabla_r`, then the second equation is equivalently
`D_pair v+nabla_c delta p=nu Delta_c v`.

Translation invariance gives, for every fixed `r`,

\[
\boxed{\langle U,v\rangle_{L^2_c}=0,\qquad
\|U\|_2^2+\frac14\|v\|_2^2=E.}
\]

Writing

\[
T(r):=\int (v\otimes v):S_U\,dc,
\]

the two roads have the exact energy exchange

\[
\boxed{
\frac12\frac d{dt}\|v\|_2^2=-T(r)-\nu\|\nabla_cv\|_2^2,
\qquad
\frac12\frac d{dt}\|U\|_2^2=\frac14T(r)-\nu\|\nabla_cU\|_2^2.}
\]

Euler only exchanges common and relative endpoint energy; heat only removes it.

For the relative density `q=|v|^2/2`, subtracting the endpoint equations also gives the exact local
pair conservation law

\[
\boxed{
\partial_tq
+\nabla_c\!\cdot(qU+\delta p\,v)
+\nabla_r\!\cdot(qv)
=\nu\left(\frac12\Delta_c+2\Delta_r\right)q
-\nu\big(|\nabla u(x)|^2+|\nabla u(y)|^2\big).}
\]

Pressure therefore has only the **common-coordinate** flux in total pair energy.  After integrating
`c`, define

\[
Q(r,t):=\int q(c,r,t)\,dc,\qquad
J(r,t):=\int q(c,r,t)v(c,r,t)\,dc.
\]

Then pressure disappears:

\[
\boxed{\partial_tQ+\nabla_r\cdot J=2\nu\Delta_rQ-2\nu Z.}
\]

On continuum `R^3`, critical stock, critical heat and Euler work are all readings of this same
relative field.  Put

\[
\boxed{W(c,r,t):=\frac{v(c,r,t)}{|r|^2}.}
\]

Then

\[
\boxed{
K=\frac1{2\pi^2}\|W\|_{L^2_{c,r}}^2,
\qquad
M_3=\frac1{2\pi^2}\|\nabla_cW\|_{L^2_{c,r}}^2,}
\]

while

\[
\boxed{
K=\frac1{\pi^2}\int\frac{Q(r)}{|r|^4}\,dr,
\qquad
\kappa(0)=-\frac2{\pi^2}\int\frac{J(r)\cdot r}{|r|^6}\,dr.}
\]

Thus positive `kappa(0)` is weighted net inward relative transport toward `r=0`.  If
`lambda=-(v dot r)/|r|^2` is logarithmic pair compression, the material-pair equation becomes

\[
\boxed{
D_{pair}W
=2\lambda W-\frac{\nabla_c\delta p}{|r|^2}+\nu\Delta_cW.}
\]

Pair advection and pressure have zero global `L^2` work, so this single field gives directly

\[
\boxed{
\frac12\frac d{dt}\|W\|_2^2
=2\iint\lambda|W|^2-\nu\|\nabla_cW\|_2^2,}
\]

which is exactly `K'=2 kappa(0)-2 nu M_3`.  Thus critical heat is the **center-Dirichlet cost of the
inverted-pair velocity**; `D I(r)` differs from `W` only by the radial reflection
`I-2 nhat tensor nhat`.

For `mathcal I(r)=r/|r|^2`,
`D mathcal I=|r|^-2(I-2 nhat tensor nhat)` is a scaled reflection and
`D mathcal I(r)v dot D^2 mathcal I(r)[v,v]=-2|v|^2(v dot r)/|r|^6`; its pressure work integrates to
zero.  Hence material labels `a,b`, `r_ab=Phi_t(a)-Phi_t(b)`, give exactly

\[
\boxed{
K(t)=\frac1{2\pi^2}\iint
\left|\frac d{dt}\mathcal I(r_{ab}(t))\right|^2\,da\,db.}
\]

Moreover `K^2<=EZ` and `-E'=2 nu Z` imply the genuine finite-time history bound

\[
\boxed{
\|\mathcal I(r_{ab}(T))-\mathcal I(r_{ab}(0))\|_{L^2_{a,b}}
\le \sqrt{2\pi^2}\,E(0)^{1/4}T^{3/4}
\left[\frac{E(0)-E(T)}{2\nu}\right]^{1/4}.}
\]

This is finite Hilbert path length, **not** endpoint-velocity control or regularity.

### 11.9 Incompressibility fixes the radial/transverse critical split

The separation itself supplies the only further split.  Write

\[
n=\frac r{|r|},\qquad \delta u=a n+b,\qquad b\perp n.
\]

With the same critical kernel define

\[
K_\parallel=\frac1{2\pi^2}\iint\frac{a^2}{|r|^4},\qquad
K_\perp=\frac1{2\pi^2}\iint\frac{|b|^2}{|r|^4}.
\]

Plancherel plus `k dot u_hat(k)=0` gives the exact continuum identity

\[
\boxed{K_\parallel=\frac14K,\qquad K_\perp=\frac34K.}
\]

Indeed the radial `|r|` integration leaves the angular weight `|khat dot n|`; for a transverse
Fourier polarization the longitudinal angular fraction is exactly `1/4`.  Thus incompressible 3D
critical motion cannot be pure pair collapse.  The relative velocity is itself separation-space
incompressible,

\[
\boxed{\nabla_r\cdot\delta u=\tfrac12\operatorname{div}u(x)+\tfrac12\operatorname{div}u(y)=0.}
\]

Put

\[
A_3:=\iint\frac{a^3}{|r|^5},\qquad
B_3:=\iint\frac{a|b|^2}{|r|^5}.
\]

Then `kappa(0)=-pi^-2(A_3+B_3)`.  Since the `1:3` split holds at every time,
`K'_{parallel,E}=kappa(0)/2`.  Differentiating `a^2|r|^-4` along the material pair law gives the
exact pressure exchange

\[
\boxed{
P_\parallel=\frac{3}{2\pi^2}(A_3-B_3),\qquad
P_\perp=-P_\parallel.}
\]

There is no conflict with the preceding center-road law: pressure has no direct `r`-flux for the
total pair energy, but inside the critical kinetic metric it can convert radial and transverse
motion with zero net critical work.  Under inversion the same split is literal kinetic geometry,
since `D mathcal I(r)` only reflects the radial direction and multiplies both channels by `|r|^-2`.

The primitive grammar is therefore shorter: relative advection supplies the two-way inward/outward
critical work; pressure is a radial/transverse converter; viscosity removes total critical kinetic
energy one way.  Exact `K_perp=3K_parallel` is a compulsory transverse **stock** identity, not yet a
time-action theorem.

### 11.10 Affine defect and the collision boundary

The common/relative system already contains its own relay variable as a derivative, not a new
state.  Put

\[
A:=\nabla_cU=\frac{\nabla u(x)+\nabla u(y)}2,
\qquad
\boxed{C:=\nabla_cv=\nabla u(x)-\nabla u(y).}
\]

Thus `C` is exactly center non-affinity of the relative field and

\[
\boxed{M_3=\frac1{2\pi^2}\iint\frac{|C(c,r)|^2}{|r|^4}\,dc\,dr.}
\]

Subtract the two gradient-NS equations.  Since
`tr[(A+C/2)^2]-tr[(A-C/2)^2]=2 tr(AC)` and
`-Delta p=tr((nabla u)^2)`, with
`R_c:=nabla_c^2(-Delta_c)^-1`, one obtains

\[
\boxed{
D_{pair}C+AC+CA+2\mathcal R_c[\operatorname{tr}(AC)]
=\nu\Delta_cC.}
\]

The whole equation is homogeneous in `C`: pressure can redistribute affine defect but has no source
independent of it.  Global `C=0` means `nabla u(x)=nabla u(y)` for every pair, hence `u=Ax+b`.
Nontrivial such fields have infinite energy on `R^3`; the affine blowup countermodel is therefore the
exact finite-dimensional null geometry that the finite-energy problem must leave by paying `C`.

There is also no fixed-separation escape.  The exact road energy gives
`||v(.,r)||_2^2<=4E`, hence for every `R>0`

\[
\boxed{K_{|r|\ge R}\le\frac{8E}{\pi R}.}
\]

Consequently, if `K(t_n)->infinity`, then for every fixed `R>0`

\[
\boxed{\frac{K_{|r|\ge R}(t_n)}{K(t_n)}\longrightarrow0.}
\]

The normalized critical mass must concentrate on the physical collision diagonal `r=0`; no shell
selector is needed.  For smooth flow define its radial density by

\[
\mathfrak k(\rho,t):=\frac1{2\pi^2\rho^2}
\int_{S^2}\!\int |v(c,\rho n,t)|^2\,dc\,d\Omega,
\qquad K=\int_0^\infty\mathfrak k(\rho,t)\,d\rho.
\]

The symmetric endpoint expansion `v(c,r)=nabla u(c) r+O(|r|^3)` yields

\[
\boxed{\mathfrak k(0,t)=\frac{2}{3\pi}Z(t),}
\]

and therefore the boundary evolution is exactly the enstrophy law,

\[
\boxed{
\partial_t\mathfrak k(0,t)
=\frac{4}{3\pi}\int\omega\cdot S\omega
-\frac{4\nu}{3\pi}\|\operatorname{curl}\omega\|_2^2.}
\]

Thus critical transport toward `r=0` does not enter an unknown boundary dynamics.  It lands on the
existing vortex-stretching/heat law.  But concentration at the diagonal is only a **necessary
location statement**, not critical growth by itself.  Under the exact NS dilation

\[
u_\lambda(x,t)=\lambda u(\lambda x,\lambda^2t),
\]

\[
\boxed{E_\lambda=\lambda^{-1}E,\quad Z_\lambda=\lambda Z,\quad
K_\lambda=K,\quad M_{3,\lambda}=\lambda^2M_3,\quad
\kappa_\lambda=\lambda^2\kappa.}
\]

So a fixed critical profile may shrink toward `r=0` while carrying the **same** `K`.  Blowup of `K`
requires regeneration/accumulation of additional critical mass, not scale motion alone.  A
scale-covariant critical unit of size `rho` has the benchmark

\[
E_\rho\sim\rho,\qquad Z_\rho\sim\rho^{-1},\qquad
M_{3,\rho}\sim\rho^{-2},\qquad \Delta t_\rho\sim\rho^2.
\]

Hence its physical-energy bill over one parabolic lifetime is `~nu rho`, while its critical-heat
bill `nu M_3 Delta t` is scale-neutral (`~nu`).  This is a scaling benchmark, **not** a universal
per-visit lower bound.

The state-selected coefficient makes the immediate necessary history law completely transparent:

\[
\frac d{dt}\log K
=2\frac{M_3}{K}(\nu_E-\nu).
\]

Therefore

\[
\boxed{
K(t)\to\infty
\Longrightarrow
\int^T\frac{M_3}{K}(\nu_E-\nu)_+\,dt=\infty.}
\]

The sharper invariant currency already proved above is, with `N^2=Z/E`,

\[
\boxed{
K/E\to\infty
\Longrightarrow
\int^T\frac{\kappa(0,t)^2}
{N(t)^2[E(t)Z(t)-K(t)^2]}\,dt=\infty.}
\]

A useful falsifier shows why energy and instantaneous rigidity cannot close this by themselves.
Consider only a hypothetical affine-like core, **not an NS solution**, with `tau=T-t`,
`|grad u|~tau^-1` and radius `R~tau^alpha`.  Pure scaling gives

\[
E_{core}\sim\tau^{-2+5\alpha},\quad
K\sim\tau^{-2+4\alpha},\quad
Z\sim\tau^{-2+3\alpha},\quad
M_3\sim\tau^{-2+2\alpha},\quad
\kappa\sim\tau^{-3+4\alpha}.
\]

For

\[
\boxed{\frac25\le\alpha<\frac12}
\]

this geometry has bounded core energy and finite `int Z dt`, yet `K` diverges and
`kappa/(nu M_3)~tau^{2alpha-1}` diverges.  It is not a counterexample to NS; it is a countergeometry
to any proof using only the energy law, the instantaneous heat bill, and local affine rigidity.
Global dynamical compatibility must rule it out.

The latest falsification layer points the same way rather than opening another branch.  The static
`1/2` stretching/`Gamma_u` bridge is false and factor `1` remains unproved; the simple
`A(F_E)=-2 d_op^* V` codifferential closure has an order-one numerical residual; differentiating the
affine defect only climbs the Sobolev ladder; actual Galerkin states show that critical growth does
not force enstrophy growth; and a high-radius `H^{-1}_{curl}` state-speed referee makes the required
scale-free factor grow beyond `5e4`.  The exact low-frequency Galilean catalyst suppression also
does not close many-mode coherence, because multiple admissible contributions may add before the
square.  These are **guards**, not theorem inputs: they remove static shortcuts and leave persistence.

The remaining theorem is therefore the finite-persistence statement that the actual NS current
cannot regenerate `nu_E` productively often enough to make the scale-invariant action diverge.  The
exact equality rigidity says perfect productivity destroys genuine 3D geometry; the open step is to
exclude arbitrarily near-perfect productive returns **infinitely often** on finite physical history.
No global-regularity claim is made here.

---

## 12. The full `Q`-Hodge square contains the velocity gradient and regenerates Lamb

Define

\[
\boxed{\mathbb H:=Q^*Q+QQ^*,\qquad L_u^Q:=-\nu^2\Delta+|u|^2.}
\]

On one-forms and, after Hodge identification, on two-forms,

\[
\boxed{
\mathbb H^{(1)}=L_u^Q I+2\nu S,\qquad
\star^{-1}\mathbb H^{(2)}\star=L_u^Q I-2\nu S.
}
\]

Hence

\[
\boxed{\mathbb H^{(1)}-\star^{-1}\mathbb H^{(2)}\star=4\nu S.}
\]

The antisymmetric part of `nabla u` is encoded by `Q^2`; the symmetric part is the adjacent-degree imbalance of `mathbb H`.  The full velocity gradient is therefore generated by `Q,Q^*`.

### 12.1 Lamb is the adjacent-degree intertwining defect

Associativity gives the graded commutator identity

\[
\boxed{\mathbb H Q-Q\mathbb H=Q^*Q^2-Q^2Q^*=[Q^*,Q^2].}
\]

Let `eta=star b^flat` be closed, equivalently `div b=0`.  Then

\[
Q^2\eta=\nu\,\omega\cdot b,\qquad
Q^*\eta=(u\cdot b)dV,
\]

and therefore

\[
\boxed{
\mathbb H^{(1)}(Q\eta)-Q(\mathbb H^{(2)}\eta)
=\nu^2d(\omega\cdot b)+\nu[\,b\times(u\times\omega)\,]^\flat.
}
\]

The exact first term is a gradient.  The second is the curvature-transverse Lamb interaction.  For `b=omega`,

\[
\boxed{
\mathbb H^{(1)}(Q\beta)-Q(\mathbb H^{(2)}\beta)
=\nu^2d|\omega|^2+\nu[\omega\times(u\times\omega)]^\flat.
}
\]

If

\[
\Xi:=\bigl(\mathbb H^{(1)}(Q\beta)-Q(\mathbb H^{(2)}\beta)-\nu^2d|\omega|^2\bigr)^\sharp,
\]

then, wherever `omega != 0`,

\[
\boxed{u\times\omega=-\frac{\omega\times\Xi}{\nu|\omega|^2}.}
\]

Thus Lamb is not an independent forcing object.  It is reconstructed from the failure of `Q` to intertwine the adjacent Hodge degrees generated by its own curvature.  The two primitive feedback loops are

\[
\boxed{Q\to Q^2\to\beta\to Q\beta\to Q_t,}
\qquad
\boxed{Q\to\mathbb H\to S\to[\mathbb H,Q]\to\text{Lamb interaction}.}
\]

---

## 13. The midpoint `M` is canonical, but not fundamental

Critical/Poynting reflection selects the midpoint between pure Hodge heat and the actual current:

\[
\boxed{
M:=\nu\delta+\frac12\iota_u=\frac12(\nu\delta+Q),\qquad
M^*=\nu d+\frac12\alpha\wedge.
}
\]

It obeys

\[
\boxed{M^2=\frac\nu2(\beta\wedge)^*.}
\]

For

\[
L_u^M=-\nu^2\Delta+\frac{|u|^2}{4},\qquad \mathbb M=M^*M+MM^*,
\]

one has

\[
\boxed{
\mathbb M^{(1)}=L_u^M I+\nu S,\qquad
\star^{-1}\mathbb M^{(2)}\star=L_u^M I-\nu S,
}
\]

so `2nu S` is their degree imbalance, while

\[
\boxed{MM^*|_{\Omega^0}=-\nu^2\Delta+\frac{|u|^2}{4}}
\]

is exactly the scalar Schrödinger normal operator.

### 13.1 Poynting/Gauss is the midpoint self-return

On vorticity,

\[
\boxed{
(M\beta)^\sharp=\nu\operatorname{curl}\omega-rac12u\times\omega=-\frac12G,
\qquad
G=u\times\omega-2\nu\operatorname{curl}\omega.
}
\]

But the same midpoint maps that current back to positive curvature:

\[
\boxed{M(M\beta)=M^2\beta=\frac\nu2|\beta|^2,}
\]

or equivalently

\[
\boxed{
\left(\nabla\cdot-\frac{u}{2\nu}\cdot\right)G=|\omega|^2.
}
\]

Thus the current which critical growth tries to make small is not an arbitrary error.

### 13.2 Mixed Gauss law and the curvature floor

For any smooth vector field `b`,

\[
\boxed{
\mathfrak G_\nu[u;b]=u\times b-2\nu\operatorname{curl}b,\qquad
\left(\nabla\cdot-\frac{u}{2\nu}\cdot\right)\mathfrak G_\nu[u;b]=\omega\cdot b.
}
\]

For divergence-free `b`,

\[
\boxed{
\nu\langle b,Sb\rangle
=\nu^2\|Cb\|_2^2+rac14\|u\times b\|_2^2-rac14\|\mathfrak G_\nu[u;b]\|_2^2.
}
\]

Hence, with `A_u=PS` and `U_u b=u cross b`,

\[
\boxed{4\nu A_u=P(4\nu^2C^2+U_u^*U_u-\mathfrak G_u^*\mathfrak G_u)P.}
\]

The helicity-odd block is a difference of Gram currents, not a negative square.  If

\[
D_u=\nabla\cdot-\frac{u}{2\nu}\cdot,\qquad M_\omega b=\omega\cdot b,
\]

then `D_u mathfrak G_u=M_omega` and

\[
\boxed{D_uD_u^*=-\Delta+\frac{|u|^2}{4\nu^2}.}
\]

Therefore, wherever the inverse is defined,

\[
\boxed{
\mathfrak G_u^*\mathfrak G_u
\ge M_\omega^*\left(-\Delta+\frac{|u|^2}{4\nu^2}\right)^{-1}M_\omega.
}
\]

---

## 14. The midpoint gives the pressure-free critical tax, and static coercivity stops there

Let

\[
F_E=P(u\times\omega),\qquad
j_0=\nu C^2u,\qquad
j_1=\nu C^2u-F_E=-u_t.
\]

Their midpoint is

\[
\boxed{j_{1/2}=\frac{j_0+j_1}{2}=\nu C^2u-\frac12F_E=-\frac12PG.}
\]

For every `K_s=||Lambda^s u||_2^2`,

\[
\boxed{
K_s'=\frac1{2\nu}
\left(\|\Lambda^{s-1}F_E\|_2^2-\|\Lambda^{s-1}PG\|_2^2\right).
}
\]

At `s=1/2`, Euler energy nullity implies

\[
\boxed{
\|\Lambda^{-1/2}PG\|_2^2\ge\frac{4\nu^2Z^2}{K},
}
\]

and therefore

\[
\boxed{
K'\le\frac1{2\nu}\|\Lambda^{-1/2}F_E\|_2^2-2\nu\frac{Z^2}{K}.
}
\]

Helicity nullity supplies the second invariant-normal projection and the sharpened two-null Gram-plane tax recorded in the theorem certificate; no new stock is introduced.

Static Gauss coercivity cannot close the projected gap.  For

\[
q_\theta^*\phi=\nu d\phi+\theta\alpha\phi,
\]

requiring this scalar test to be divergence-free gives

\[
\nu\Delta\phi+\theta u\cdot\nabla\phi=0.
\]

On periodic/decaying domains incompressibility yields

\[
\boxed{\phi=\text{constant}.}
\]

Thus every nonconstant Schrödinger refinement necessarily reads the exact/pressure sector.

---

## 15. Material Hodge geometry measures the same turning in another gauge

Let `Phi_t` be the incompressible flow and set

\[
\widetilde\beta=\Phi_t^*\beta,\qquad g=\Phi_t^*g_0,\qquad L_g=d\delta_g+\delta_gd.
\]

Then

\[
\boxed{\partial_t\widetilde\beta=-\nu L_g\widetilde\beta,}
\]

while Euler moves the Hodge frame:

\[
\boxed{
\partial_tL_g=[\mathcal L_v,L_g],\qquad
\partial_t\delta_g=[\mathcal L_v,\delta_g].
}
\]

For `c=delta_g tilde beta`,

\[
\boxed{
\partial_t c=[\mathcal L_v,\delta_g]\widetilde\beta-\nu L_g^{(1)}c.
}
\]

With deformation gradient `F=D_a Phi`, Maurer--Cartan fields

\[
\Gamma=F^{-1}dF,\qquad B=F^{-1}F_t
\]

obey

\[
\boxed{d\Gamma+\Gamma\wedge\Gamma=0,\qquad \partial_t\Gamma=D_\Gamma B,}
\]

and the turnover speed is locked to the Hodge current:

\[
\boxed{
\|D_\Gamma B\|_{L^2_g}^2=\|\delta_g\widetilde\beta\|_{L^2_g}^2,\qquad
\|\nabla^g g_t\|_{L^2_g}^2=2\|\delta_g\widetilde\beta\|_{L^2_g}^2.
}
\]

The separate material theorem adds the transverse determinant and Minkowski memory laws.  Thus critical/current reorientation and material heat memory are two gauges of the same state-generated turning, not independent escape mechanisms.

---

## 16. Falsification guards after the `Q` collapse

1. **No global-regularity claim.**  Concentrated many-mode states can still produce large critical action and rapid reorientation.
2. **No static Gauss/Schrödinger closure.**  Nonconstant scalar tests enter the pressure/exact sector.
3. **No negative square after helicity projection.**  The critical odd block is a Gram difference with no fixed inherited sign.
4. **No moving-projector proof without spectral-gap information.**  Reynolds eigenspaces remain a representation, not the ontology.
5. **No free `L^4_omega` damping.**  The apparent quartic Gauss term cancels in the exact current-energy calculation.
6. **No higher-Sobolev control from `Q_t^2=0`.**  Nilpotent tangents do not prevent high-frequency concentration.
7. **No spectral-only monotone beyond energy/helicity**, no low-Krylov closure and no automatic viscous dephasing.
8. **No eigenvalue-only Riccati law** and no unrestricted Nambu--Poisson theorem; the surviving structure is lossless Cartan/Jacobi current geometry.
9. **No local finite-dimensional closure.**  The missing theorem must couple full spatial current persistence to material Hodge memory.
10. **No uniform gap from the new rigidity.**  `{mathcal A,Q*}=0 => u=0` is an exact zero-set statement; it does not bound the ratio between regeneration and incompatibility away from saturation on arbitrary large states.

---

## 17. Minimal ontology and the remaining primitive theorem

The primitive state algebra is

\[
\boxed{Q=\nu\delta+\iota_u,
\qquad Q^*=\nu d+u^\flat\wedge,}
\]

\[
\boxed{
1\xrightarrow{Q^*}\alpha\xrightarrow{Q^*}\nu\beta,
\qquad
Q^2=\nu(\beta\wedge)^*,
\qquad
[\mathbb H,Q^*]1=\nu e.}
\]

The critical operator layer from Sections 11.4--11.7 remains exact but is now a compressed reading,
not extra ontology.  With `R=Lambda^-1`, `A=[R,alpha wedge]`, `B=[R,beta wedge]` and
`V=(alpha wedge)R(alpha wedge)`, one has

\[
\boxed{\mathcal G_c=\mathcal V+2\nu\mathcal B=2\nabla\mathcal A,
\qquad \nabla\mathcal G_c=\nu(\beta\wedge)\mathcal A,
\qquad K'=\frac{\pi^2}{2\nu}(\|\mathcal V\|^2-\|\mathcal G_c\|^2).}
\]

The same `Lambda` generates the positive `Gamma_u`, whose trace is `2K` and whose vorticity
contraction is the Gauss source.  The detailed Gauss floor and null/rank geometry remain in those
sections; they are not repeated in the minimal ontology.

Beneath these operator readings, Section 11.8 gives the lower physical transport statement

\[
\boxed{
K=\pi^{-2}\int Q(r)|r|^{-4}dr,\qquad
\kappa(0)=-2\pi^{-2}\int J(r)\cdot r\,|r|^{-6}dr.}
\]

At every fixed separation, Euler can only exchange common and relative endpoint energy, while heat
only removes it.  The relative critical field `W=v/|r|^2` satisfies

\[
\boxed{K=(2\pi^2)^{-1}\|W\|_2^2,\qquad
M_3=(2\pi^2)^{-1}\|\nabla_cW\|_2^2.}
\]

The same endpoint compatibility gives `Delta_c v=4 Delta_r v`; center relay and separation motion
do not carry independent heat bills.  Incompressibility additionally fixes `K_parallel:K_perp=1:3`,
while the affine defect `C=nabla_c v` obeys a homogeneous pressure-coupled equation.  Since
`K_{|r|>=R}<=8E/(pi R)`, any critical escape must concentrate at the collision diagonal, whose
smooth boundary law is exactly the enstrophy stretching/heat balance.

The material all-scale speed lock already proved in the companion law identifies the same two
critical endpoints of the metric velocity,

\[
\boxed{K=\frac12\|g_t\|_{\dot H_g^{-1/2}}^2,\qquad
M_3=\frac12\|g_t\|_{\dot H_g^{1/2}}^2,}
\]

with `-E'=nu ||g_t||_g^2` at the middle rung.  Scale motion itself is neutral at criticality:
`K_lambda=K`, whereas `M3_lambda=lambda^2 M3` and the parabolic clock scales `lambda^-2`.

The remaining scalar race can now be written without a new scale or owner.  For `M3>0`,

\[
\boxed{
\nu_E:=\frac{\kappa}{M_3},
\qquad
K'=2M_3(\nu_E-\nu).}
\]

The same `nu_E` is the pair least-squares coefficient of `A_u` onto `D omega` and the actual-state
Rayleigh quotient of the lossless critical Euler operator.  Its orthogonal remainder `C_perp` is
exactly the 3D reconfiguration road: `C_perp=0` collapses the velocity image to affine-line/shear
geometry, yet `||C_perp||^2` cancels identically from the instantaneous `K'` square.  Equality
rigidity is therefore not a static coercive gap; it is information for persistence.

The exact closure target remains

\[
\boxed{
\int^T\frac{\kappa(0,t)^2}
{N(t)^2[E(t)Z(t)-K(t)^2]}\,dt<\infty,}
\]

which would contradict the already proved necessary divergence for critical escape.  In the
primitive grammar: **Euler chooses the signed two-way coefficient `nu_E`; heat supplies the fixed
one-way coefficient `nu`; a proof must show that NS cannot return `nu_E` productively with infinite
successful persistence on finite history.**  This theorem is **not proved**, and no global-regularity
claim is made.

Material Hodge turnover remains the natural history gauge through

\[
\partial_tL_g=[\mathcal L_v,L_g],
\qquad
\partial_t\Gamma=D_\Gamma B,
\qquad
\|D_\Gamma B\|_{L^2_g}^2=\|\delta_g\widetilde\beta\|_{L^2_g}^2.
\]

Reynolds, Krein, Poisson, midpoint, Hilbert--Schmidt and Gauss--Bianchi formulations remain exact
compressed readings above this lower pair/current law.  No shell, owner, packet, moving projector or analyst
clock is needed to state the remaining problem.
