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

## 11. The primitive NS object is the actual current operator

Put

\[
\alpha=u^\flat,\qquad \beta=d\alpha,\qquad d\beta=0,
\]

and define the actual graded current operator

\[
\boxed{Q:=\nu\delta+\iota_u,\qquad Q^*=\nu d+\alpha\wedge.}
\]

This is not auxiliary notation: the physical electromotive current and rotational momentum law are

\[
\boxed{e=Q\beta,\qquad \alpha_t+dB=-Q\beta,}
\]

hence, after Leray projection,

\[
\boxed{u_t=-P(Q\beta)^\sharp,\qquad Q_t=\iota_{u_t}=-\iota_{P(Q\beta)^\sharp}.}
\]

The current made by `Q` is therefore the velocity of the state variable which determines `Q`.

### 11.1 Curvature is the failure of nilpotence

Exterior Leibniz and adjointness give

\[
\{d,\alpha\wedge\}=\beta\wedge,\qquad
\boxed{\{\delta,\iota_u\}=(\beta\wedge)^*.}
\]

Since `delta^2=i_u^2=0`, on the whole graded exterior algebra

\[
\boxed{Q^2=\nu(\beta\wedge)^*.}
\]

Thus vorticity curvature is the obstruction to nilpotence of the actual NS current differential.  In particular

\[
\boxed{Q^2\beta=\nu|\beta|^2.}
\]

Cartan's formula gives the vorticity generator

\[
\boxed{\{d,Q\}=\mathcal L_u+\nu L,\qquad \beta_t=-dQ\beta=-\{d,Q\}\beta,}
\]

and associativity gives, because `d beta=0`,

\[
\boxed{[\mathcal L_u+\nu L,Q]\beta=\nu d|\beta|^2.}
\]

Concentration therefore forces noncommutation between the current and the generator made by that current.

### 11.2 Nilpotent chords and finite native chains

The state dependence is affine:

\[
\boxed{Q(u)-Q(v)=\iota_{u-v},\qquad (Q(u)-Q(v))^2=0.}
\]

Every tangent is likewise square-zero and any two tangents anticommute,

\[
\boxed{Q_t^2=0,\qquad \iota_a\iota_b+\iota_b\iota_a=0.}
\]

So non-nilpotence comes only from the interaction of the fixed codifferential with physical contraction, not from nonlinear chords in state space.  This does not by itself control frequency concentration.

Two finite chains expose the same point.  Since `Q^* alpha=nu beta`,

\[
\boxed{
\alpha\xrightarrow{Q^*}\nu\beta\xrightarrow{Q}\nu e
\xrightarrow{Q}\nu^2|\beta|^2\xrightarrow{Q}0.
}
\]

Adjointness of the first two arrows is exactly

\[
\boxed{
\frac12E'=-\langle Q\beta,\alpha\rangle
=-\langle\beta,Q^*\alpha\rangle
=-\nu\|\beta\|_2^2.
}
\]

With the orientation convention of this note,

\[
\boxed{
 dV\xrightarrow{Q}\star\alpha\xrightarrow{Q}\nu\omega^\flat
\xrightarrow{Q}\nu(u\cdot\omega)\xrightarrow{Q}0,\qquad Q^4=0.
}
\]

Velocity, vorticity, helicity density and energy dissipation are readings of one lowering operator.

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

---

## 17. Minimal ontology and the remaining primitive theorem

The deepest state-dependent object reached here is

\[
\boxed{Q=\nu\delta+\iota_u,\qquad Q^*=\nu d+u^\flat\wedge.}
\]

Its core closure laws are

\[
\boxed{Q^2=\nu(\beta\wedge)^*,\qquad
\alpha_t+dB=-Q\beta,\qquad
Q_t=-\iota_{P(Q\beta)^\sharp},}
\]

\[
\boxed{
\mathbb H=Q^*Q+QQ^*,\qquad
\mathbb H^{(1)}-\star^{-1}\mathbb H^{(2)}\star=4\nu S,
}
\]

and, for closed `eta=star b^flat`,

\[
\boxed{
\mathbb H^{(1)}(Q\eta)-Q(\mathbb H^{(2)}\eta)
=\nu^2d(\omega\cdot b)+\nu[\,b\times(u\times\omega)\,]^\flat.
}
\]

So one endogenous operation closes the loop

\[
\boxed{
\text{operation}\to\text{curvature}\to\text{degree imbalance}
\to\text{self-force}\to\text{operation}.
}
\]

The midpoint

\[
M=\nu\delta+\frac12\iota_u
\]

is the canonical critical/Poynting reading of `Q`, not a second fundamental operator, and

\[
\boxed{M(M\beta)=\frac\nu2|\beta|^2.}
\]

Sections 4--10 show that the Krein/sech/Poisson/Reynolds/HS language is a higher representation of this same strain/current law.  Critical escape still requires

\[
\boxed{
\int^T\frac{\kappa(0,t)^2}{N(t)^2[E(t)Z(t)-K(t)^2]}
\,dt=\infty.
}
\]

The primitive persistence question is lower than a moving eigenspace: can the **projected** midpoint current `P(M beta)` remain near-null in the shifted critical metric while the unprojected current satisfies

\[
M(M\beta)=\frac\nu2|\beta|^2
\]

and the actual current continuously changes the differential through `Q_t`?

The exact candidate turning operator is already present:

\[
\boxed{[\mathbb H,Q]=Q^*Q^2-Q^2Q^*.}
\]

Its degree-two-to-degree-one non-exact component is the curvature-transverse Lamb interaction, from which the Lamb field is locally reconstructed when `omega != 0`.  The next theorem should therefore be sought in the form

\[
\boxed{
\text{persistent critical near-kernel of }P M
\Longrightarrow
\text{non-summable objective graded turning generated by }[\mathbb H,Q].
}
\]

That implication is **not proved**.  To close regularity it must be quantitatively coupled to

\[
\partial_tL_g=[\mathcal L_v,L_g],\qquad
\partial_t\Gamma=D_\Gamma B,\qquad
\|D_\Gamma B\|_{L^2_g}^2=\|\delta_g\widetilde\beta\|_{L^2_g}^2
\]

and the transverse heat-memory law, so that the turnover forced by critical persistence cannot be sustained on the relevant finite positive-energy history.  Establishing that bridge is the remaining theorem, not an already available action bound.

The methodological rule is now:

\[
\boxed{
\text{do not add an object unless it is generated by }Q,Q^*,Q^2
\text{ or by their natural material pullback.}
}
\]

No shell, owner, packet, moving projector, entropy budget or analyst clock is needed to state the remaining problem.
