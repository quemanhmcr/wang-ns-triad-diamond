# Certified single-edge helical stability

This module is a finite-dimensional theorem for the scale-normalized helical
transfer edge.  It is **not** the missing PDE flux bridge and is not a proof of
Navier--Stokes regularity.

## 1. Exact reduction over helicity signs

Normalize the child magnitude to one and order the parent magnitudes

\[
0<x\le y<1,\qquad x+y>1.
\]

For helicities \(s_x,s_y,s_z\in\{\pm1\}\), the only sign-dependent factor in
\(\mathcal J\) is

\[
P=|s_xx-s_yy|\,|s_xx+s_yy+s_z|.
\]

There are only two cases up to simultaneous sign reversal.

* If \(s_x=-s_y\), then \(P=(x+y)(1\pm(x-y))\), so the larger value is
  \((x+y)(1+y-x)\).
* If \(s_x=s_y\), the larger possibility is
  \((y-x)(1+x+y)\), and
  \[
  (x+y)(1+y-x)-(y-x)(1+x+y)=2x>0.
  \]

Therefore the full sign problem reduces **exactly** to

\[
\boxed{P\le (x+y)(1+y-x).}
\]

The maximizing orbit has opposite-helicity parents.  At \(x=y\), the two child
helicities have equal magnitude, so there is no fictitious discrete child-sign
gap.

Put

\[
s=x+y,\qquad d=y-x,\qquad L=\log(1/y).
\]

Using Heron's formula, the sign-reduced envelope is

\[
\boxed{
J_{\rm env}^2(x,y)=
\frac{L^2s^2(s^2-1)(1+d)^3(1-d)}
{8(s^2-d^2)^2}.
}
\]

This identity is exact.

## 2. The symmetric critical point is exact and unique

On \(x=y=r\),

\[
J_{\rm sym}(r)=
\frac{\sqrt{4r^2-1}}{4\sqrt2\,r}\log(1/r).
\]

Its critical equation is

\[
-\log r=4r^2-1.
\]

The function \(f(r)=-\log r-4r^2+1\) has

\[
f'(r)=-1/r-8r<0,
\]

so it has at most one root on \((1/2,1)\).  The Arb certificate verifies a
sign change on the rational bracket

\[
0.6109041018306<r_*<0.6109041018307,
\]

hence this root exists and is unique.  Define

\[
\gamma_*=-\log r_*,\qquad J_*=J_{\rm sym}(r_*).
\]

## 3. Coordinates forced by scale holonomy

For an ordered edge define

\[
u=\log(y/x)\ge0,
\qquad
v=-\frac12(\log x+\log y)-\gamma_*.
\]

Thus

\[
x=r_*e^{-v-u/2},\qquad y=r_*e^{-v+u/2}.
\]

The variable \(u\) is the parent-scale imbalance and \(v\) is exactly the
mean parent-to-child scale residual used by the Hodge module.  If
\(r_x,r_y\) are the two preferred arc residuals, then the already-proved exact
identity is

\[
\boxed{r_x^2+r_y^2=\frac{u^2}{2}+2v^2.}
\]

These are therefore not arbitrary optimization coordinates: they are the
coordinates seen by the physical scale-transfer graph.

## 4. Local analytic certificate

Let

\[
U_0=V_0=\frac{2}{25},\qquad A=\frac1{50},\qquad B=1.
\]

On the whole rectangle

\[
0\le u\le U_0,\qquad |v|\le V_0,
\]

write \(R=r_*e^{-v}\),
\(S=2R\cosh(u/2)\), \(\delta=2R\sinh(u/2)\), and
\(L=\gamma_*+v-u/2\).  Direct differentiation of the exact envelope gives

\[
\partial_u\log J_{
m env}
=-\frac1{2L}
+\frac12\tanh(u/2)
+\frac{S\delta}{2(S^2-1)}
+\frac{S(1-2\delta)}{2(1-\delta^2)}.
\]

Arb ball arithmetic certifies on the entire rectangle that

\[
\partial_u\left(1-\frac{J_{\rm env}}{J_*}\right)>\frac1{50}.
\]

Along \(u=0\), put \(a=R^{-2}\) and \(q=\sqrt{4-a}\).  Exact
differentiation gives

\[
\partial_v^2\left(1-\frac{J_{\rm sym}}{J_*}\right)
=
\frac{2a/q+(\gamma_*+v)(2a/q+a^2/q^3)}{4\sqrt2\,J_*}.
\]

The root equation gives \(D(0)=D'(0)=0\), and Arb certifies the second
derivative is strictly larger than \(2\) for \(|v|\le2/25\).  Integrating the
two derivative bounds yields the theorem-level mixed stability inequality

\[
\boxed{
1-\frac{J_{\rm env}}{J_*}
\ge \frac1{50}u+v^2.
}
\]

Because \(0\le u\le2/25\),

\[
\frac1{50}u\ge\frac14u^2.
\]

Consequently

\[
\boxed{
1-\frac{J_{\rm env}}{J_*}
\ge
\frac12\left(\frac{u^2}{2}+2v^2\right)
=\frac12(r_x^2+r_y^2).
}
\]

This is the desired single-edge-to-Hodge coefficient on the near-extremal
network:

\[
\boxed{c_{\rm stab}=1/2.}
\]

## 5. Global exclusion

The complement of the local rectangle is certified separately.

For \(y\ge0.9\), the exact geometric form

\[
J=\frac{\log(1/y)}{4\sqrt2}\,\sin\theta\,(x+y)(1+y-x)
\]

and the triangle inequalities give

\[
J\le\frac{\log(1/y)}{\sqrt2}
\le\frac{\log(10/9)}{\sqrt2}<0.99J_*.
\]

For \(1/2\le y\le0.9\), parameterize the whole ordered triangle domain by

\[
a=2y-1,\qquad
s=1+\lambda a,\qquad
d=(1-\lambda)a,\qquad0\le\lambda\le1.
\]

Then \(x=(s-d)/2\), and the exact formula for \(J_{\rm env}^2\) is evaluated by
Arb on a dyadic branch-and-bound partition.  A box is accepted only if either
it lies wholly in the local rectangle or Arb proves

\[
J_{\rm env}^2<(0.99J_*)^2.
\]

Therefore

\[
\boxed{
(u>2/25)\ \text{or}\ (|v|>2/25)
\quad\Longrightarrow\quad
1-J/J_*\ge1/100.
}
\]

The proof uses inclusion-preserving Arb ball arithmetic; random tests are only
adversarial regression checks and are not part of the certificate.

## 6. How this enters the no-escape architecture

For every retained near-extremal interaction edge,

\[
\operatorname{Def}_e\ge\frac12\mathcal D_e,
\qquad
\mathcal D_e=r_{e,1}^2+r_{e,2}^2.
\]

Thus for transfer/Hodge conductances \(w_e\),

\[
\sum_e w_e\operatorname{Def}_e
\ge\frac12\sum_e w_e\mathcal D_e,
\]

so the local multiplier deficit pays at least one half of the Hodge residual
energy on the good-edge network.  Edges outside the local rectangle each have
a certified \(1/100\) deficit; if a block has small average deficit, their
total transfer weight is therefore automatically small and they can be removed
before Hodge synchronization at a quantitatively controlled cost.

This closes the **finite-dimensional single-edge stability gap** that was
previously only numerical.  It does not close the PDE bridge: one still has to
show that a genuine Navier--Stokes near-extremal scale-flux block produces the
weighted atomic edge ledger to which this certificate applies.
