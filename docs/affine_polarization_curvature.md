# Spatial strain/polarization curvature on an affine grain

The scalar Gaussian forcing theorem sees only the fully symmetric part of the
grain-normalized velocity Hessian.  The complementary vector channel is the
spatial variation of the **physical symmetric strain** across the grain.

Let

\[
A(x)=A_0+H[x-X]+\cdots,\qquad x-X=Lz,
\]

and define

\[
\boxed{
C_{ijc}=\frac12(H_{ijk}+H_{jik})L_{kc}.
}
\]

Then

\[
S(z)-S(0)=C_c z_c.
\]

Differentiated incompressibility gives `tr C_c=0`.  For the transfer-relevant
extremal observable

\[
Q_{rel}(S)=\|D_\Pi\|_F^2+\|D_1-D_2\|_F^2+\|D_3\|_F^2,
\]

the Arb-certified theorem already gives

\[
Q_{rel}(S)\ge\frac12\|S\|_F^2.
\]

Since the normalized grain coordinate has `E[z_c z_d]=delta_cd`, quadratic
homogeneity yields the exact integrated identity and bound

\[
\boxed{
\mathbb E_z Q_{rel}(S(z)-S(0))
=\sum_c Q_{rel}(C_c)
\ge\frac12\|C\|_F^2.
}
\]

This is the vector/helical complement of the third-Hermite scalar forcing

\[
\|F_\perp\|_2^2/\|\psi\|_2^2=(3/8)\|\operatorname{Sym}B\|_F^2.
\]

In particular the quadratic swirl kernel `Sym B=0` is not declared free: it is
routed to spatial variation of the physical strain/polarization generator.  The
strength of that channel relative to the affine-normalized `||B||` norm may
depend on ellipsoid aspect; this note does not insert a false aspect-independent
coercivity constant.
