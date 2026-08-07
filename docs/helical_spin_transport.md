# Helical spin transport and the triad-normal gauge

The helical eigenspace at each nonzero Fourier wavevector is a complex line.  A
single global smooth choice of phase is impossible: in the standard spherical
gauge

\[
h_s=(e_\theta+i s e_\phi)/\sqrt2,
\qquad i k\times h_s=s|k|h_s,
\]

the Berry connection and curvature are

\[
\mathcal A_s=i h_s^*dh_s=s\cos\theta\,d\phi,
\qquad
\mathcal F_s=d\mathcal A_s=-s\sin\theta\,d\theta\wedge d\phi,
\]

so

\[
\frac1{2\pi}\int_{S^2}\mathcal F_s=-2s.
\]

Thus the helicity line has the spin-one Chern charge.  This topology does **not**
mean that a single rotating Navier--Stokes triad pays a Berry-phase deficit.
That would incorrectly charge an exact rotational symmetry.

For a nondegenerate triad choose its oriented unit normal `n`.  Every carrier in
the plane has the SO(3)-covariant gauge

\[
\boxed{
h_s(k;n)=\frac{n\times\widehat k+i s n}{\sqrt2}.
}
\]

If the whole triad is rotated by `R`, then exactly

\[
h_s(Rk;Rn)=R h_s(k;n).
\]

In this gauge the Waleffe triple product is purely imaginary.  Therefore, away
from a coupling zero, its phase is locally constant (`+pi/2` or `-pi/2`).  A
single-triad geometric phase can be absorbed into this moving physical frame.

The genuinely observable geometric phase appears when the **same Fourier mode**
is reused by two triads whose planes have different normals.  If `n2` is the
rotation of `n1` by signed dihedral angle `psi` about `+k`, then

\[
\boxed{
h_s(k;n_2)=e^{-is\psi}h_s(k;n_1).
}
\]

This is the spin-one transition function.  Products of these transition phases
around a reuse cycle are gauge invariant and cannot be removed by choosing a
phase convention for individual modes.

There is a second physical consequence.  Let the real symmetric trace-free
strain seen in a transverse frame be

\[
D=\begin{pmatrix}\delta&\beta\\\beta&-\delta\end{pmatrix}.
\]

In the circular/helical basis `(h_+,h_-)`,

\[
\boxed{
D_{hel}=
\begin{pmatrix}
0&\delta-i\beta\\
\delta+i\beta&0
\end{pmatrix}.
}
\]

Thus nonconformal transverse strain is exactly an opposite-helicity conversion
operator, with conversion rate

\[
|\delta-i\beta|=\|D\|_F/\sqrt2.
\]

For the frozen incompressible strain `D=diag(d,-d)` with carrier normal fixed,
a pure positive-helicity Kelvin amplitude evolves as

\[
u_+(t)=\cosh(dt),\qquad u_-(t)=-\sinh(dt).
\]

So the same deformation that drives scalar triad-shape rigidity also creates
opposite-helicity content.  In a genuine packet theorem nonlinear forcing and
spatial variation must still be controlled, but helicity conversion is not an
independent mysterious error: it is another direct readout of nonconformal
strain.
