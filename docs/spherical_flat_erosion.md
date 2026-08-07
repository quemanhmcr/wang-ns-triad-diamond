# Spherical flat-network erosion

This is a finite-dimensional geometric module. It is not a proof of Navier--Stokes regularity.

## 1. Exact midpoint barrier

Let `p,q,n` be unit vectors with `n.p,n.q>0`, let

\[
\theta=d_{S^2}(p,q)<\pi,
\qquad
m=\frac{p+q}{|p+q|},
\]

and define the open-hemisphere barrier

\[
\Phi_n(x)=-\log(n\cdot x).
\]

Since `|p+q|=2 cos(theta/2)`, the arithmetic--geometric mean inequality gives

\[
n\cdot m
=\frac{n\cdot p+n\cdot q}{2\cos(\theta/2)}
\ge
\frac{\sqrt{(n\cdot p)(n\cdot q)}}{\cos(\theta/2)}.
\]

Therefore

\[
\boxed{
\Phi_n(m)
\le
\frac{\Phi_n(p)+\Phi_n(q)}2
-\kappa(\theta),
\qquad
\kappa(\theta)=-\log\cos(\theta/2).
}
\]

For the optimal helical triad angle

\[
\theta_*\approx70.1383^\circ,
\qquad
\kappa_*\approx0.2003318956.
\]

This is the spherical replacement for the one-dimensional angular-diameter erosion ledger.

## 2. Distinguished-lineage ledger

Let `x_{j+1}` be the midpoint of the distinguished parent `x_j` and a companion `q_j`, with every parent angle equal to `theta_*`. Then

\[
\Phi_n(q_j)
\ge
2\Phi_n(x_{j+1})-\Phi_n(x_j)+2\kappa_*.
\]

Summing gives

\[
\boxed{
\sum_{j=0}^{L-1}\Phi_n(q_j)
\ge
2L\kappa_*-\Phi_n(x_0)
+
\sum_{j=1}^{L-1}\Phi_n(x_j)
+2\Phi_n(x_L).
}
\]

In particular,

\[
\sum_{j<L}\Phi_n(q_j)
\ge2L\kappa_*-\Phi_n(x_0).
\]

Thus a flat lineage contained in one open hemisphere must import companion barrier at a rate linear in depth.

## 3. Cap contraction and finite lifetime

If every parent lies in the cap `n.x >= mu>0`, then every exact child satisfies

\[
n\cdot m\ge\frac{\mu}{\cos(\theta_*/2)}.
\]

After `L` no-fresh generations,

\[
\mu_L\ge\mu_0\cos(\theta_*/2)^{-L}.
\]

Since `mu_L<=1`, a nonempty no-fresh cascade must satisfy

\[
\boxed{
L\le
\frac{\log(1/\mu_0)}{-\log\cos(\theta_*/2)}.
}
\]

Any infinite flat cascade must therefore leave every open hemisphere infinitely often.

## 4. The balanced exception

For a finite direction set `S`, let `K=conv(S)`. If `0` is outside `K`, the nearest-point projection gives an open hemisphere containing `S`. Hence failure of the hemispherical branch implies

\[
0\in\operatorname{conv}(S).
\]

By Caratheodory in `R^3`, at most four directions already certify this. Moreover any such certificate contains two directions separated by at least

\[
\boxed{
\theta_{\rm tet}=\arccos(-1/3)\approx109.4712^\circ.
}
\]

Indeed, if `0=sum lambda_i x_i` with at most four positive weights and all pairwise dot products were larger than `-1/3`, then

\[
0=\left|\sum_i\lambda_i x_i\right|^2
>
\frac43\sum_i\lambda_i^2-\frac13
\ge0,
\]

a contradiction. Equality is achieved by a regular tetrahedron.

Thus the only escape from hemispherical erosion is a quantitatively broad, balanced configuration. Such a configuration cannot be one narrow triad grain; it must be routed back into the branching/percolation modules.

## 5. Barycenter amplification

Let a coupling of parent directions have equal marginals `mu` and be supported on pairs at angle `theta_*`. Push the coupling forward by spherical midpoint to obtain the child measure `nu`. Their barycenters satisfy exactly

\[
\boxed{
b(\nu)=\frac{b(\mu)}{\cos(\theta_*/2)}.
}
\]

Consequently a depth-`L` no-fresh balanced chain obeys

\[
|b(\mu_0)|\le\cos(\theta_*/2)^L.
\]

A stationary flat cell must have zero barycenter. Long-lived flat networks are therefore forced toward the balanced exceptional class, not merely toward arbitrary spherical spread.

## 6. Current trichotomy

The finite-dimensional architecture is now

\[
\boxed{
\begin{cases}
\text{hemispherical flat component}
&\Rightarrow \text{barrier erosion / fresh companion cost},\\
\text{non-hemispherical component}
&\Rightarrow \text{balanced simplex of diameter at least }\theta_{\rm tet},\\
\text{curved component}
&\Rightarrow \text{Hodge cycle energy}.
\end{cases}}
\]

The remaining PDE bridge is to convert the barrier or balanced-simplex certificates into a lower bound for fresh critical mass in the transfer-weighted Gaussian grain decomposition.

## 7. Balanced flat chains pay collision entropy

The barycenter rigidity has a direct entropy consequence. Let

\[
\mu=\sum_i w_i\delta_{x_i},
\qquad
b=\sum_iw_ix_i,
\qquad
Q=\sum_iw_i^2.
\]

If `w_*` is the largest atomic weight, then

\[
|b|\ge w_*-(1-w_*)=2w_*-1.
\]

Hence

\[
w_*\le\frac{1+|b|}{2},
\qquad
Q\le w_*\le\frac{1+|b|}{2}.
\]

Therefore the collision entropy satisfies the exact bound

\[
\boxed{
H_2(\mu)=-\log Q
\ge
\log\frac{2}{1+|b|}.
}
\]

Now consider a depth-`L` flat chain built from equal-marginal couplings. Since

\[
b_{j+1}=b_j/\cos(\theta_*/2),
\]

boundedness of the terminal barycenter gives

\[
|b_j|\le\cos(\theta_*/2)^{L-j}.
\]

Thus

\[
\boxed{
H_2(\mu_j)
\ge
\log\frac{2}{1+\cos(\theta_*/2)^{L-j}}.
}
\]

Summing over levels,

\[
\sum_{j=0}^{L-1}H_2(\mu_j)
\ge
L\log2
-
\sum_{r=1}^{L}\log(1+c_*^r),
\qquad c_*=\cos(\theta_*/2).
\]

The correction is uniformly bounded because

\[
\sum_{r\ge1}\log(1+c_*^r)<\infty.
\]

Numerically the limiting correction is about `3.76730`. Consequently the balanced exception does not provide a zero-cost infinite branch: under the equal-marginal hypothesis it pays asymptotically `log 2` of atomic collision entropy per generation.

This is the spherical counterpart of the transfer-weighted Bellman cost. The remaining task is to prove that the PDE grain coupling is sufficiently balanced, or else quantify the imbalance itself as fresh mass.
