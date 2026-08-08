# Physical transfer-weighted defect moat: summable cross-cell `Xi` without packet counting

The coherent localization operators reconstruct the velocity and the trilinear form exactly, so the remaining global interface is not a synthesis error. It is the **actual physical transfer carried by cell triples that cross the selected ancestry components**.

This can be controlled directly by the transfer-weighted single-edge defect. No Gaussian coefficient tail or coordinate-grid boundary is needed.

## 1. Physical transfer controls the defect moment

Let \(w_e\) be normalized unsigned capacity weights and

\[
\mathcal D_e=r_{e,1}^2+r_{e,2}^2
\]

be the scale/Hodge edge defect. On the local certified box,

\[
\operatorname{Def}_e\ge\frac12\mathcal D_e.
\]

For block deficit \(\epsilon=1-R\),

\[
\sum_ew_e\mathcal D_e\le2\epsilon.
\]

Fix the certified good core \(G_{10^{-4}}\). If \(\epsilon<1/20000\), its capacity mass is at least \(1/2\). Therefore the normalized capacity law on the core has

\[
\mathbb E_{cap,G}\mathcal D\le4\epsilon.
\]

The physical positive child-transfer law \(\widetilde w\), normalized on the same core, obeys

\[
\frac{50}{53}w_G\le\widetilde w\le\frac{53}{50}w_G.
\]

Hence

\[
\boxed{
\mathbb E_{phys,G}\mathcal D
\le\frac{106}{25}\epsilon.
}
\]

This is the key input: the moat is controlled by the **actual positive Navier--Stokes child-energy transfer law**.

## 2. Put the moat in transfer defect, not in a coordinate grid

At one coherent-cell depth choose \(R>0\). Divide

\[
[R/2,R]
\]

into \(M\) equal bins. One bin has physical transfer mass at most \(1/M\); remove it. Let its endpoints be \(\ell<u\).

Connect coherent/packet vertices using every active triad edge with

\[
\mathcal D_e<\ell.
\]

Consider an edge joining two different resulting components. It cannot have defect below \(\ell\), because then it would have joined its three vertices. Therefore every cross-component edge is either

1. in the selected moat \([\ell,u)\); or
2. in the tail \(\mathcal D_e\ge u\ge R/2\).

The moat costs at most \(1/M\). Markov gives

\[
\widetilde w\{\mathcal D\ge R/2\}
\le
\frac{2\mathbb E_{phys}\mathcal D}{R}.
\]

Thus

\[
\boxed{
\eta_{cross}
\le
\frac1M
+
\frac{2\overline{\mathcal D}}R.
}
\]

No packet count appears.

## 3. A summable recursive schedule

Choose

\[
M_j=M_0(j+2)^2,
\qquad
R_j=R_0(j+2)^2.
\]

Since

\[
\sum_{j\ge0}\frac1{(j+2)^2}
=\frac{\pi^2}{6}-1
<\frac{13}{20},
\]

one has

\[
\boxed{
\sum_j\eta_{cross,j}
\le
\frac{13}{20}
\left(
\frac1{M_0}
+
\frac{2\overline{\mathcal D}}{R_0}
\right).
}
\]

On the near-extremal branch, \(\overline{\mathcal D}\le(106/25)\epsilon_*\). Both \(M_0\) and \(R_0\) are selection parameters, so the total omitted physical-transfer mass is finite and can be made arbitrarily small.

For the logarithmic master penalty

\[
\xi_j=\log(1+\eta_je^{C_j}),
\]

if the recursive extraction is being performed on the low-cost branch \(C_j\le C_{cap}\), then

\[
\boxed{
\sum_j\xi_j
\le
e^{C_{cap}}
\sum_j\eta_j<\infty.
}
\]

Costly blocks are already stopped in the multiplicative branch and need not be refined by this moat.

## 4. Percolation is not an interface escape

A large \(R_j\) may make the low-defect graph very connected. This is not a failure of the construction. After deleting the summable cross edges, every retained connected 3-uniform interaction component obeys

\[
\boxed{(n-1)+\beta=2m.}
\]

Therefore it is fresh-rich or cycle-rich. The moat theorem and the fresh/cycle theorem are complementary:

\[
\boxed{
\text{summably decoupled components}
\quad\lor\quad
\text{connected sticky fresh/reuse ancestry}.
}
\]

The construction never assumes spatial disjointness of Gaussian packets.

## 5. Why this is the canonical global `Xi`

Previous nested-Gaussian notes used the explicit envelope `exp(-D^2)` to obtain a tail. The present theorem is stronger for the actual PDE architecture: it uses only

- the exact physical positive transfer law;
- the certified single-edge stability defect;
- Markov and a transfer-weighted moat;
- exact incidence topology after the cut.

It is therefore compatible with the continuous coherent localization operators and does not require a separate discrete Gaussian synthesis.

What remains after this theorem is no longer an interface-summability mechanism. It is the branch-by-branch PDE audit ensuring every efficient block reaches the same signed-good physical transfer law, or has already paid a source/sideband/high-strain/fresh/reuse cost.
