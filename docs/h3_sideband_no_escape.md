# H3 curvature sideband: a local no-escape theorem

This note combines the affine curvature connection, first-Duhamel coherence, the critical Young footprint, and the odd-Hermite triad selection rule into one master-facing statement.

Let

\[
I_3=\int_0^T\|\operatorname{Sym}\widetilde B(t)\|_Fdt,
\qquad
J_3=\int_0^T\|\operatorname{Sym}\widetilde S(t)\|_Fdt,
\]

in the curvature interaction frame, so `dot Btilde=Stilde`.  The normalized `H3` forcing has `L2` norm

\[
\|f_3(t)\|_2/\|g\|_2=\sqrt{3/8}\|\operatorname{Sym}\widetilde B(t)\|_F.
\]

The interval-variation theorem gives

\[
\boxed{J_3\ge I_3/T}
\]

or a first-Duhamel impulse

\[
\boxed{\delta_1^2\ge\frac3{32}I_3^2.}
\]

If subsequent nonlinear sideband feedback has norm at least `delta1/2`, it is already a genuine high--high/sideband interaction event.  Otherwise the surviving odd daughter satisfies

\[
\boxed{\delta^2\ge\frac3{128}I_3^2.}
\]

Let `sigma^2` denote its second moment in the critical `|G|^(3/2)` Gaussian measure.  For an `H3` Gaussian, the variance change from the `L2` measure to the critical Young measure gives `sigma>=delta`.

If

\[
\sigma\ge1/80,
\]

there is a definite daughter-capacity event.  Otherwise the odd-sideband norm theorem gives the single-role normalized transfer loss

\[
d_0\ge\sigma^2/16\ge\frac3{2048}I_3^2.
\]

The only way another odd sideband can offset this is through the quadratic pair-sideband terms.  Splitting `d0` in half gives the clean local no-escape constant

\[
\boxed{c_{H3}=\frac3{4096}.}
\]

Hence every nontrivial `H3` curvature impulse enters at least one branch:

\[
\boxed{
\begin{array}{ll}
J_3\ge I_3/T, & \text{acceleration-Hessian dephasing source},\\[1mm]
\|R_{Duhamel}\|\ge\delta_1/2, & \text{nonlinear sideband feedback},\\[1mm]
\sigma\ge1/80, & \text{definite daughter capacity},\\[1mm]
\operatorname{Def}_{transfer}\ge(3/4096)I_3^2, & \text{base-edge cost},\\[1mm]
R_{pair}\ge(3/4096)I_3^2, & \text{pair-sideband interaction}.
\end{array}}
\]

The last branch is deliberately **not** declared Bellman cost yet: it must be inserted into the spacetime component/ancestry graph without assuming that overlapping Hermite modes are disjoint spatial packets.  That is now a sharply identified interface problem rather than an unidentified forcing term.
