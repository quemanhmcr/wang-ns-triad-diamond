# Physical mild-aspect `H1`/swirl local no-escape theorem

This note corrects two tempting but invalid shortcuts in the first H1 no-escape draft.

First, the polarization-only observable

\[
Q_{pol}=\sum_c(\|D_1-D_2\|_F^2+\|D_3\|_F^2)
\]

has an auxiliary relative-coordinate action energy `Q_pol/2`, but the odd-Hermite Young theorem acts on the **three physical roles separately**. For unit role spinors,

\[
\|D u\|^2=\frac12\|D\|_F^2,
\]

so

\[
\sum_{i=1}^3\|F_i^{H1}\|^2
\ge\frac14Q_{pol}
\ge\frac1{100}\|B_{hook}\|^2.
\]

Second, homogeneous polarization transport is `SL(2)`, not unitary. Let each physical role propagator solve `Mdot=G M`. On the existing low-strain packet branch,

\[
\boxed{K:=\int_0^T\|G(t)\|_{op}dt\le1/30.}
\]

Hence both `M` and `M^-1` have singular values in `[e^{-K},e^K]`. If

\[
I_1=\int_0^T\|B_{hook}(t)\|dt,
\]

the interaction-picture forcing has `L1` norm at least `e^{-K}I_1/10`. The interval-variation theorem gives

\[
\boxed{J_1\ge I_1/(11T)}
\]

or an interaction Duhamel impulse containing at least half that forcing mass. Pushing it back to physical roles costs a second factor `e^{-K}`. Arb verifies

\[
e^{-1/30}/10>1/11,
\qquad e^{-2/15}>5/6,
\]

so the **physical three-role** first-Duhamel daughter obeys

\[
\boxed{\delta_1^2\ge I_1^2/480.}
\]

If physical nonlinear sideband feedback is at least half of `delta_1`, it is a high--high interaction branch. Otherwise the surviving total role-daughter energy is at least

\[
\delta^2\ge I_1^2/1920.
\]

At least one of the three physical roles has

\[
\delta_i^2\ge I_1^2/5760.
\]

Below critical sideband size `1/80`, the odd-Hermite one-role theorem and the pair-rescue split give

\[
\boxed{
\operatorname{Def}_{transfer}\ge I_1^2/184320
\quad\text{or}\quad
R_{pair}\ge I_1^2/184320.
}
\]

Thus

\[
\boxed{c_{H1}^{phys}=1/184320.}
\]

Using `I_B<=sqrt(6)I_3+I_1`, the common physical mild-aspect full-curvature constant is

\[
\boxed{c_{mild,curv}^{phys}=1/737280.}
\]

These constants supersede the idealized relative-coordinate/isometric-pullback numbers `1/25600` and `1/102400`. The pointwise Arb bridge `Q_pol>=1/25||B_hook||^2` remains valid and unchanged.

The remaining open H1 issue is source calculus for the interaction forcing variation `J_1`. High-aspect grains are still routed to affine ancestry/reuse, never charged by aspect itself.
