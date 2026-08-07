# Resistance-to-Bellman stopping without a log loss

This is a finite weighted-tree theorem used as a stopping module inside the grain-cascade program; it is not a Navier--Stokes regularity theorem.

Let a tree have normalized conductances `c_e>0`, `sum c_e=1`, terminal probability `p`, edge resistances `r_e=1/c_e`, and tree resistance distance `R(i,j)`.

For a parameter `Lambda>0`, cut edges independently with

\[
q_e=1-e^{-r_e/\Lambda}.
\]

For a fixed pair `(i,j)`, the probability that no edge of its unique path is cut is exactly

\[
\prod_{e\in P(i,j)}e^{-r_e/\Lambda}=e^{-R(i,j)/\Lambda}.
\]

Hence if `Q(P)=sum_C p(C)^2` is the collision probability of the random partition,

\[
\mathbb E Q=\mathbb E_{i,j}e^{-R(i,j)/\Lambda}=1-F_\Lambda,
\qquad
F_\Lambda:=\mathbb E_{i,j}(1-e^{-R(i,j)/\Lambda}).
\]

The expected deleted conductance is

\[
\bar C_\Lambda=\sum_e c_e(1-e^{-1/(c_e\Lambda)})
\le \rho_\Lambda,
\qquad
\rho_\Lambda:=\sum_e\min(c_e,\Lambda^{-1}).
\]

A one-line probabilistic stopping argument gives one deterministic cut with

\[
Q\le1-F_\Lambda/2,
\]

and

\[
C_{\rm cross}\le \frac{2(2-F_\Lambda)}{F_\Lambda}\,\bar C_\Lambda
\le \frac{2(2-F_\Lambda)}{F_\Lambda}\rho_\Lambda.
\]

Thus a resistance-separated pair distribution either yields a Bellman partition with positive collision-entropy gain and small cross error (when `rho_Lambda` is small), or `rho_Lambda` itself is large.

The latter case is not an escape. For the conductance distribution itself,

\[
Q_e:=\sum_e c_e^2\le \max_e c_e.
\]

For an edge attaining the maximum, all other edges contribute at least `rho_Lambda-1/Lambda` to their total mass. Therefore

\[
\max_e c_e\le1-\rho_\Lambda+\Lambda^{-1},
\]

and hence

\[
\boxed{H_2(c)=-\log\sum_e c_e^2\ge -\log(1-\rho_\Lambda+\Lambda^{-1}).}
\]

So the soft-complexity term created by resistance routing becomes atomic collision entropy, which is handled by the preceding atomic-to-component chain rule.

If `Lambda` is a median of the pair resistance distribution, then

\[
F_\Lambda\ge \tfrac12(1-e^{-1}),
\]

so both the Bellman entropy gain and the simultaneous-cut constant are universal. This gives a no-log stopping dichotomy:

\[
\text{high median resistance}
\Rightarrow
\begin{cases}
\text{cheap Bellman cut},\\
\text{or conductance atomic entropy}.
\end{cases}
\]

The low-median branch is sent to the multicommodity Hodge/synchronization module.
