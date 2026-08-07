# Annular pressure collision: multipole decay beats packet packing

This note sharpens the pressure-cancellation branch of the localized SGS ledger.
The exact local-energy identity already showed that pressure is spatial
transport rather than interscale production.  The remaining question was
whether pressure work in a packet moat could repeatedly cancel forward SGS
transfer without paying a fresh/critical-mass cost.

The key point is physical and dimension-specific: **incompressibility removes
the constant pressure mode, upgrading the far pressure kernel from monopole
order `|x-y|^{-3}` to a dipole difference of order `r |x-y|^{-4}`.**  In three
space dimensions shell packing grows only like `R^3`, so one full inverse power
of separation remains summable.

The continuum estimates below are exact elementary consequences of the Newtonian
pressure kernel.  The final fresh-packet statement is presently theorem-level
inside a standard finite packet model (Bernstein, bounded overlap/packing and a
near/far pressure split); constructing that packet model from an arbitrary PDE
block remains conditional.

## 1. Constant pressure cannot do boundary work

For a compactly supported spatial window `chi` and divergence-free resolved
velocity `U`,

\[
\int U\cdot\nabla\chi=0.
\]

The pressure kernel in `R^3` is

\[
K_{ij}(z)=\frac1{4\pi}
\left(3\frac{z_i z_j}{|z|^5}-\frac{\delta_{ij}}{|z|^3}\right).
\]

For every source point `y` one may therefore subtract an arbitrary constant
`K_ij(x_0-y)` inside the pressure boundary work.  Thus far sources appear only
through

\[
K_{ij}(x-y)-K_{ij}(x_0-y).
\]

Direct differentiation gives, componentwise,

\[
|\partial_kK_{ij}(z)|
\le \frac{24}{4\pi}|z|^{-4}.
\]

After summing the tensor components against `u_i u_j`, the clean estimate

\[
\boxed{
\sum_{ij}|\nabla K_{ij}(z)|\,|u_i u_j|
\le 10 |z|^{-4}|u|^2
}
\]

is valid.  Mean value along the packet window then gives the dipole gain

\[
\boxed{
|P_{far}(x)-P_{far}(x_0)|
\lesssim
r\int_{far}\frac{|u(y)|^2}{|y-x_0|^4}\,dy.
}
\]

This is the pressure analogue of using a neutral multipole rather than a raw
Coulomb field.

## 2. Fourth-power decay versus three-dimensional packing

At packet scale `r`, split remote sources into shells with radius comparable to
`2^n r`.  Let `mu_{n,a}` be the scale-critical `L^2` mass of packet `a` in shell
`n`, and suppose a packet frame has the standard geometric packing bound

\[
\#\{a\hbox{ in shell }n\}\le C_{geom}2^{3n}.
\]

The normalized far pressure work is bounded by a dimensionless boundary factor
times

\[
\boxed{
\mathfrak P_{far}
=\sum_{n\ge n_0}2^{-4n}\sum_a\mu_{n,a}.
}
\]

If no packet has critical mass exceeding `mu_*`, then

\[
\mathfrak P_{far}
\le C_{geom}\mu_*
\sum_{n\ge n_0}2^{-n}
= C_{geom}\mu_*2^{1-n_0}.
\]

Hence

\[
\boxed{
W_{far}\le C_{far}\mu_*.
}
\]

The exponent `4-3=1` is the decisive spare power.  Without subtraction of the
constant pressure mode one would have `3-3=0`, a logarithmic/non-summable shell
ledger and no such collision theorem.

## 3. Local pressure work forces critical velocity mass

For the local pressure source, Calderon--Zygmund gives schematically

\[
\|P_{loc}\|_{3/2}\le C_R\|V\|_3^2.
\]

Thus, with `g_chi=r ||grad chi||_infinity`, normalized local pressure work obeys

\[
r|W_{loc}|
\le g_\chi C_R\|V\|_3^3.
\]

If the active local field is frequency localized at `N~r^{-1}`, Bernstein gives

\[
\|V\|_3\le C_B r^{-1/2}\|V\|_2.
\]

Therefore a normalized local pressure cancellation `rho` forces

\[
\boxed{
r^{-1}\|V\|_2^2
\ge
\frac1{C_B^2}
\left(\frac{\rho}{g_\chi C_R}\right)^{2/3}.
}
\]

So the near pressure branch is already a scale-critical local `L^2` mass event,
precisely the currency required by the fresh/reuse energy ledger.

## 4. Packet collision theorem

Abstract the no-fresh estimates into

\[
W_{cancel}
\le C_{near}\mu_*^{3/2}+C_{far}\mu_*.
\]

If the pressure cancellation that must be explained is at least `rho>0`, define

\[
\boxed{
\mu_{fresh}
=
\min\left\{
\left(\frac{\rho}{2C_{near}}\right)^{2/3},
\frac{\rho}{2C_{far}}
\right\}.
}
\]

Then

\[
\boxed{
W_{cancel}\ge\rho
\quad\Longrightarrow\quad
\text{some packet has critical mass }\mu\ge\mu_{fresh},
}
\]

unless one of the packet-frame hypotheses used to obtain the near/far bounds
fails.  Those failures are themselves structural: excessive local mass is the
critical-mass/fresh branch, and excessive shell population is replication to be
charged by the Bellman ledger.

Thus pressure cancellation has been reduced to the same alternatives already
present in the no-escape architecture:

\[
\boxed{
\text{pressure cancellation}
\to
\text{fresh critical mass / replication}
\quad\text{or}\quad
\text{summable dipole tail}.
}
\]

## 5. What remains PDE-level

The unresolved step is no longer an estimate for pressure itself.  One must
construct a transfer-adapted wave-packet/frame decomposition of a genuine SGS
block for which:

1. the active near field has the Bernstein/bounded-overlap estimates used above;
2. remote packet shells satisfy the geometric packing ledger;
3. the pressure source split has a summable cross-interaction remainder; and
4. the critical `L^2` mass produced by the local branch is entered into the
   existing fresh/reuse ancestry graph over a packet lifetime.

The fourth-power pressure multipole decay and its `4-3` summability are exact
continuum facts; only their packet bookkeeping remains conditional.
