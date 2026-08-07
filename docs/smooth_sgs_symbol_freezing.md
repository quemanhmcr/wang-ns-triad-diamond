# Smooth SGS symbol freezing: finite transfer edges with summable error

The current PDE frontier asks for a transfer-preserving wave-packet or frequency
atomic realization of one smooth SGS crossing block.  The multiplier part of
that problem is easier than a general Coifman--Meyer decomposition because the
SGS transfer has an explicit Fourier trilinear form and the common spectral moat
keeps all active frequencies away from zero and from filter singularities.

This note isolates the exact elementary reduction.  Spatial packet tracking is
still a separate PDE problem.

## 1. Global smooth SGS work is a trilinear Fourier multiplier

For a smooth resolved filter `G_N`, the spatial average of

\[
-\nabla\bar u:\tau_N(u,u)
\]

is a Fourier trilinear form

\[
T_{M_N}(f,g,h)
=\iint M_N(\xi,\eta)
 f(\xi)g(\eta)\overline{h(\xi+\eta)}\,d\xi d\eta.
\]

After factoring the dimensional power of `N`, the normalized symbol is
scale-covariant: it is a smooth function of `xi/N` and `eta/N`.  The Leray
symbol and helical projectors are smooth on the crossing shells because those
shells are bounded away from zero.  Hence on every certified compact crossing
block there is a finite dimensionless Lipschitz constant `L_M`.

No statistical turbulence closure is used here; this is simply the Fourier
representation of the resolved nonlinear work.

## 2. Freeze the symbol on relative frequency cells

Partition the dimensionless crossing domain into cells of Euclidean diameter at
most `h`.  On each cell choose a representative `z_Q` and define the frozen
symbol

\[
M_h(z)=M(z_Q),\qquad z\in Q.
\]

If `M` is `L_M`-Lipschitz, then pointwise

\[
\|M-M_h\|_\infty\le L_Mh.
\]

For the critical Young exponents `3/2,3/2,3/2`, the scalar sharp constant is

\[
A_3=\left(\frac{\sqrt3}{2}\right)^3.
\]

Therefore

\[
\boxed{
|T_M(f,g,h)-T_{M_h}(f,g,h)|
\le
A_3 L_M h
\|f\|_{3/2}\|g\|_{3/2}\|h\|_{3/2}.
}
\]

This is not a heuristic symbol expansion.  It is the pointwise multiplier error
followed by sharp Young.

The frozen form is a finite/countable sum of cell interactions.  Those cells
are precisely the atomic frequency edges needed by the transfer/Hodge ledger;
phase polarization can then be applied to their signed contributions.

## 3. The multiscale synthesis error is summable

Choose relative cell diameters

\[
h_j=(j+3)^{-2}.
\]

The integral test gives

\[
\sum_{j\ge0}h_j
\le \frac19+\frac13=\frac49.
\]

More importantly, if each generation is allowed its own prefactor `epsilon_j`
one can simply choose `h_j` so that

\[
A_3L_{M,j}h_j\le\epsilon_j,
\qquad \sum_j\epsilon_j<\infty.
\]

Thus smooth multiplier freezing contributes directly to the existing summable
cross-error ledger `Xi`; there is no packet-count factor.

## 4. Why the common moat matters physically

Without scale separation the filter multiplier can change rapidly across a cell
that straddles the generation boundary, and the helical chart can approach the
zero-frequency singularity.  The new crossing-to-moat theorem removes exactly
that problem: a one-quarter physical-transfer/Hodge subcore lies in parent and
child shells of halfwidth `<2/25` with a smooth moat `>1/25`.

On this compact set the normalized SGS/Leray/helical symbol has a uniform finite
`L_M`.  The finite edge model is therefore a controlled approximation of the
actual smooth physical flux, not a free-standing graph ansatz.

## 5. What remains

This closes only the **frequency multiplier freezing** part of packetization.
The remaining hard construction is spatial and temporal:

1. localize the one-shot Christ Gaussian profile by moving physical windows;
2. control the overlap between different spatial grains and the `RU`, viscous
   and window-transport leakage;
3. preserve the frozen frequency-edge identities over a packet lifetime; and
4. synchronize the resulting grains across generations.

The multiplier/synthesis error itself can now be made summable by a direct
relative-cell schedule.
