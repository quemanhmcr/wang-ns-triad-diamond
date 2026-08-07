# Ancestor reservoir synchronization: a spectral half-life from Kelvin transport

The Onsager increment collision can return a dominant low/base-frequency packet rather than immediate high-frequency dissipation.  That packet is an **ancestor reservoir candidate**, not automatically a selected child.  This note proves that one materially reused reservoir cannot service infinitely many successive near-extremal generations for free.

The mechanism is physical: velocity increments are Galilean/sweeping neutral, while a Kelvin covector can only move through frequency space at the rate supplied by symmetric strain.

## 1. Near-extremal generations outrun a low-strain reservoir

On the signed-good core with threshold `eta_0=10^-4`, the physical log progress satisfies

\[
|g-\gamma_*|\le {1\over80},
\qquad g=\log(q/p).
\]

Using the certified `r_*` bracket, Arb gives

\[
\boxed{
 e^{\gamma_*-1/80}>\frac85.
}
\]

Thus a sticky child-to-parent lineage advances its reference/filter frequency by more than `8/5` per generation.

For a Fourier covector transported by the same affine background,

\[
\dot k=-A^Tk,
\qquad S={A+A^T\over2},
\]

hence

\[
{d\over dt}\log|k|=-\widehat k^TS\widehat k
\le\|S\|_{op}.
\]

On the existing low-strain lifetime branch

\[
\int\|S\|_{op}dt\le {1\over30},
\]

Arb gives

\[
\boxed{
{M_{q+1}\over M_q}\le e^{1/30}<\frac{21}{20}.
}
\]

So the cascade scale advances much faster than one old reservoir can Kelvin-stretch on a low-strain episode.

## 2. Galilean cancellation gives a service half-life

The exact LP increment estimate assigns to a low band `M<=N` the square-function contribution

\[
(\beta R_G)^2\left({M\over N}\right)^2\mu_M,
\qquad
\mu_M=M\|u_M\|_2^2.
\]

Therefore one materially reused reservoir has a **per-unit-critical-mass** service coefficient whose ratio between consecutive low-strain generations is at most

\[
\boxed{
\left({21/20\over8/5}\right)^2
=\left({21\over32}\right)^2
={441\over1024}
<\frac12.
}
\]

This is the spectral version of the statement that uniform large-scale sweeping does not create a small-scale velocity increment.

Writing `E_M=||u_M||_2^2`, the same contribution equals

\[
(\beta R_G)^2{M^3\over N^2}E_M.
\]

Hence the service coefficient **per unit physical energy** contracts by

\[
\boxed{
\left({21\over20}\right)^3
\left({5\over8}\right)^2
={231525\over512000}
<\frac12.
}
\]

Consequently, even under the adversarial assumption that this one reservoir may own the entire conserved global energy at every service time,

\[
\mathsf C_q
\le
\mathsf C_0
\left({231525\over512000}\right)^q
<2^{-q}\mathsf C_0,
\]

and

\[
\boxed{
\sum_{q\ge0}\mathsf C_q<2\mathsf C_0.
}
\]

Thus a single low-strain reservoir has finite total future service capacity.  If every efficient block requires a uniform positive increment/source service threshold, one material reservoir can satisfy only finitely many generations.

## 3. "Same reservoir" has a material phase-space fingerprint

Reuse is not declared by putting the same ancestry label on two packets.  For the material affine frame

\[
\dot L=AL,
\qquad
\dot k=-A^Tk,
\]

one has the exact Kelvin/Liouville invariant

\[
\boxed{q=L^Tk=\text{constant}.}
\]

This is the affine-grain covector already used in the Gaussian forcing module.  A later packet whose carrier cannot be obtained from the old `q` by the transported frame is a **relinking/source event**, not free reuse.  The common helical `SL(2)` propagator and covariance ancestry similarly compose rather than reset.

The theorem therefore distinguishes:

1. **material reuse**: same transported phase-space state, to which the half-life applies;
2. **spectral relinking**: nonlinear/new-frequency content, routed to cross/source/fresh bookkeeping;
3. **fragmentation**: many packets, already routed by the `log 2` Bellman / `1/4` cycle theorem;
4. **high-strain tracking**: if a reservoir moves in frequency fast enough to avoid the half-life, the low-strain hypothesis itself fails and the existing strain/source ledger is activated.

## 4. What this closes and what remains

This closes the feared loophole in which one old low-frequency reservoir services an infinite chain of progressively smaller scales while paying only once.  On low-strain material reuse its maximum service capacity decays geometrically, regardless of how its amplitude is replenished, because the global energy cap is inserted only after the kinematic coefficient has decayed.

It does **not** yet control many spatially separated high-frequency grains sampling the same band-limited reservoir at one generation.  That is a band-limited sampling/packing theorem.  Nor does it yet register the near-field pressure-third source.  Those are the next continuum bridges.
