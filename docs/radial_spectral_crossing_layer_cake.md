# Radial spectral crossing layer cake

Status candidate:

`EXACT_RADIAL_SPECTRAL_CROSSING_LAYER_CAKE__PHYSICAL_HELICAL_DONOR_RECIPIENT_RADII__TAIL_INWARD_OUTWARD_FLUX__TRUNCATED_LOG_RADIUS_LAYER_CAKE__NO_REHAHN_DYADIC_PROGRESS_OR_GROSS_CROSSING_BUDGET`

This theorem is downstream of the certified cyclic donor kernel and helical mode-set energy continuity theorem.  It introduces **no new work law**.  It asks one narrower physical question: when the already-existing same-time energy flow goes from one actual helical Fourier mode to another, which Fourier spheres does that energy transfer really cross?

## 1. Physical input

A cyclic donor atom is

\[
(d,r,m),\qquad m>0,
\]

where `d=(k_d,s_d)` is the actual helical donor mode, `r=(k_r,s_r)` is the actual helical recipient mode, and `m` is the physical donor/recipient work mass inherited from the canonical `dW^- -> dW^+` disintegration.

Nothing here re-Hahn-splits a tail or shell.  The positive flow atoms already exist before a radius is chosen.

Write

\[
\rho_d=|k_d|,\qquad \rho_r=|k_r|.
\]

For a radial exterior

\[
H_R=\{(k,s):|k|>R\},
\]

an atom has exactly four physical fates:

- low→low: `rho_d<=R` and `rho_r<=R`;
- high→high: `rho_d>R` and `rho_r>R`;
- low→high: `rho_d<=R<rho_r`;
- high→low: `rho_r<=R<rho_d`.

The last two are the actual radial boundary crossings.

## 2. Exact tail boundary law

Define

\[
\Phi_\uparrow(R)=\mathcal M\{\rho_d\le R<\rho_r\},
\]

and

\[
\Phi_\downarrow(R)=\mathcal M\{\rho_r\le R<\rho_d\}.
\]

Let `I_>(R)` be high→high internal traffic.  The already-canonical recipient and donor marginals restricted to `H_R` are

\[
W_{>R}^+=I_>(R)+\Phi_\uparrow(R),
\]

\[
W_{>R}^-=I_>(R)+\Phi_\downarrow(R).
\]

Therefore

\[
\boxed{
W_{>R}^+-W_{>R}^-=\Phi_\uparrow(R)-\Phi_\downarrow(R).
}
\]

The high→high traffic is physically real, but cancels exactly from tail divergence.  Low→low traffic never enters the tail balance.

This is the radial specialization of the certified mode-set divergence theorem.  It does not make `[W_{>R}^+-W_{>R}^-]_+` into another causal Hahn law.

## 3. Actual Navier–Stokes tail continuity

For the physical tail energy

\[
E_{>R}(t)=\sum_{|k|>R}\sum_{s=\pm1}|a_{k,s}(t)|^2,
\]

Navier–Stokes gives

\[
\boxed{
E_{>R}(t_1)
+D_{>R}[t_0,t_1]
+\int_{t_0}^{t_1}\Phi_\downarrow(R,t)\,dt
=
E_{>R}(t_0)
+\int_{t_0}^{t_1}\Phi_\uparrow(R,t)\,dt,
}
\]

with

\[
D_{>R}=2\nu\int_{t_0}^{t_1}\sum_{|k|>R,s}|k|^2|a_{k,s}(t)|^2dt.
\]

This says precisely what increases or decreases radial tail stock:

- low→high crossing supplies the tail;
- high→low crossing removes nonlinear tail stock;
- viscosity removes tail stock;
- high→high nonlinear circulation does neither.

It does **not** say the gross crossing traffic is globally finite.

## 4. Truncated logarithmic layer cake

Fix `0<R0<R1` and define the clipped log potential

\[
\phi_{R_0,R_1}(\rho)
=
\log\frac{\min(\max(\rho,R_0),R_1)}{R_0}.
\]

For one atom,

\[
\int_{R_0}^{R_1}{\bf1}_{\rho_d\le R<\rho_r}\frac{dR}{R}
=
[\phi(\rho_r)-\phi(\rho_d)]_+,
\]

and

\[
\int_{R_0}^{R_1}{\bf1}_{\rho_r\le R<\rho_d}\frac{dR}{R}
=
[\phi(\rho_d)-\phi(\rho_r)]_+.
\]

After integrating the already-existing physical flow measure,

\[
\boxed{
\int_{R_0}^{R_1}\Phi_\uparrow(R)\frac{dR}{R}
=
\int [\phi(\rho_r)-\phi(\rho_d)]_+\,d\mathcal M,
}
\]

\[
\boxed{
\int_{R_0}^{R_1}\Phi_\downarrow(R)\frac{dR}{R}
=
\int [\phi(\rho_d)-\phi(\rho_r)]_+\,d\mathcal M.
}
\]

Subtracting gives

\[
\boxed{
\int_{R_0}^{R_1}(\Phi_\uparrow-\Phi_\downarrow)\frac{dR}{R}
=
\int\phi(|k|)\,dW^+
-
\int\phi(|k|)\,dW^-.
}
\]

This is a layer-cake identity of the same physical donor/recipient measure.  `dR/R` appears because a transfer crosses exactly the logarithmic interval between its two physical radii; it is not an analyst-created clock.

## 5. Infinite-range finite-atom identity and continuum caution

For a finite physical flow law, every mode radius is positive and the full log moment is finite, so

\[
\boxed{
\int_0^\infty\Phi_\uparrow(R)\frac{dR}{R}
=
\sum_{d\to r}m_{dr}\log_+\frac{|k_r|}{|k_d|},
}
\]

\[
\boxed{
\int_0^\infty\Phi_\downarrow(R)\frac{dR}{R}
=
\sum_{d\to r}m_{dr}\log_+\frac{|k_d|}{|k_r|}.
}
\]

The signed identity is

\[
\boxed{
\int_0^\infty(\Phi_\uparrow-\Phi_\downarrow)\frac{dR}{R}
=
\sum_{d\to r}m_{dr}\log\frac{|k_r|}{|k_d|}.
}
\]

A continuum infinite-range extension requires an explicit finite logarithmic moment.  The theorem does **not** infer that moment merely from local Radon variation.  The truncated identity is the unconditional local statement.

## 6. No minimum radial progress

Physical energy transfer does not imply radial scale progress.

There is a regular closed helical triad with three distinct wavevectors of equal magnitude.  With a nontrivial helicity assignment and physical Waleffe phase, one mode loses energy and another gains it while

\[
|k_d|=|k_r|.
\]

The donor work is positive, yet

\[
\log\frac{|k_r|}{|k_d|}=0,
\]

and it crosses no Fourier sphere.

Thus there is no universal dyadic step, no positive minimum log increment, and no event-count interpretation of radial crossing.

## 7. Radial displacement is not single-edge Young progress

The existing one-edge quantity is

\[
g_e=\log_+\frac{|k_{child}|}{\max(|k_{parent,1}|,|k_{parent,2}|)}.
\]

Radial donor displacement is instead

\[
\log\frac{|k_r|}{|k_d|}
\]

for one particular physical energy donor `d` and recipient `r` in the cyclic closed triad.

They answer different questions:

- `g_e` asks whether one child is above **both interaction parents** and enters the certified `A J c` progress identity;
- radial displacement asks which Fourier spheres the actual energy transfer between one donor and one recipient crosses.

A positive recipient can have two energy donors, so there can be two distinct radial displacements attached to one child edge.  The theorem therefore never identifies radial action with `J_e`, Young/Christ saturation, or the signed-good generated-HH scale ratio.

## 8. Relationship to the high-tail route

The theorem clarifies but does not replace the existing high-tail machinery.

For a fixed radial exterior, high→high circulation does not contribute to nonlinear tail divergence.  Therefore a genuine nonlinear increase of tail stock requires low→high boundary crossing at that radius.

But this theorem alone does not supply:

- the high-tail supplier's `M/N>=2` progress;
- the common causal unit `N dW` used by that route;
- ultraviolet locality;
- a natural-time window;
- a critical-shell reentry lower bound.

Those remain separate certified PDE theorems with their own hypotheses and units.

## 9. Numerical certification policy

The actual-PDE probe keeps two readings separate on the same evolved Fourier–Galerkin state.

1. **Full tail reading.**  It computes the actual radial-tail energy, signed nonlinear tail work, and viscous dissipation directly from the Fourier state at every RK4 output time and checks the integrated Navier–Stokes tail stock identity.
2. **Selected cyclic crossing reading.**  It reads the actual evolving helical coefficients of the certified closed triad from the same state, reconstructs the three-root cyclic donor law, and restricts those donor atoms by the same radial boundary.

The selected triad is not claimed to exhaust the full tail.  Cross-FFT comparison uses one common Galerkin cutoff represented on grids `24` and `28`; native stock/flow throughputs are numerical error envelopes only.

The sign-reversed physical initial condition is retained.  At the canonical radius `R=8`, the ordinary signed-good triad begins low→high, while global sign reversal reverses cubic work and produces high→low crossing on the same geometry.

## 10. Scope

If certified, this theorem closes **radial spectral-crossing registration and its logarithmic layer cake**.

It does not prove that total upward crossing is finite over a singular-time history.  It does not prove a positive average radial drift on a recursive lineage.  It does not convert a layer-cake integral into a Bellman clock.  It does not close degenerate Young/Christ margin, generic HH recurrence, mixed-owner recurrence, or Navier–Stokes global regularity.

The physically honest next question after certification would be whether this radial boundary law can be interfaced with the already-certified hard-tail dissipation/locality route so that true upward crossing—not high→high circulation—is the quantity entering ultraviolet continuation.
