# Radial spectral crossing layer cake

Status: **CERTIFIED**.

`EXACT_RADIAL_SPECTRAL_CROSSING_LAYER_CAKE__PHYSICAL_HELICAL_DONOR_RECIPIENT_RADII__TAIL_INWARD_OUTWARD_FLUX__TRUNCATED_LOG_RADIUS_LAYER_CAKE__NO_REHAHN_DYADIC_PROGRESS_OR_GROSS_CROSSING_BUDGET`

Exact certified implementation/provenance SHA: `667e687cb740a77df944753c575f581abad14199`.

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
The finite full-log marginal identity and uniform-dilation covariance are likewise guarded on the total physical log-variation scale, never on a signed action that may vanish by upward/downward cancellation.

The sign-reversed physical initial condition is retained.  At the canonical radius `R=8`, the ordinary signed-good triad begins low→high, while global sign reversal reverses cubic work and produces high→low crossing on the same geometry.

## 10. Scope

This theorem closes **radial spectral-crossing registration and its logarithmic layer cake**.

It does not prove that total upward crossing is finite over a singular-time history.  It does not prove a positive average radial drift on a recursive lineage.  It does not convert a layer-cake integral into a Bellman clock.  It does not close degenerate Young/Christ margin, generic HH recurrence, mixed-owner recurrence, or Navier–Stokes global regularity.

The physically honest next question is whether this radial boundary law can be interfaced with the already-certified hard-tail dissipation/locality route so that **true low→high boundary supply**, not high→high circulation, is the quantity entering ultraviolet continuation.  That next theorem must preserve the existing high-tail common causal unit and distinguish incoming crossing, outgoing crossing, inherited tail stock, internal high-frequency circulation, and viscous dissipation before any estimate is attempted.

---

## 11. Certification evidence

Exact certified theorem/provenance SHA: `667e687cb740a77df944753c575f581abad14199`.

Dedicated GitHub Actions run `31606829692` completed successfully.  The full theorem suite passed `884` tests; focused radial tests passed `11`, donor-flow dependency tests passed `26`, and focused evolved-NS tests passed `9`.  The `75,000`-triad physical radial stress recorded `56,351` upward and `56,306` downward crossing cases, with `37,657` laws exhibiting both directions across sampled radii.  Worst radial partition, tail-divergence, truncated layer-cake, full-log-marginal and dilation residuals were respectively `0`, `1.977e-16`, `3.460e-16`, `2.036e-11`, and `2.731e-11` on their native physical work/log-variation scales.  The equiradial anti-theorem carried physical work `0.0777635049765` with exactly zero radial action; `131,168` sampled donor atoms explicitly disagreed with the single-edge Young/progress observable.

The same dedicated run evolved the actual `2/3`-dealiased incompressible Fourier--Galerkin Navier--Stokes system at common cutoff `7` on FFT grids `24,28`.  At `R=8`, the ordinary `64`-step branch had integrated selected upward/downward crossing `0.00848754122477 / 0`, while the sign-reversed branch had `0 / 0.00848518440777`.  Their full-tail stock/work/viscosity interval residuals were `3.121e-11` and `1.006e-11`; cross-FFT tail observables remained below `5.6e-16`.  Amplitude adversaries `0.5` and `2.0` also passed.

Independent audit `31606829776` completed successfully.  Its algebra lane passed `37` focused theorem/dependency tests and ran `100,000` physical radial laws: `75,047` upward cases, `74,959` downward cases, and `50,006` laws carrying both directions across sampled radii.  Worst tail-divergence and truncated layer-cake residuals were `2.095e-16` and `3.454e-16`; full-log-marginal and dilation residuals were `3.969e-11` and `5.947e-11`.  The equiradial work/action remained `0.0777635049765 / 0`, and `174,972` donor atoms separated radial displacement from single-edge progress.

The independent evolved-NS audit used amplitude `1.3`, `96` RK4 steps, duration `0.0012`, common cutoff `7`, and FFT grids `24,28`.  On the upward branch the full tail changed `4.0053 -> 4.01741597661`, with signed tail nonlinear work `0.0453411766896`, viscous loss `0.0332252003104`, interval residual `5.592e-11`, and selected integrated upward crossing `0.0223643661196`.  Under sign reversal the same physical geometry produced downward crossing `0.0223546449575`; the tail changed `4.0053 -> 3.92730005969`, signed nonlinear work was `-0.0451218090429`, viscosity `0.0328781311252`, and interval residual `3.526e-11`.  Radius adversaries `R=6.5,8,10` also passed.  The full-tail reading and selected-triad crossing remained distinct observables throughout.

Full physical-energy causal integration `31606829708` completed successfully on the same SHA.  It passed the same `884`-test suite and ran `50,000` radial laws immediately after mode-set continuity and before complex Young.  The integration stress recorded `37,651` upward cases, `37,510` downward cases, `25,161` both-direction laws, worst tail-divergence residual `1.975e-16`, truncated layer-cake residual `3.307e-16`, and the same equiradial zero-action anti-theorem; the complete downstream causal spine remained green.

Stored artifact trees:

- `recorded-results/31606829692/`, deterministic tree digest `sha256:66c925e0b9ebafc91e72a9097a943ef437eaed701620ef6d0af7c09ec80299f4`;
- `recorded-results/31606829776/`, deterministic tree digest `sha256:3652cd5834ddd3aa75319408093fc2d2c31a5c3cdace1bdc480a5d95fc9bc174`;
- `recorded-results/31606829708/`, deterministic tree digest `sha256:a4bf0a8377f818f2628e25d9d0d850652c842ad86adc02181b1587ef46197117`.

Failure lineage is retained.  Initial candidate `fb8e95132b9ce8da566e3b338cee0676b6be5dc9` already passed both dedicated and independent actual-NS radial jobs, but one of `11` focused algebra tests failed before stress: the equiradial physical fixture called the Waleffe coupling with `(k_1,k_2,-k_0)` instead of the actually closed triple `(k_1,k_2,k_0)`.  Repair `697bb77d406d28d51efc396af1c92b948bf8930d` changed that one fixture orientation and then passed dedicated `31605985715`, independent audit `31605985621`, and integration `31605985585` without changing any PDE law or theorem identity.  Final SHA `667e687...` added hard source guards for the finite full-log marginal and dilation covariance on the total physical log-variation envelope; it changed no crossing law, sample count, PDE evolution, threshold semantics, or causal unit, and reran the full three-gate certification successfully.

This theorem does not prove finiteness of total upward crossing over a singular-time history, positive average radial drift, or recurrence termination.  Its next sharp interface is the actual hard-tail energy law: separate low→high supply from high→high circulation and then ask what the existing physical dissipation/locality machinery forces.  No Navier--Stokes global-regularity claim is made.
