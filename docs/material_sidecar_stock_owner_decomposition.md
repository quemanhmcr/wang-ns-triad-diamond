# Material-sidecar / inherited-stock owner decomposition

Exact theorem SHA: `c39f70ed3e9f4edfc08e23b996ba872d91c63ed6`.

## Physical question

A same physical Fourier--helical carrier can persist across a no-first-stop interval while material bookkeeping changes.  The inherited energy on that carrier is already known to be between-time stock when the inherited-stock gate is met.  The remaining question is what, physically, the material sidecars are allowed to mean.

They are not all one currency.  The theorem keeps the physical carrier stock, signed nonlinear work, and material observation-set bookkeeping separate.

## Intrinsic membership rereading

For an intrinsic material-membership update, the positive service edge weights are held fixed and only the OO/ON/NN ownership labels are reread.  Their partition changes, but their total does not:

`service_before = service_after`,

and the service created by relabeling is exactly zero.

Thus intrinsic membership change is retained as **zero-charge provenance**.  It neither creates nonlinear work nor clones the inherited carrier stock.

## Selected-family switching

For a selected family changing from `S_old` to `S_new`, the exact sidecar currency is the Moyal symmetric-difference energy

`R_switch = sum_{C in S_old Δ S_new} E_C`,

with the native selected-energy jump law

`|E(S_new)-E(S_old)| <= R_switch`

up to the already-certified floating representation tolerance of the source theorem.  The inherited-stock certificate therefore carries the numerical `selected_family_switch_energy` forward and downstream use must bind it back to this exact Moyal certificate.

`R_switch` is **selected-service boundary / ancestry currency**.  It is not any of the following:

- signed or positive Navier--Stokes nonlinear work `dW` / `dW+`;
- inherited physical carrier stock `E0`;
- a second carrier coefficient impulse;
- smooth physical relink `K_phys`;
- a reason to Hahn-split again downstream.

The inherited stock remains one physical charge, not one copy per sidecar.

## Same-state anti-theorem

The decisive anti-theorem does not require a surrogate evolution.  Take the same coherent state at both observations and change only the selected family.  Every cell energy increment is then exactly zero **before** any positive/negative split:

`Delta E_C = 0` for every cell,

so

`P_plus = P_minus = 0`

and total coherent state energy is unchanged.  Nevertheless, whenever the two selected families differ on cells of positive energy,

`R_switch > 0`.

The same `R_switch` is reconstructed directly from the Moyal cell energies by summing over `S_old Δ S_new`.  Therefore a positive selected-family boundary charge by itself cannot certify physical work or a physical generation event.  **The boundary charge itself has zero generation depth.**

This is deliberately narrower than saying that every material event is zero-depth.  A genuine material/source recurrence remains possible, but it requires independent physical service/source evidence.

## Relation to the selected-cell no-escape law

The existing selected-cell identity is

`P_plus <= E_final + P_minus + R_switch`.

Its `relink_symmetric_difference` candidate reads exactly this Moyal boundary currency.  The name is historical selected-set bookkeeping terminology.  It must not be identified with smooth physical `K_phys` relink.  Deciding how this zero-generation-depth boundary term and any independently witnessed material/source service enter the central/joint first-stop logic is a separate integration theorem.

## Actual Navier--Stokes referee

No artificial “material PDE” is introduced.  The sidecar theorem is exact Moyal/service algebra, while the physical carrier ontology is checked against the already-certified dealiased Fourier--Galerkin Navier--Stokes trajectory on the same helical carrier `(7,6,5)`, helicity `+1`, at FFT representations `24,28`.

Dedicated amplitude `1`:

- stock `1 -> 1.00240956249`;
- gross nonlinear `dW+ / dW- = 0.0100088487847 / 0.000991283878665`;
- viscous dissipation `0.00660800248318`;
- maximum modal continuity residual `6.253e-11`;
- maximum global energy-balance residual `2.963e-12`.

Independent amplitude `0.5`, opposite phase:

- stock `0.25 -> 0.245853600709`;
- gross nonlinear `dW+ / dW- = 0.000182591040316 / 0.00187455879673`;
- viscous dissipation `0.00245443154604`;
- maximum modal continuity residual `4.761e-11`;
- maximum global energy-balance residual `2.614e-11`.

Thus inherited stock remains compatible with real nonlinear exchange and viscosity.  The theorem changes representation typing, not the Navier--Stokes dynamics.

## Certification

Exact-SHA gates on `c39f70ed3e9f4edfc08e23b996ba872d91c63ed6`:

- dedicated `31678779916`: success; `75,000` exact sidecar/stock states and actual NS referee;
- independent audit `31678779930`: success; `100,000` states and opposite-phase actual NS audit;
- physical-energy causal integration v2 `31678801395`: success; `948` tests, master `20,000` episodes with worst margin `0.0`, `61` integration nodes with `58` exact-content reuse and `3` execution, all `61` outputs present and all `58` reused outputs byte/digest verified.

The `75k / 100k` stresses reported zero membership-charge violations, zero Moyal-charge residual, zero inherited-stock mass residual, zero same-state increment residual, and zero smooth-`K_phys` identifications.  They exercised `25,000 / 33,333` selected-family `relink_symmetric_difference` service branches and `50,000 / 66,666` positive-`R_switch` same-state anti-theorem cases.

No global-regularity claim is made.
