# Wang–NS Triad Diamond

A theorem-driven research repository exploring a rigidity/no-escape architecture for the 3D incompressible Navier–Stokes equations.

The project began as a four-edge helical Fourier “triad diamond” experiment inspired by extremizer/rigidity ideas.  That experiment is now only the historical seed.  The active programme works directly with physical Navier–Stokes observables: positive nonlinear work, smooth moving Fourier roles, objective source ownership, coherent/Moyal service, material transport, critical hard shells, physical tail dissipation, and recursive first-stop laws.

This repository **does not claim a proof of global regularity for 3D Navier–Stokes**.  It contains exact identities, conditional continuum theorems, certified constants, counterexamples to tempting false reductions, and GitHub-Actions stress/certificate evidence for the current proof architecture.

## Start here

1. [`RESEARCH_LEDGER.md`](RESEARCH_LEDGER.md) — the current canonical map of the programme (~1,100 lines), organized by physical structure rather than discovery order.
2. [`docs/smooth_quadratic_carrier_interface.md`](docs/smooth_quadratic_carrier_interface.md) — the native `Q^2` energy law for genuinely smooth propagated PDE carriers.
3. [`docs/resolved_interface_donor_quotient.md`](docs/resolved_interface_donor_quotient.md) — the complementary hard-event donor/circulation quotient.
4. [`docs/high_tail_natural_window_reentry.md`](docs/high_tail_natural_window_reentry.md) — the completed high-tail continuation theorem.
5. [`docs/physical_energy_causal_bridge.md`](docs/physical_energy_causal_bridge.md) — why causal weights are actual positive child-energy work, not normalized Duhamel mass.
6. [`docs/objective_source_routing_compiler.md`](docs/objective_source_routing_compiler.md) — current source-owner routing.
7. [`docs/critical_shell_service_reentry.md`](docs/critical_shell_service_reentry.md) — generic critical-shell first-stop / own-scale service reentry.
8. [`docs/master_no_escape.md`](docs/master_no_escape.md) — the finite-dimensional/master-side no-escape architecture and its continuum interface.

For the full chronological development, superseded routes, fixture failures, and detailed CI provenance, read [`docs/history/RESEARCH_LEDGER_history_through_2026-08-10.md`](docs/history/RESEARCH_LEDGER_history_through_2026-08-10.md).

## Current physical spine

The canonical high-level chain is

```text
actual positive nonlinear work
  -> exact physical owner / first-hit cause
  -> critical hard shell or coherent service
  -> own-scale recursive first stop
  -> renewed service / regeneration / t=0 / named physical owner
```

The high-frequency branch is now more explicit:

```text
high coherent service
  -> physical tail dissipation/energy
  -> inherited critical shell OR actual regeneration work
  -> common-unit HH/interface causal ownership (N dW)
  -> Fourier UV locality paid by D_tail
  -> sliding M-natural-time concentration
  -> actual critical hard-shell reentry
```

For the sliding-window branch,

\[
\frac{\sqrt{\mu_{\rm win}}}{p_{\rm scale}p_{\rm time}}
\ge
\frac{\nu D_{\rm tail}}
{48c\sqrt\pi\,R\,N E_{\rm global}}.
\]

The selected high-tail shell also has genuine support geometry

\[
M/N\ge2,\qquad T_M/T_N\le1/4.
\]

No packet persistence or observer-chosen time bins are used in this step.

## Rules that matter

- **Causal probability** comes from actual positive physical work.  Duhamel supplies support/adjoint identities, not an interchangeable probability law.
- **Representation freedom is quotiented before physics is charged.** Common affine/Kelvin motion, cutoff repartition, coherent-cell refinement, and pure material-label changes do not create new physical currencies by themselves.
- **Exact ties remain joint causes.** There is no canonical lexicographic theorem priority.
- **Critical scale-normalized events recurse.** `NE`, `D_V`, and similar `O(1)` scale-critical quantities are not automatically finite global reset budgets.
- **Hard projectors are event readers; smooth roles propagate between events.**
- **Smooth carrier energy is read at `Q^2`.** A non-idempotent `Q` is completed by a square partition; `I-Q` is not its energy complement.
- **Coefficient obstruction locates an interval, not an owner.** Actual carrier energy and native physical work must reenter the energy gate before causality is assigned.
- **Entropy symbols are measure-specific.** Pressure-pair, fresh-scale, high-tail scale/time concentration, and causal ancestry entropies must not be interchanged.

The anti-theorem/forbidden-inference list is maintained in [`RESEARCH_LEDGER.md`](RESEARCH_LEDGER.md#28-anti-theorems-and-forbidden-inferences).

## Current frontier

The main remaining programme-level task is **final continuum master assembly**: combine the now-certified source, shell, service, material, high-strain, and high-tail supplier routes into one measurable recursive no-escape theorem without double counting physical owners, without a synthetic common clock, and without promoting scale-critical currencies into artificial finite resets.

Supplier-specific signed-good geometry and the exact global closure of all recursive routes remain separate issues.  No global-regularity conclusion is claimed.

## Validation and reproducibility

Serious theorem validation in this project is performed through GitHub Actions.  The smooth quadratic-carrier theorem was certified by dedicated run `31401197668` and full causal integration `31401197364` on exact theorem SHA `54bf268fd363a3778627e281756482c64f37ae97`.

Key CI artifacts are stored under `recorded-results/<run-id>/`.

For exploratory historical experiments, see the older modules/workflows and the history ledger.  For current theorem work, follow the workflow described in [`RESEARCH_LEDGER.md`](RESEARCH_LEDGER.md#34-certification-discipline).
