# Wang–NS Triad Diamond

A theorem-driven research repository exploring a rigidity/no-escape architecture for the 3D incompressible Navier–Stokes equations.

The project began as a four-edge helical Fourier “triad diamond” experiment inspired by extremizer/rigidity ideas.  That experiment is now only the historical seed.  The active programme works directly with physical Navier–Stokes observables: positive nonlinear work, smooth moving Fourier roles, objective source ownership, coherent/Moyal service, material transport, critical hard shells, physical tail dissipation, and recursive first-stop laws.

This repository **does not claim a proof of global regularity for 3D Navier–Stokes**.  It contains exact identities, conditional continuum theorems, certified constants, counterexamples to tempting false reductions, and GitHub-Actions stress/certificate evidence for the current proof architecture.

## Start here

1. [`RESEARCH_LEDGER.md`](RESEARCH_LEDGER.md) — the current canonical map of the programme, organized by physical structure rather than discovery order.
2. [`docs/canonical_positive_edge_work_routing.md`](docs/canonical_positive_edge_work_routing.md) — canonical `dW+` routing, stage-zero bad-work loss, signed hard-cell inheritance, and the mixed-fate Young seam.
3. [`docs/continuum_helical_edge_measure_registration.md`](docs/continuum_helical_edge_measure_registration.md) — the signed unordered Fourier/helicity NS edge measure that routing inherits.
4. [`docs/physical_energy_causal_bridge.md`](docs/physical_energy_causal_bridge.md) — why causal weights are actual positive child-energy work, not normalized Duhamel mass.
5. [`docs/critical_shell_service_reentry.md`](docs/critical_shell_service_reentry.md) — generic critical-shell first-stop / own-scale service reentry.
6. [`docs/master_no_escape.md`](docs/master_no_escape.md) — the finite-dimensional/master-side no-escape architecture and its continuum interface.

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
- **Entropy symbols are measure-specific.** Pressure-pair, fresh-scale, high-tail scale/time concentration, and causal ancestry entropies must not be interchanged.

The anti-theorem/forbidden-inference list is maintained in [`RESEARCH_LEDGER.md`](RESEARCH_LEDGER.md#28-anti-theorems-and-forbidden-inferences).

## Current frontier

The programme-level frontier is **mixed genuine-owner recurrence** after the representation/checkpoint/interface quotients and the finite pure high-strain / signed-good generated-HH epoch theorems.  Generic HH is now resolved one physical layer deeper: canonical bad `dW+` terminates as same-law transfer loss, while good `dW+` is only Young-eligible.  The sharp local HH seam is the fate-pure/mixed-fate signed hard-cell Young/Christ handoff (or a genuinely positive coherent kernel), not construction of a new causal law.

The eventual master must combine the surviving physical owners without a synthetic clock, without re-Hahn after analysis coarsening, and without promoting scale-critical quantities into artificial finite resets.  No global-regularity conclusion is claimed.

## Validation and reproducibility

Serious theorem validation in this project is performed through GitHub Actions.  Canonical positive edge-work routing was certified on exact SHA `d78a3f4223c16f541548e7628cd73f70e6fdee6c` by dedicated run `31559241272`, independent adversarial/actual-NS audit `31559242700`, and full causal integration `31559248529`, with continuum/helical compatibility gates on the same SHA.  The same SHA also passed a `22/22` directly affected current/historical workflow sweep after the global performance pass.

Key CI artifacts are stored under `recorded-results/<run-id>/`.

### Experimental engineering discipline

The project does **not** trade experimental severity for runtime.  Serious candidates should be challenged by large stress ensembles, adversarial fixtures, exact provenance checks, and actual Navier--Stokes/Galerkin probes whenever the theorem calls for them.  The engineering obligation is instead to make those experiments as cheap as the mathematics permits.

Profile before buying compute.  Separate invariant geometry/provenance from the parameters that genuinely vary, precompute only what is mathematically immutable, remove repeated replay inside one verified transaction, and keep a slow reference path for differential/equivalence checks.  A performance change must not silently change the sampled law, sample count, tolerances, physical units, signed-before-Hahn order, causal provenance, or routing thresholds.  Small fixed-dimensional kernels should not pay large-array framework overhead merely for convenience; use scalar/vectorized/native C++ or Rust kernels only when profiling shows they are the remaining bottleneck, and bind them to the reference implementation by regression tests.

Performance benchmarks measure the research instrument, not the theorem: record wall time, CPU use and peak memory, but keep theorem certification on the same dedicated/adversarial/actual-PDE GitHub Actions gates.  Faster code is valuable because it permits **more** hard experiments and quicker falsification, not because it licenses fewer of them.  The 2026-08-12 global pass is the reference implementation of this rule: it preserved exact/reference equivalence while cutting the Actions actual-Galerkin routing step from `309s` to `26s` and the integration helical `50,000` stress from `436s` to `65s`, without reducing either experiment.

For exploratory historical experiments, see the older modules/workflows and the history ledger.  For current theorem work, follow the workflow described in [`RESEARCH_LEDGER.md`](RESEARCH_LEDGER.md#34-certification-discipline).
