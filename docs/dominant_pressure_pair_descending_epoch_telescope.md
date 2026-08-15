# Dominant resolved-pressure-pair descending epoch telescope

Status: **DRAFT THEOREM CANDIDATE — independent PR, not certified, not merge-ready.**

This note closes one genuine source-service subepoch using only two physical statements already certified for the **same resolved pressure-pair owner**:

1. supplier-specific shell progress,
   \[
   N_{j+1}\le \frac14 N_j;
   \]
2. on the canonical quarter-dominant face `q_max>=1/4`, an actual full-velocity hard shell with
   \[
   \mu_{j+1}:=N_{j+1}\|P_{N_{j+1}}u\|_2^2
   \ge \frac{80}{c}\Sigma_{P,j}.
   \]

No generic shell registration, event count reset, pressure entropy cost, or artificial clock is used.

## 1. Uniform source floor from the objective first stop

At an objective-source event with variation action lower `A_*`, the exact four-owner cover gives every qualifying pressure owner

\[
\Sigma_{P,j}\ge \sigma_*:=\frac{A_*}{4c}.
\]

If that pressure owner is carried by the resolved positive pair law and an actual pair has normalized mass at least one quarter, the certified pressure atomization theorem gives

\[
\boxed{
\mu_{j+1}\ge \mu_*:=\frac{80\sigma_*}{c}
=\frac{20A_*}{c^2}.
}
\]

The shell is an actual `u` hard shell inherited by contraction from `V=S_{N/4}u`; no inverse low-pass estimate is used.

## 2. Global energy turns shell criticality into a physical frequency floor

Let

\[
\|u(t)\|_2^2\le E_{\rm global}
\]

on the physical interval under consideration. For every selected pressure-pair child shell,

\[
\mu_*
\le N_{j+1}\|P_{N_{j+1}}u\|_2^2
\le N_{j+1}E_{\rm global}.
\]

Therefore

\[
\boxed{
N_{j+1}\ge N_{\min}:=\frac{\mu_*}{E_{\rm global}}.
}
\]

This is not an imposed infrared cutoff. It is a consequence of an actual shell carrying a fixed critical mass while the entire flow has finite kinetic energy.

## 3. Consecutive dominant pressure-pair recursion cannot continue indefinitely

Consider a maximal consecutive epoch in which each pressure-pair child hard shell is the next recursive state's incoming physical frequency. Then the supplier-specific scale theorem gives

\[
N_L\le 4^{-L}N_0
\]

after `L` pressure-pair transitions. But every selected child also obeys `N_L>=N_min`. Hence

\[
4^L\le \frac{N_0}{N_{\min}}
=\frac{N_0E_{\rm global}}{\mu_*},
\]

and therefore

\[
\boxed{
L\le
\left\lfloor
\log_4\!\left(\frac{N_0E_{\rm global}}{\mu_*}\right)
\right\rfloor.
}
\]

With the objective owner floor,

\[
\boxed{
L\le
\left\lfloor
\log_4\!\left(
\frac{N_0E_{\rm global}c^2}{20A_*}
\right)
\right\rfloor.
}
\]

If the argument of the logarithm is below `4`, even one such renewal transition is impossible under the same source floor.

The mechanism is intrinsic: pressure-pair geometry tries to drive the recursive shell downward by at least a factor four per event, while the actual critical shell mass and global kinetic-energy bound forbid the frequency from falling below a fixed physical floor.

## 4. Why this is not a finite-reset argument

The theorem does **not** say that every critical shell consumes a scale-independent amount of energy. That statement is false in general and is explicitly forbidden by the multi-currency anti-reset theorem. Instead, along this one supplier-specific descending lineage, the same shell criticality furnishes a **frequency lower bound**, while the pressure-pair support theorem furnishes a **frequency upper contraction**. The contradiction is geometric, not additive.

It also does not use the generic carrier registration `A=3M/4`; that registration is not progress. The only scale ratio used is the pressure theorem's actual resolved-pair support `N_next/N<=1/4`.

## 5. Diffuse pressure is deliberately excluded

If `q_max<1/4`, the pressure theorem still gives the unconditional entropy-shell conjugacy

\[
\mu_{child}e^{H_2^{pair}}\ge 320\Sigma_P/c,
\]

but the same theorem explicitly states that `H_2^{pair}` is **not** a child-energy causal probability, not a new causal fate, and not a terminal transfer cost. This draft therefore does not convert diffuse pressure entropy into Bellman cost.

The dominant theorem closes only the quarter-dominant resolved-pair subepoch. Diffuse pressure renewal remains in the genuine source/service frontier.

## 6. Consequence for mixed genuine-owner recurrence

The native recursive source class can now be refined further:

- local/viscous objective source -> existing strain/critical-dissipation root;
- quarter-dominant resolved pressure-pair source -> finite descending pressure epoch by this theorem;
- SGS high-frequency -> physical hard-tail owner;
- selected-interface `Xi` -> existing transfer cost;
- old-pool capacity -> inherited reservoir state, never a new causal budget;
- **remaining hard source-service cases:** diffuse resolved pressure-pair renewal and fresh SGS critical-scale renewal (plus any realized pressure-SGS service not yet descended through its own owner theorem).

This does not yet close arbitrary alternation with strain and HH/high-tail owners. It narrows the source component by one more native theorem.

No Navier--Stokes global-regularity claim is made.
