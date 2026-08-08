# Single-charge physical branch compiler: causal quotient before ledger projection

This note addresses the current frontier: every transfer-selected efficient
Navier--Stokes block must be routed into the physical multi-currency master
without paying the same physical work twice under different theorem names.

It does **not** claim that the continuum recursive extraction is closed.  The
main new point is structural: the object which can be canonically partitioned is
the **positive physical child-transfer measure after causal quotient**, not a
naive list of overlapping block predicates.

## 1. Why a raw mutually-exclusive case tree is not physical

The current theorem modules produce observables which are not algebraically
exclusive.  A block may simultaneously have resolved strain, differentiated SGS
source, Hermite dephasing and viscous dissipation.  Those are often different
views of one causal event.  Conversely, two genuinely independent events can
occur at the same physical time.

Therefore the rule

\[
\text{``test theorem A, else theorem B, else theorem C''}
\]

is not canonical unless the order is itself derived from the evolution.
The correct preliminary operation is a **causal quotient**:

two theorem manifestations are identified whenever the existing PDE identity
shows that one is the certified cause/consequence of the other and no new
physical work has been introduced.

Examples:

- H1/H3 dephasing caused by pressure/SGS/viscous source is one source event;
- high resolved strain and the forced lower bound for `D_V` are one critical
  dissipation event;
- a large physical covariance deformation and the strain/source which causes it
  are one event; changing only the nearby analysis covariance is instead a
  representation `Xi`;
- pair rescue and the reuse/cycle currency to which it is routed are one causal
  chain, not two payments;
- a Duhamel `R_class` label delegates to its underlying source/interface term and
  is never a new currency;
- reaching `t=0` is an absorbing boundary event, not a fresh interior packet.

After this quotient, the compiler acts on **causal roots**, not theorem names.

## 2. Stage zero: exact physical-transfer/interface split

Let `dT_B` be the positive physical child-transfer measure of one selected
block on the signed-good coherent core.  The coherent localization operators
satisfy

\[
\sum_C A_C=I
\]

and reconstruct the trilinear work exactly.  Hence there is no
continuous-to-discrete synthesis measure.

At recursive depth `j`, the physical defect moat deletes a selected set of
cross-component interactions.  This gives the exact first split

\[
\boxed{
 d\mathcal T_B
 =d\mathcal T_{\Xi,B}+d\mathcal T_{B}^{\rm ret},
}
\]

where `dT_Xi` is **actual omitted physical cross-cell transfer**, with

\[
\eta_{{\rm cross},j}
\le M_j^{-1}+2\overline{\mathcal D}/R_j.
\]

The retained measure is the only measure passed to the causal branch compiler.
Thus coherent relinking in a retained connected component and omitted cross-cell
`Xi` are disjoint by construction.

Frequency-symbol and small covariance-representative errors are separate
representation operators.  They may coexist with physical cross-cell transfer,
but they are added once to different `Xi` subledgers; there is no duplicated
trilinear term.

## 3. Transfer loss is an absorbing selection gate

A fixed multiplicative loss is not an interior source branch.  It is already a
deficit of the positive transfer law used to select the block.  Therefore if

\[
1-R_B\ge\delta_*>0
\]

or a certified Hodge/entropy/sideband endpoint already gives the fixed
multiplicative cost, the retained block stops immediately in

\[
\boxed{\mathcal C_{\rm mult}.}
\]

No subsequent source, strain or ancestry observable is needed to pay this same
retained transfer mass.  They may remain diagnostics, but they are not additional
primary charges.

This priority is not a proof-order convention: it is the fact that the conserved
object being recursively transported has already lost a definite fraction before
one asks how the surviving near-extremal part evolves.

## 4. First causal defect on the low-transfer-cost branch

Assume the fixed transfer-loss gate did not fire.  For every **causal root** `r`
provided by the existing physical theorems, let `tau_r` be the first physical
time at which its certified threshold is reached.  The relevant roots are

\[
\begin{aligned}
&\text{resolved source/SGS service},\\
&\text{high-strain critical dissipation},\\
&\text{material relink/new coherent ancestry},\\
&\text{causal reuse/Renyi pair or entropy},\\
&\text{intrinsic H1/H3 sideband cost not already owned by a source},\\
&\text{a genuinely uniform globally bounded resource},\\
&\text{the initial boundary}.
\end{aligned}
\]

Define

\[
\boxed{
\tau_B=\inf_r\tau_r.
}
\]

All later consequences are diagnostic only.  In particular, if a source event
occurs first and later produces H1 dephasing, coherent increment service and new
mass, the source event owns that causal chain until the compiler deliberately
recurses into its downstream sticky ancestry.  The downstream theorems certify
that the source cannot disappear; they do not create extra payment for the same
source weight.

If no causal root is hit on the slab, the block is admissible only if the
service-or-flat theorem certifies the Kelvin-flat regime.  Then the destination
is

\[
\boxed{\mathcal C_{\rm flat}}
\]

and the physical flat-episode theorem supplies barycentric erosion
`kappa_0>0.17` at `tau=1/100`.

A retained low-cost block with neither a first causal defect nor certified
Kelvin flatness is **uncompiled**.  It is a continuum-bridge gap, not a new
currency.

## 5. Exact ties: partition the measure, do not invent a theorem priority

There is no Navier--Stokes conservation law saying two independent causal clocks
cannot hit at exactly the same time.  A literal one-block/one-label theorem would
therefore require an arbitrary tie breaker unless an additional transversality
theorem were proved.

The natural symmetric replacement is measure-theoretic.  On the exact first-time
tie set let `a_r>=0` be the positive stopping weights supplied by the physical
observables.  Put

\[
\lambda_r={a_r\over\sum_s a_s}.
\]

Then

\[
\boxed{
 d\mathcal T_{r,B}=\lambda_r\,d\mathcal T_B^{\rm ret},
 \qquad
 \sum_r d\mathcal T_{r,B}=d\mathcal T_B^{\rm ret}.
}
\]

Duplicate theorem manifestations of the **same** causal root are summed before
this normalization; they cannot create a second branch.  Independent roots may
split a tie measure, but the total conserved transfer mass is still charged
exactly once.

Thus the canonical statement is stronger and more physical than a lexicographic
case list:

\[
\boxed{
 d\mathcal T_B
 =d\mathcal T_{\Xi,B}
 +d\mathcal T_{{\rm mult},B}
 +d\mathcal T_{{\rm flat},B}
 +d\mathcal T_{{\rm sticky},B}
 +d\mathcal T_{{\rm reuse},B}
 +d\mathcal T_{{\rm side},B}
 +d\mathcal T_{{\rm src},B}
 +d\mathcal T_{{\rm diss},B}
 +d\mathcal T_{{\rm reset},B}
 +d\mathcal T_{0,B},
}
\]

with every term nonnegative and the total mass exact.

The pieces are mutually singular away from exact independent tie sets.  On a tie
set they are an exact Radon--Nikodym partition rather than duplicated measures.
This is the precise sense in which the compiler is **single charge**.

## 6. The Duhamel measure is a causal kernel, not another copy of `dT`

For a generated adjoint node, the existing theorem gives

\[
 d\Gamma_c(t,\alpha)
 =\big[\Re(e^{-i\vartheta_c}dI_{HH,\alpha})\big]_+.
\]

This is positive and its total dominates the generated amplitude.  But the
adjoint theorem explicitly warns that `dGamma` is an **amplitude-generation
measure**, not automatically the positive physical child-energy transfer law.

Therefore the compiler does **not** write

\[
d\mathcal T+d\Gamma
\]

as if the two objects were one physical measure.  Instead, after a generated
node, normalized `dGamma` is a **conditional parent-pair kernel** for the next
causal layer.  Common Kelvin transport is already in the propagator and the
canonical material label

\[
\zeta=(L^{-1}X/2,L^Tk)
\]

is used on both sides.

To feed this generated layer into the transfer-weighted Shannon/Renyi master
currency, one still needs the continuum bridge

\[
\boxed{
K_\Gamma(\text{parent pair}\mid\text{child})
\equiv
K_{\rm phys}(\text{parent pair}\mid\text{child})
}
\]

on the same selected material law, up to a discrepancy already routed once to a
multiplicative/cross branch.

That identification is **not proved by this note**.  This is an important reason
not to claim a universal one-measure closure prematurely.

## 7. Compiler-to-master projection

The compiler has a fine physical destination and a coarser master disposition.
They are not the same object.

| Physical currency | Master disposition | Rule |
|---|---|---|
| multiplicative transfer loss | `N_T` | fixed transfer cost |
| Kelvin-flat barycentric erosion | `N_F` | consumes flat potential |
| coherent sticky/new ancestry | recurse | critical causal architecture |
| Renyi/Shannon reuse entropy or pair mass | `N_T` when certified fixed cost | multiplicative ancestry cost |
| terminal H1/H3 sideband curvature cost | `N_T` | only the direct deficit endpoint; source-owned dephasing is not charged again |
| resolved source / SGS service | recurse | source weight must collide into ancestry/dissipation/transfer |
| critical `D_V` or critical mass | recurse | **not** a finite reset from energy alone |
| globally bounded uniform resource | `N_A` | only with scale-independent threshold in that same global resource |
| initial boundary | absorbing terminal | not an interior reset |
| `Xi` | `Xi` | summable omitted/representation measure |

In particular the compiler refuses the implication

\[
D_V\ge d_0\quad\Longrightarrow\quad N_A\lesssim D_{global}/d_0,
\]

because the physical dissipation price is `nu d_0/N` and is geometrically
summable along a scale cascade.  The same objection applies to fresh critical
mass `NE`.

Once every `recurse` branch has itself been causally expanded until it reaches a
multiplicative endpoint, a boundary endpoint, a summable `Xi`, or a genuinely
uniform reset, the existing multi-currency theorem applies with no duplicate
resource use.

## 8. Initial boundary is an absorbing state

Asynchronous Duhamel synchronization gives finite backward progress, so an
interior causal chain either stops earlier or reaches `t=0` after finite depth.
At `t=0` the adjoint interval is truncated.  There is no negative-time parent and
no fresh-interior charge.

For a dyadic initial band `M`, if every selected root has critical mass
`ME_a>=eta`,

\[
\boxed{
\#\mathcal R_M(0)
\le
\eta^{-1}M^{1-2m}\|u_0\|_{\dot H^m}^2.
}
\]

For any fixed geometric band sequence and `m>1/2`, these boundary counts are
summable in `M`.  Smooth initial data therefore supplies a genuine high-frequency
absorbing tail.  This is still conceptually different from a uniform interior
reset: the ancestry **terminates** at the physical boundary.

## 9. Forbidden double-charge matrix

The compiler uses the following canonical pair decisions.

| Pair | Relation | Compiler rule |
|---|---|---|
| H1/H3 dephasing vs causing pressure/SGS/viscous source | downstream | source owns the event; dephasing is diagnostic unless it independently reaches a direct transfer endpoint first |
| high strain vs forced `D_V` lower bound | downstream | one critical-dissipation event |
| physical covariance change vs causing strain/source | downstream | charge the physical cause, not a second covariance currency |
| nearby covariance **representative** change vs physical covariance deformation | primary/diagnostic | small representative update is `Xi_cov`; large physical deformation stays source/strain/relink |
| coherent relinking vs omitted cross-cell `Xi` | mutually exclusive | retained graph vs deleted moat support |
| new coherent Moyal mass vs fresh affine radius | downstream | Moyal new mass is the ancestry event; affine radius certifies its physical size |
| pair rescue vs reuse-cycle/Renyi currency | downstream | rescue is an intermediate route, not a second payment |
| transfer deficit vs backscatter/cancellation | downstream | cancellation is a physical mechanism of the same transfer loss |
| physical cross-cell transfer vs symbol-freezing approximation | independent | distinct error operators may both enter `Xi`, each once |
| Duhamel classified residual vs underlying source/interface term | downstream | residual delegates; it is never a new currency |
| initial-boundary termination vs “fresh packet” | mutually exclusive | `t=0` absorbs ancestry and forbids an interior-fresh relabel |

Two truly independent conserved resources may both be recorded as side ledgers,
but only one copy of any given physical transfer/source measure may enter the
primary telescope.

## 10. The theorem actually proved here

**Causal-quotient transfer-measure partition theorem.**  Suppose an extracted
retained efficient block supplies:

1. the positive physical child-transfer law and the certified cross-cell `Xi`
   excision;
2. the fixed transfer-loss gate;
3. measurable causal-root stopping witnesses, with duplicate theorem
   manifestations identified by physical provenance;
4. a Kelvin-flat certificate when no stopping witness occurs;
5. positive stopping weights on exact independent tie sets.

Then the formulas above define a nonnegative, exhaustive, exactly mass-preserving
single-charge partition of the block's physical transfer measure.  The master
projection is unique except that scale-critical source/mass/dissipation branches
remain marked `recurse`; they cannot be silently inserted into the finite reset
count.

The theorem also identifies the remaining causal-measure bridge precisely:
normalized adjoint `dGamma` may be used as a conditional ancestry kernel, but the
actual PDE extraction must still identify that kernel with the same physical
transfer-weighted parent law used by Shannon/Renyi/Hodge, or route their
discrepancy once.

This is a compiler theorem for an **actual physical block state once its PDE
witnesses exist**.  It is not a proof that every continuum Navier--Stokes block
already supplies those witnesses.  No global-regularity claim follows.
