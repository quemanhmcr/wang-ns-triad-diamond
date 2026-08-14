# Intrinsic Navier–Stokes law set and mixed genuine-owner reduction

Status: **DRAFT THEOREM CANDIDATE — the local laws are native; the final master-depth composition statement remains conditional on one explicit wiring theorem.  This PR is not merge-ready.**

The important reduction is now **below** the earlier ontology

\[
\text{STOCK}\;|\;\text{SIGNED POWER}\;|\;\text{LOCAL ACTION}.
\]

Those three types are useful recurrence manifestations, but they are not the smallest mechanisms.  The present candidate picture is that a much smaller internal law set generates them.

## I. Closed-triad current — the primitive beneath nonlinear POWER

For one regular closed helical triad, put

\[
a_i=s_i|k_i|.
\]

Its three cyclic modal energy works are not independent.  They are one scalar current \(R\) seen through the triad geometry:

\[
\boxed{
T_0=(a_1-a_2)R,\qquad
T_1=(a_2-a_0)R,\qquad
T_2=(a_0-a_1)R.
}
\]

Consequently the same one-dimensional current has two simultaneous null laws,

\[
\boxed{T_0+T_1+T_2=0},
\qquad
\boxed{\sum_i s_i|k_i|T_i=0}.
\]

Thus energy conservation and signed-helicity conservation are not two unrelated recurrence mechanisms.  They are two constraints on the same physical triad current.  The phase/amplitude description can be complicated while the modal energy routing of the actual closed triad remains one-dimensional.

### No perfect free forward transfer

For an already-certified signed/geometry-good positive forward child the parent/child ratios satisfy

\[
\frac35<D,S<\frac58,
\]

with opposite parent helicities.  The parent sharing the child's helicity is the unique energy donor; the other cyclic root is a simultaneous positive nonforward recipient with \(J=0\), hence existing `TRANSFER_WORK_LOSS`.

Normalize the child frequency to one and write

\[
C=W_{\rm child}^+,\qquad
D_w=W_{\rm donor}^-,\qquad
B=W_{\rm side}^+.
\]

The one-current law gives

\[
\boxed{
C=(D+S)R,\qquad
D_w=(1+S)R,\qquad
B=(1-D)R,
}
\]

so in particular

\[
\boxed{D_w=C+B},
\]

and the clean geometry window yields

\[
\boxed{\frac34<\frac{C}{D_w}<\frac{10}{13}},
\qquad
\boxed{\frac3{13}<\frac{B}{D_w}<\frac14}.
\]

A good high child therefore never receives all of the donor work.  Nature itself creates a same-event nonforward sibling.  This is not entropy, packet splitting, a hard-cell convention, or an analyst-imposed branch.

### Physical restriction does not remove the branching

Let \(q\ge0\) be one common measurable weight applied to the **lifted closed-triad occurrence before coarsening**.  For any measurable family of good occurrences set

\[
G=\int qC,\qquad
D_*=\int qD_w,\qquad
B_*=\int qB.
\]

Then atomwise integration preserves

\[
\boxed{D_*=G+B_*},
\qquad
\boxed{G\le\frac{10}{13}D_*},
\qquad
\boxed{B_*\ge\frac3{13}D_*}.
\]

This is compatible with the existing downstream rules whenever the same physical positive restriction is inherited: time restriction, hard-cell restriction, resolved-contact weights, donor restrictions, and analogous routes that inherit canonical \(dW^+\) without a new Hahn split.

Hence the native same-event good-continuation factor is

\[
\mathcal R_{\rm cyc}=\frac{G}{D_*}\le\frac{10}{13},
\]

or

\[
\boxed{
C_{\rm cyc}=-\log \mathcal R_{\rm cyc}
\ge \log\frac{13}{10}
=0.262364264\ldots .
}
\]

The side work is real energy and may continue evolving later.  The theorem does **not** say that the good child terminates.  It says that selected free good-child continuation cannot retain more than \(10/13\) of the actual same-event donor-row work.  If cyclic continuation cost is used, the same side mass is not charged again by Young/entropy/hard-cell bookkeeping.

The mixed-fate hard-cell and Young/Christ seams are therefore not load-bearing for this POWER law: the branching is already present before those coarsenings.

A non-load-bearing algebra audit over 100,000 randomly weighted signed-good geometry families (seed `2026081401`) gave

- maximum retained fraction `0.769211723826083 < 10/13`;
- minimum side fraction `0.230788276173917 > 3/13`;
- minimum local log cost `0.262389023800093 > log(13/10)`.

The audit is only a regression check; the theorem is the atomwise current identity plus common positive restriction.

## II. Mode-stock continuity — the primitive beneath STOCK

For a physical Fourier/helical mode set \(A\), the exact continuity law has the form

\[
\boxed{
E_A(t_1)+D_A+\Phi_{\rm out}
=
E_A(t_0)+\Phi_{\rm in}.
}
\]

Equivalently for one stock segment,

\[
\boxed{E_{\rm out}+D=E_{\rm in}+W_{\rm phys}.}
\]

Nonlinearity transfers stock, viscosity is the sink, and representation/owner changes introduce no reset term.  Material rereads, same-carrier inheritance, cutoff repartition and conservative relink can change the description of the state but cannot mint generation depth.  Fresh stock requires inherited physical energy or actual signed physical work.

This law is why same-time donor provenance must not be promoted to a between-time wallet.

## III. Local-action speed lock — the primitive beneath ACTION

The native high-strain face is reached when

\[
K_N=\int\|S_V\|\,dt=\frac1{30}.
\]

On every compact pre-singular interval \(I\), smoothness gives the scale-free rate bound

\[
\dot K_N\le g_1\|\nabla u\|_{L^\infty(I)},
\]

hence a positive physical-time price independent of \(N\).

For the objective-source action, the same physical observable obeys simultaneously

\[
\dot A_{{\rm obj},N}\le aN^3+bN^{5/2}
\]

from the resolved source calculus and

\[
\dot A_{{\rm obj},N}\le dN^{-2}
\]

from compact pre-singular smoothness.  Therefore

\[
\boxed{
\sup_N \dot A_{{\rm obj},N}
\le
(2a)^{2/5}d^{3/5}+(2b)^{4/9}d^{5/9}<\infty.
}
\]

The fixed native face \(A_{\rm obj}=\tau/60\) consequently also has a positive physical-time price.  Thus

\[
\boxed{N_{\rm ACTION}(I)<\infty}
\]

on every compact pre-singular interval.  No analyst scale cutoff, fixed event-gap axiom, or critical-energy reset budget is introduced.

## The ontology is now a manifestation layer

The earlier classification

\[
\boxed{\text{STOCK}\quad|\quad\text{SIGNED POWER}\quad|\quad\text{LOCAL ACTION}}
\]

remains useful, but it should now be read as the recurrence ontology generated by the smaller laws above.  Labels such as

\[
\text{material/source/strain/HH/contact/relink/fresh}
\]

are still higher-level manifestations or event readers, not additional primitive mechanisms.

In this reading:

- STOCK does not mint recursive depth;
- ACTION-containing vertices are finite on compact pre-singular intervals;
- geometry-bad/nonforward POWER is already the existing transfer-loss route;
- geometry-good POWER pays the intrinsic same-event cyclic branching cost.

## Temporal-ancestry guard

There is an exact Eulerian stock/flow generator identity, but it must not be converted into a temporal matching rule for individual energy quanta.  In particular, a Radon–Nikodym quotient of same-time donor work against modal energy must **not** be interpreted as a canonical Markov ancestry deciding which earlier deposit pays which later withdrawal.

No FIFO, no LIFO, no proportional temporal mixing, and no hidden old-stock-to-later-work matching enters the cyclic theorem.  The \(10/13\) law is entirely same-time.

Signed helicity is used here as a structural null law of the same closed-triad current.  Absolute-helicity magnitude is **not** promoted to a finite transfer budget or charged as an additional cost.

## Candidate master-facing consequence

After at most \(N_{\rm ACTION}(I)\) ACTION-containing vertices, a hypothetical infinite free survivor must be a POWER continuation.  Bad POWER has already left through transfer loss.  Every remaining good POWER junction has the local ceiling \(10/13\).

Let \(N_{\rm ACTION}(I)\) denote an **upper bound** on the number of ACTION-containing vertices in \(I\).  A depth-\(L\) survivor therefore contains at least \(L-N_{\rm ACTION}(I)\) free POWER vertices once \(L\ge N_{\rm ACTION}(I)\).  Thus, **conditional on the remaining master composition theorem**, the natural master-facing estimate is

\[
\boxed{
\mathcal R_L
\lesssim
\left(\frac{10}{13}\right)^{L-N_{\rm ACTION}(I)}
\qquad (L\ge N_{\rm ACTION}(I)).
}
\]

For all \(L\), the exponent is read as the positive part \([L-N_{\rm ACTION}(I)]_+\).

Equivalently, after separating the already-existing finite/summable prefactor,

\[
\boxed{
-\log \mathcal R_L
\ge
\log\frac{13}{10}\,[L-N_{\rm ACTION}(I)]_+
}
\]

up to those finite/summable terms.

Exact ACTION/POWER ties remain **joint physical causes**.  For the cost count, an ACTION-containing tie may be included among the finite ACTION exceptions and its simultaneous POWER co-charge omitted.  This preserves the joint cause while preventing double charging.

The symbol \(\lesssim\) above does not hide a new physical budget.  It only allows the finite/summable representation terms already present in the master architecture.

### The one remaining load-bearing seam

The local statement

\[
G/D_*\le10/13
\]

is a certified same-event physical law on the closed-triad Radon level.  What is **not yet certified** is that the master's recursive continuation variable is exactly the object on which these local factors multiply across arbitrary depth.

That wiring theorem must prove composition without:

- treating \(dW^-\) as between-time stock;
- performing temporal matching;
- replacing the canonical \(dW^+\) cause;
- re-Hahn splitting after restriction;
- charging the side sibling twice.

If that theorem passes review, the mixed genuine-owner pillar no longer has an independent pure-POWER frontier: the apparently complicated owner recurrence is generated by the three smaller internal laws above, and free good POWER pays a branching price imposed by energy and signed-helicity conservation themselves.

No global-regularity claim is made here.
