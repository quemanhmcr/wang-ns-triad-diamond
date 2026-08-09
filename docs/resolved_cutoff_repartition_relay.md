# Resolved-cutoff repartition is a decomposition gauge

## Status

**EXACT_RESOLVED_CUTOFF_REPARTITION_GAUGE_AND_PARENT_SCALE_RENEWAL__GENERATED_SURVIVOR_REENTRY_NO_NEW_INTERFACE**

A recursive parent lives at a lower frequency than its child.  It is therefore natural to change the resolved transporter from `S_(N/4)u` to `S_(N_p/4)u` when the parent becomes the new carried scale.  The apparent difficulty is created only if the two cutoffs are regarded as defining different physics.

They do not.

For any resolved field `V`, write `h=u-V` and `L_V f=B(V,f)+B(f,V)`.  For the same smooth role `w=Q u`, the cutoff-dependent nonlinear part of the outer-role equation, after putting the transporter on the right, is

`G_V=-L_V(Q u)+Q B(V,V)-Q B(h,h)+(L_V Q-Q L_V)u`.

The first and fourth terms cancel at the `Q u` level, leaving

`G_V=Q[B(V,V)-B(h,h)-L_Vu]`.

Since

`B(h,h)=B(u,u)-B(u,V)-B(V,u)+B(V,V)`,

one obtains identically

`G_V=-Q B(u,u)`.

Hence `G_V` does not depend on the cutoff at all.  Changing `V` only changes which part of the same quadratic interaction is described as resolved transport, HH source, or Heisenberg interface.  There is no instantaneous cutoff-switch forcing and no new currency.

This identity is especially useful at common-slice relay.  A signed-good parent scale satisfies

`3/5 < N_p/N < 5/8`.

The smooth event envelope has lower edge `11N/20`.  One low-strain transport to the common slice and one renewed low-strain parent slab cost at most `e^(-1/15)` together.  Thus in units of the new parent scale,

`|xi|/N_p >= (11/20)e^(-1/15)/(5/8)=(22/25)e^(-1/15)>1/2`.

So after the cutoff is changed to `V_p=S_(N_p/4)u`, its low-low output remains below `N_p/2` and is excluded by the same relayed carrier for the entire renewed slab.  The lifetime transforms parabolically:

`64/25 < T_p/T < 25/9`,

exactly the already certified asynchronous window.

The lesson is structural: **scale renewal changes the resolution at which we describe the flow, not the nonlinear law being described.**

Together with smooth material-carrier relay, the generated-survivor route now needs neither a fresh hard packet at the common slice nor an artificial cutoff-interface charge.  The remaining work is to prove universal physical slab renewal/exhaustion for every recursive route, especially source/critical-dissipation/material-relink routes.  No global-regularity claim is made.
