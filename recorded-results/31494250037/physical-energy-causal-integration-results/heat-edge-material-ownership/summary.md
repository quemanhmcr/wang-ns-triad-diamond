# Material ownership of positive heat-increment edges

Status: **EXACT_POSITIVE_HEAT_EDGE_MATERIAL_OWNERSHIP_PARTITION__AFFINE_ENDPOINT_LABELS_INVARIANT__CAPACITY_ROUTING_REMAINS**.

The high-strain heat theorem already supplies a positive coherent edge measure.  Material ownership should be assigned to that **measure**, not to a decomposition of the velocity field.

For one coherent heat edge, translation covariance identifies its endpoints as `(X,k)` and `(X-r,k)`.  In the intrinsic material coordinate

`zeta=(L^-1 X/2,L^T k)`,

let `O` be the fixed transported old material pool and put `chi_i=1_O(zeta_i)`.  Pointwise on every nonnegative edge weight `s`,

`s_OO=chi_0 chi_1 s`,

`s_ON=[chi_0(1-chi_1)+(1-chi_0)chi_1]s`,

`s_NN=(1-chi_0)(1-chi_1)s`,

and exactly

`s_OO+s_ON+s_NN=s`.

There is no quadratic cross term because no identity `V=V_old+V_new` is used.  The partition is performed **after** the Moyal/heat edge density has become positive.  It is also unoriented: swapping the two physical endpoints leaves OO, interface and NN weights unchanged.

Under common affine/Kelvin transport

`L->ML,  X->MX,  r->Mr,  k->M^-T k`,

both endpoint labels are individually invariant:

`zeta(ML,MX,M^-T k)=zeta(L,X,k)`,

`zeta(ML,M(X-r),M^-T k)=zeta(L,X-r,k)`.

Thus common affine motion cannot convert an old--old edge into an interface or new--new edge.  A change of endpoint membership along a continuous nonaffine material evolution must meet the boundary of the old material set.  Half-open dyadic conventions merely assign that null boundary and cannot manufacture service mass.

The exact coherent increment identity also yields, class by class,

`|e^(-ik.r)A_1-A_0|^2 <= 2(|A_0|^2+|A_1|^2)`.

Hence each ownership class is locally supported by the Moyal energy of its own two endpoints.  This is the correct starting point for old-pool capacity; the theorem does **not** import the signed-good low-strain half-life through a high-strain slab.

Stress: `50000` material-edge states
- worst affine two-endpoint residual: `2.681e-15`
- worst positive-measure partition residual: `1.705e-13`
- minimum ownership-specific endpoint-capacity margin: `0.000e+00`
- orientation failures: `0`
- affine membership-invariance failures: `0`

This closes **material ownership classification** of the high-strain heat-edge seed.  The remaining high-strain problem is quantitative routing: bound repeated OO heat service using the physically valid history of the transported old pool; route ON service as actual material-interface/relink provenance; and show what NN service creates without postulating a packet mass floor.  Universal slab renewal remains open, and no global-regularity claim is made.
