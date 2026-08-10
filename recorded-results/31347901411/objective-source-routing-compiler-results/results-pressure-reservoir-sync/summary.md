# Pressure reservoir synchronization: derivative order sets the pair lifetime

Status: **EXACT_OBJECTIVE_PRESSURE_HESSIAN_AND_PRESSURE_THIRD_PAIR_EROSION_GIVEN_SIGNED_GOOD_LOW_STRAIN_LINEAGE**.

There are two distinct pressure source objects and they must not be conflated.

For the **resolved objective pressure Hessian** `N^-4||grad^2 P||_inf`, one low--low pair `(a,b)` has, up to the fixed Riesz/Bernstein constant, per-unit `sqrt(E_a E_b)` coefficient

`M_max^4 sqrt(M_a M_b) / N^4`.

On a signed-good low-strain lineage, reservoir frequencies grow by less than `21/20` while the block scale grows by more than `8/5`.  Hence the objective-Hessian pair coefficient contracts by

`(21/20)^5 (5/8)^4 = 4084101/20971520 < 1/5`.

Even allowing both reservoirs the entire global energy cap at every future generation, one fixed objective-Hessian pair has total future capacity less than `5/4` times its generation-zero capacity.

For the separate **H1 pressure-third** source, the coefficient is

`M_max^3 sqrt(M_a M_b) / N^3`,

and its previously certified contraction remains

`(21/20)^4 (5/8)^3 = 194481/655360 < 1/3`,

with total future fixed-pair capacity `<3/2` of generation zero.

Thus the derivative order is physical provenance: objective Hessian and pressure-third share the same material-pair erosion mechanism but have different exact lifetimes.  Persistent objective pressure service must relink pairs, fragment into pair/component entropy, leave the supplied signed-good low-strain lineage, or use its SGS-stress alternative.

Stress: `50000`
- maximum sampled objective-Hessian pair ratio: `0.188320537`
- minimum margin below `1/5`: `1.168e-02`
- minimum objective-Hessian one-fifth envelope margin: `0.000e+00`
- maximum sampled pressure-third pair ratio: `0.288962935`
- minimum margin below `1/3`: `4.437e-02`
- minimum pressure-third one-third envelope margin: `0.000e+00`
