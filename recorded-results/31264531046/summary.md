# Coherent averaged strain source — fixture failure provenance

Run `31264531046` on `3d508a3` failed in pytest before theorem stress execution.

Classification: **fixture / floating equality semantics**, not a mathematical countermodel.

The failing assertion compared `190*0.01`, represented as `1.9000000000000001`, with the literal `1.9` using exact `==`. The same clean pressure/SGS route constants were already theorem-level in the resolved objective-strain module. The fixture was corrected to use `math.isclose`; no theorem formula or constant was changed.
