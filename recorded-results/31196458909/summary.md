# Physical mild-aspect H1/swirl local no-escape theorem

Status: **CERTIFIED_PHYSICAL_LOW_STRAIN_CONDITIONING** + **EXACT_AFTER_PHYSICAL_ROLE_SPLIT**.

Two corrections are built into this final physical theorem. First, `Q_pol/2` is an auxiliary relative-coordinate forcing norm; the three physical Young roles only satisfy `sum_i||F_i^H1||^2 >= Q_pol/4 >= ||B_hook||^2/100`. Second, the base polarization propagators are non-unitary. On the existing low-strain lifetime branch each physical role has action budget `K<=1/30`, so pullback and pushforward singular values are controlled by `exp(+-K)` rather than treated as isometries.

Let `I1=int||B_hook||dt`. The conditioned interaction-picture variation theorem gives

`J1 >= I1/(11 T)`

or a **physical three-role** first-Duhamel daughter with

`delta1^2 >= I1^2/480`.

If nonlinear physical feedback is less than half, total surviving daughter energy is at least `I1^2/1920`; one of the three roles therefore has energy at least `I1^2/5760`. Below critical sideband size `1/80`, odd-Hermite convexity and the pair-rescue split yield

`net transfer deficit >= I1^2/184320`

or

`pair-sideband rescue >= I1^2/184320`.

Combining with the H3 branch and `I_B<=sqrt(6)I3+I1` gives the clean physical mild-aspect full-curvature cost

`pair rescue or transfer deficit >= I_B^2/737280`

outside dephasing/source, nonlinear-feedback and large-daughter branches.

This supersedes the idealized relative-coordinate/isometric-pullback constants `1/25600` and `1/102400` from run `31195130386`. The pointwise mild-aspect bridge `Q_pol>=1/25||B_hook||^2` is unchanged. H1 covariant dephasing source calculus remains open, and high-aspect grains remain ancestry/reuse rather than an aspect defect.

Stress: `50000`
- branch counts: `{'H1_covariant_dephasing': 9978, 'nonlinear_sideband_feedback': 10144, 'pair_sideband_rescue': 10006, 'transfer_deficit': 9926, 'large_daughter_capacity': 9946}`
- minimum H1 pair-cost margin: `1.546e-11`
- minimum H1 deficit-cost margin: `5.585e-12`
- minimum full-channel margin: `5.238e-09`
