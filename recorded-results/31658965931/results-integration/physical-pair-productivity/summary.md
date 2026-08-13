# Physical pair-weighted amplitude productivity

Status: **EXACT_PHYSICAL_TRANSFER_WEIGHTED_LOG_PRODUCTIVITY__NO_DUHAMEL_PAIR_WEIGHT_IDENTIFICATION**.

The parent-product law should be averaged under the same positive child-energy work that drives causal Shannon/Renyi.  No Duhamel pair reweighting is necessary.

Let `r_e(t) dt` be actual positive work in one hard parent-pair cell.  Sharp Young gives

`r_e(t) <= C_Y N a_c(t) a_1,e(t) a_2,e(t)`.

Normalize the physical measure by its total `W` and compare it with uniform normalized physical time times the `M` hard pair cells.  Nonnegativity of relative entropy gives

`E_(dT/W) log(a_1 a_2) >= log[W/(T M C_Y N)] - E_(dT/W) log a_c`.

On the generated low-strain branch, `E0,W_R<=E1/5`, `W>=8E1/15`, and energy Gronwall imply `sup E <= C_E W` with `C_E=1.87064343506<2`.  Finite shell volume gives `a_c<=C_Omega sqrt(N E)`.  Since `T=cN^-2` and every L2-normalized terminal coefficient satisfies `alpha_c<=sqrt(N E1)`, all scale powers cancel:

`E_(dT/W) log(a_1 a_2) >= log alpha_c + log Lambda_role,M`.

Dual-Gaussian event marking and the conservative common-slice `1/4` factor on each parent give

`E_(dT/W) log(alpha_p1 alpha_p2) >= log alpha_c + log Lambda_M`,

with default `Lambda_1=0.00185770938993` and `Lambda_M=Lambda_1/M`.

The `1/M` factor is harmless under the actual symbol-freezing refinement.  If `M_j<=M0(j+3)^p`, binary log recursion weights its depth-j constant by `2^(-j-1)`, and

`sum_j 2^(-j-1) log M_j <= log M0 + p log 6 < infinity`.

Thus polynomially refining physical pair cells changes only a finite amplitude-entropy offset; it cannot change the linear reuse slope.  This is precisely the homogeneity expected from a quadratic cascade.

Stress: `50000` random physical-work/time/pair laws and variable-depth recursions
- minimum KL margin: `1.161e-04`
- minimum physical log-product margin: `1.665e+00`
- minimum variable-recursion identity margin: `-1.066e-14`
- maximum sampled tempered-penalty/upper ratio: `0.999945`

The remaining continuum assembly is now local: on the retained signed-good hard event cells, verify the sharp-Young work-density bound with the same physical normalization already used by the SGS transfer theorem, and stop at the existing transfer/phase/relink/source branches whenever complex Young or common-slice registration fails.  Duhamel remains a support/adjoint identity but no Duhamel-to-physical **pair-weight** theorem is required.  No global-regularity claim is made.
