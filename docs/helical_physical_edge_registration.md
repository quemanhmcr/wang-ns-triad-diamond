# Helical physical edge registration: read the capacity from Navier–Stokes itself

Status: **certified on exact implementation SHA `6d2b721be579e299c42f4ad370448e5336f3f43c`**.

## 1. Do not start from a Young norm

The missing bridge in the signed-good architecture is not another stability
estimate.  The downstream single-edge and polarization theorems already know what
to do with a physical edge carrying

`capacity A_e`, `geometry J_e/J_*`, `phase c_e`.

What was missing was a production theorem saying what those quantities are on an
**actual Navier--Stokes high--high Fourier interaction**.

The capacity must therefore be read from the nonlinear term itself, not chosen as
a convenient norm product.

## 2. One unordered physical parent orbit

Take two helical parent modes at `x,y` and child `z=x+y`,

`u_x=a_x h_(s_x)(x)`, `u_y=a_y h_(s_y)(y)`, `u_z=a_z h_(s_z)(z)`.

For incompressible Navier--Stokes,

`-P(u.grad u)=P(u x omega)`.

The unordered convolution orbit consists of the two physical parent orders

`u_x x omega_y + u_y x omega_x`.

Because `i k x h_s(k)=s|k|h_s(k)`, their sum is

`(s_y|y|-s_x|x|) a_xa_y h_x x h_y`.

Pairing with the divergence-free child helical vector makes the Leray projector
invisible.  With the repository's closed-triad Waleffe convention

`g_e=coupling_g(x,y,-z,s_x,s_y,s_z)`,

reality of the helical frame gives exactly

`conj(h_z).(h_x x h_y)=-2 conj(g_e)`.

Hence

`<h_z,F_z>=2(s_x|x|-s_y|y|)conj(g_e)a_xa_y`.

The implementation computes both the direct curl/Leray vector and this Waleffe
coefficient and requires equality; the Waleffe formula is not used as the
definition of the direct NS source.

## 3. Physical work supplies the factor four

The production energy law uses

`T_e=2 Re[conj(a_z)<h_z,F_z>]`.

Therefore

`T_e=4(s_x|x|-s_y|y|) Re[conj(a_z)conj(g_e)a_xa_y]`.

The factor four has two physical pieces already present elsewhere in the
repository:

- factor two from the two ordered parent terms in one unordered convolution orbit;
- factor two from differentiating child kinetic energy.

No packet count or observer orientation enters.

The same physical edge also obeys the pointwise capacity domination

`|T_e| <= A_e`.

Indeed, in divergence form `u_x.y=u_x.z` and `u_y.x=u_y.z` by incompressibility
and `x+y=z`.  The two ordered parent terms therefore have total source norm at
most `2|z||a_x a_y|`; pairing with the child and differentiating child energy
adds the remaining factor two.  This remains physical on nonforward/backscatter
edges where the upper-progress multiplier is zero.

## 4. Native modal capacity, geometry and phase

Define

`A_e=4|z||a_xa_ya_z|`.

This is the **available modal interaction amplitude** of the same edge.  It is
not the sharp-Young `L^(3/2)` capacity used later to compare a block with a
near-extremizer.

Let `p_top=max(|x|,|y|)` and

`gscale=log_+(|z|/p_top)`.

The existing single-edge geometry is

`J_e=gscale |s_x|x|-s_y|y|| |g_e|/|z|`.

Let

`c_e=sign(s_x|x|-s_y|y|)
 Re[conj(a_z)conj(g_e)a_xa_y] / |g_e a_xa_ya_z|`,

with `c_e=0` on a zero interaction factor.  Then `c_e in [-1,1]` and

`m_e=J_e/J_* in [0,1]`

by the already certified sign-exhausted single-edge envelope.

The exact physical registration identity is

`T_e gscale = A_e J_e c_e`.

Thus the synthetic arrays accepted by the old polarization theorem now have a
natural candidate meaning on one actual Fourier edge.

## 5. The representation quotients are exact

Swapping `x` and `y` changes both

`s_x|x|-s_y|y| -> -(s_x|x|-s_y|y|)`

and

`g_e -> -g_e`.

Their product, the direct source, child work, `A_e`, `J_e` and `c_e` are unchanged.
Parent orientation is therefore not physics.

Likewise reciprocal helical-basis phase changes transform `g_e` and modal
coefficients oppositely, leaving the physical phase alignment and edge identity
invariant.

Uniform wavevector dilation leaves the dimensionless geometry `J_e`, `J_e/J_*`,
`c_e`, and the forward ratio unchanged; both actual work and native capacity scale
linearly, as the derivative order of Navier--Stokes requires.

## 6. What this does not yet prove

This is deliberately a **one-edge theorem**.

It does not yet:

- construct the continuum unordered Fourier/helicity edge measure;
- show that the block transfer deficit used by `coherent_service_or_flat` is the
  deficit of that same continuum physical edge law;
- claim that nonforward or generic comparable HH is signed-good;
- claim a raw majority of positive HH work lies on the signed-good capacity core.

Those are measure-registration questions for the next theorem.  Backscatter and
nonforward positive child work remain physical and are not erased by setting their
upper-progress multiplier to zero.

No Navier--Stokes global-regularity conclusion is asserted.

## 7. Certification

Dedicated GitHub Actions run `31472591138` passed the `748`-test theorem suite,
`50,000` random physical helical triads/amplitudes, the global single-edge
`J_*` dependency, smooth-flux polarization, and physical-pair factor-four/productivity.
The dedicated stress recorded worst direct-Leray/Waleffe relative residual
`3.0472764930177753e-14`, worst `T log=A J c` residual
`1.2645440250480533e-13`, maximum sampled `J/J_*=0.9990427280708386`,
`16543` positive forward-work samples, and `16950` nonforward samples retained
with zero upper-progress multiplier.

The dedicated artifact is stored under
`recorded-results/31472591138/helical-physical-edge-registration-results/` and has
GitHub digest
`sha256:0aa0e684d210dc8af109dc9cfe31b0b312967468fb3a5bd79d4c29a6143fda3a`.

Full physical-energy causal integration run `31472591129` also succeeded on the
same exact SHA through `60` successful steps, including continuum master,
complex Young, dual/Bargmann/common-slice registration, smooth `Q^2`, event roles,
physical pair productivity, recursive witness, signed-good generated epoch,
Shannon/Rényi, branch compiler and master episode. Its artifact digest is
`sha256:79e9c35a1b4de41830483d15e4fc0f2d26305eb29a2b2ba9239b48a9954ea27d`.
