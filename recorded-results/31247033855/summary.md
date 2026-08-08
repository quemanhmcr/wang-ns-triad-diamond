# Physical flat episode: explicit master perturbation and erosion rate

Status: **CERTIFIED_PHYSICAL_BARYCENTRIC_PERTURBATION**.

For a signed-good transfer edge, the exact half-angle formula is

`c_e^2 = c_*^2 exp(2v) - sinh(u/2)^2`.

Arb on the full local box certifies

`|c_e-c_*| <= (6/5)|v|+(1/5)u^2`,

while the exact unequal-parent midpoint decomposition gives

`||n_child-m_e|| <= |u|/2`.

For a normalized positive-transfer coupling with physical Hodge energy

`H=E[2v^2+u^2/2]`

and parent-barycenter mismatch `Delta_b=|b_2-b_1|`, these identities imply

`||b_child-b_1/c_*|| <= e`,

`e <= 2 sqrt(H)+(1/2)H+(5/8)Delta_b`.

On the concentrated master branch `|b_1|>=c_*`, the barycentric potential obeys

`P_child <= P_parent-kappa_* -log(1-e)`.

For `e<=1/2`, one may take the clean perturbation

`zeta <=4 sqrt(H)+H+(5/4)Delta_b`.

Combining with the service-or-flat theorem at `tau=1/100`, `sqrt(H)<=tau/3`, and the synchronized-marginal threshold `Delta_b<=tau`, gives

`zeta <=31 tau/12+tau^2/9`,

hence

`kappa_0=kappa_*-zeta >17/100`.

Thus a physical `1%`-Kelvin-flat, parent-synchronized block is already a quantitative master flat step with a uniform barycentric erosion rate; the abstract per-step `zeta_j` is no longer free.

Stress: `50000`
- worst actual barycenter error / clean bound: `0.956538235`
- minimum potential inequality margin: `9.737e-02`
- margin above clean `kappa_0=0.17`: `4.487e-03`
