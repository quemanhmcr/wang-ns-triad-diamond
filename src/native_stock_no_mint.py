from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Sequence

STATUS='DRAFT_NATIVE_STOCK_NO_MINT__Q2_ENERGY_CONTINUITY_BEFORE_OWNER_LABELS'

@dataclass(frozen=True)
class StockStep:
    energy_in: float
    energy_out: float
    dissipation: float
    signed_physical_work: float
    def __post_init__(self):
        xs=(self.energy_in,self.energy_out,self.dissipation,self.signed_physical_work)
        if not all(math.isfinite(x) for x in xs): raise ValueError('finite stock data required')
        if min(self.energy_in,self.energy_out,self.dissipation)<0: raise ValueError('nonnegative stock/dissipation required')
        r=self.energy_out+self.dissipation-self.energy_in-self.signed_physical_work
        if abs(r)>1e-11*max(1.0,*(abs(x) for x in xs)): raise ValueError('E_out+D=E_in+W_phys required')

def stock_telescope(steps: Sequence[StockStep]):
    rows=tuple(steps)
    if not rows: raise ValueError('nonempty physical stock path required')
    for a,b in zip(rows,rows[1:]):
        if abs(a.energy_out-b.energy_in)>1e-11*max(1.0,a.energy_out,b.energy_in):
            raise ValueError('owner change attempted to reset endpoint stock')
    D=sum(x.dissipation for x in rows); W=sum(x.signed_physical_work for x in rows)
    r=rows[-1].energy_out+D-rows[0].energy_in-W
    return {'steps':len(rows),'initial':rows[0].energy_in,'final':rows[-1].energy_out,'dissipation':D,'signed_work':W,'residual':r,'owner_reset':False}

def theorem_certificate():
    return {
      'status':STATUS,
      'carrier_law':'eta=Q^2: dE_Q/dt+D_Q=<u,dot eta u>-2 Re<eta u,B(u,u)>',
      'quadratic_partition':'sum A_a^2=I: channel energies reconstruct full kinetic energy; differentiated partition motion has zero total energy work',
      'cutoff_invariance':'-L_V(Qu)+QB(V,V)-QB(u-V,u-V)+(L_VQ-QL_V)u=-QB(u,u), independent of V',
      'relink':'after common transported-gauge quotient, K_phys pair work is antisymmetric and has zero total work',
      'no_mint':'source/strain/HH may change event ontology but cannot reset physical stock ancestry; fresh stock requires inherited energy or actual signed physical work',
      'scope':'does not rule out conservative circulation and makes no global-regularity claim',
    }
