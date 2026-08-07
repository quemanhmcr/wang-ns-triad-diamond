from __future__ import annotations
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
from scipy.stats import qmc


def area(x,y,z):
    t=(x+y+z)*(-x+y+z)*(x-y+z)*(x+y-z)
    return .25*np.sqrt(np.maximum(t,0.0))

def eff(x,y,z,sx,sy,sz):
    valid=(z>np.maximum(x,y)*(1+1e-8)) & (x+y>z) & (np.abs(x-y)<z)
    g=area(x,y,z)*np.abs(sx*x+sy*y+sz*z)/(2*np.sqrt(2)*x*y*z)
    out=np.log(np.maximum(z/np.maximum(x,y),1.0))*np.abs(sx*x-sy*y)/z*g
    return np.where(valid,out,0.0)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--power',type=int,default=18); ap.add_argument('--outdir',default='results-fast')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    # Global single-edge probe in (x,y), child normalized to 1.
    uv=qmc.Sobol(2,scramble=True,seed=111).random_base2(a.power)
    x=.05+.949*uv[:,0]; y=.05+.949*uv[:,1]; lo=np.minimum(x,y); hi=np.maximum(x,y)
    single={'J':-1.0}
    for signs in itertools.product((-1,1),repeat=3):
        v=eff(lo,hi,np.ones_like(lo),*signs); i=int(np.argmax(v))
        if v[i]>single['J']: single={'J':float(v[i]),'x':float(lo[i]),'y':float(hi[i]),'signs':list(signs)}
    # Diamond geometries.
    u=qmc.Sobol(5,scramble=True,seed=222).random_base2(a.power)
    rb=np.exp(math.log(.25)+(math.log(4)-math.log(.25))*u[:,0]); rc=np.exp(math.log(.25)+(math.log(4)-math.log(.25))*u[:,1])
    th=.03+(math.pi-.06)*u[:,2]; ph=.03+(math.pi-.06)*u[:,3]; ps=2*math.pi*u[:,4]
    A=np.column_stack([np.ones(len(u)),np.zeros(len(u)),np.zeros(len(u))])
    B=rb[:,None]*np.column_stack([np.cos(th),np.sin(th),np.zeros(len(u))])
    C=rc[:,None]*np.column_stack([np.sin(ph)*np.cos(ps),np.sin(ph)*np.sin(ps),np.cos(ph)])
    M=A+B; N=B+C; D=A+B+C
    L={k:np.linalg.norm(v,axis=1) for k,v in {'a':A,'b':B,'c':C,'m':M,'n':N,'d':D}.items()}
    best3={'ratio':-1}; best4={'ratio':-1}
    for signs in itertools.product((-1,1),repeat=6):
        sa,sb,sc,sm,sn,sd=signs
        e1=eff(L['a'],L['b'],L['m'],sa,sb,sm)
        e2=eff(L['m'],L['c'],L['d'],sm,sc,sd)
        e3=eff(L['b'],L['c'],L['n'],sb,sc,sn)
        e4=eff(L['a'],L['n'],L['d'],sa,sn,sd)
        r3=np.minimum(np.minimum(e1,e2),e3)/single['J']; i3=int(np.argmax(r3))
        if r3[i3]>best3['ratio']:
            best3={'ratio':float(r3[i3]),'signs':list(signs),'index':i3,'edge_ratios':[float(e1[i3]/single['J']),float(e2[i3]/single['J']),float(e3[i3]/single['J'])], 'lengths':{k:float(v[i3]) for k,v in L.items()}}
        r4=np.minimum(r3,e4/single['J']); i4=int(np.argmax(r4))
        if r4[i4]>best4['ratio']:
            best4={'ratio':float(r4[i4]),'signs':list(signs),'index':i4,'edge_ratios':[float(e1[i4]/single['J']),float(e2[i4]/single['J']),float(e3[i4]/single['J']),float(e4[i4]/single['J'])], 'lengths':{k:float(v[i4]) for k,v in L.items()}}
    payload={'samples':2**a.power,'single':single,'best3':best3,'best4':best4}
    (out/'fast_probe.json').write_text(json.dumps(payload,indent=2))
    md=f"# Fast Sobol probe\n\nSamples: `{2**a.power}`\n\nSingle edge: `{single}`\n\nThree-edge best observed ratio: `{best3['ratio']:.9f}`\n\nFour-edge best observed ratio: `{best4['ratio']:.9f}`\n\nMonte Carlo/Sobol values are lower bounds on the unknown optima, not certified upper bounds.\n"
    (out/'summary.md').write_text(md); print(md); print(json.dumps(payload,indent=2))
if __name__=='__main__': main()
