#!/usr/bin/env python3
"""Table design study: Korf's formula on sampled lookups (examples/tabsample.rs records).

relative nodes = sum over k of 12.93^-k P(bound <= k), the bound the maximum of a design's
lookups on uniformly random positions. usage: scripts/korf.py samples.bin (15-byte records)
"""
import numpy as np, sys
B=12.93
names=['AE','AEsw','A0','A1','A2','A3','B0','B1','B2','B3','CE','CCB','CCA','EA','EAsw']
ix={n:i for i,n in enumerate(names)}
CBE=['A0','A1','A2','A3','B0','B1','B2','B3']
def korf(h, kmax=18):
    n=len(h); c=np.bincount(h,minlength=kmax+1)[:kmax+1].cumsum()/n
    terms=np.array([B**(10-k)*c[k] for k in range(kmax+1)])
    return terms.sum(), terms, c
designs=[
 ('medium',['CE','CCB','CCA','AE','AEsw'],3.4),
 ('large',['AE','AEsw']+CBE,22.0),
 ('EA/c x2',['EA','EAsw'],27.9),
 ('EA/c A-side only',['EA'],27.9),
 ('EA/c x2 + CE',['EA','EAsw','CE'],30.3),
 ('EA/c x2 + CBE x8',['EA','EAsw']+CBE,48.9),
 ('EA/c x2 + CBE A0,B0',['EA','EAsw','A0','B0'],48.9),
 ('EA/c x2 + CBE x8 + CE',['EA','EAsw','CE']+CBE,51.3),
 ('EA/c x2 + AE x2 (check: AE dominated)',['EA','EAsw','AE','AEsw'],28.9),
 ('large + CE',['AE','AEsw','CE']+CBE,24.4),
]
a=np.fromfile(sys.argv[1],dtype=np.uint8).reshape(-1,15)
print(len(a),'samples')
for n in ['EA','EAsw']:
    h=a[:,ix[n]]; print(n,'mean',h.mean().round(3),np.bincount(h,minlength=15))
print('EA >= AE always:', (a[:,ix['EA']]>=a[:,ix['AE']]).all(), (a[:,ix['EAsw']]>=a[:,ix['AEsw']]).all())
res={}
for name,cols,mem in designs:
    h=a[:,[ix[c] for c in cols]].max(axis=1)
    s,t,c=korf(h); res[name]=s
for name,cols,mem in designs:
    h=a[:,[ix[c] for c in cols]].max(axis=1)
    s,t,c=korf(h)
    print(f'{name:34s} {mem:5.1f}GB mean {h.mean():6.3f} korf {s:.5f} vs medium {s/res["medium"]:.4f} vs large {s/res["large"]:.4f} | P<=10 {c[10]:.2e} P<=11 {c[11]:.4f} P<=12 {c[12]:.4f} P<=13 {c[13]:.4f}')
    print('      terms k=9..14:', ' '.join(f'{x:.5f}' for x in t[9:15]))
