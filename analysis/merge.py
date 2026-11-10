import re,sys
# merged-variant saving estimate for iterations where s = best-1-D == 1 (critical)
txt=open(sys.argv[1]).read().split('\n')
it=None; rows={}
def flush():
    if it is None or not rows: return
    D,best=it
    s=best-1-D
    tot=sum(r[0] for r in rows.values())
    if tot<1e8: return
    save=0
    for b,(e,sv,b0,b1,b2) in rows.items():
        p=rows.get(b+1)
        if not p or p[1]==0: continue
        both = (p[3]+p[4]) if s<=1 else (p[4] if s==2 else 0)
        save+= e*both/p[1]/2
    print('iter %d best %d s=%d evals %.3g merged-saving %.1f%%'%(D,best,s,tot,100*save/tot))
for line in txt:
    m=re.match(r'stats: iteration (\d+) \(best (\d+)\)',line)
    if m:
        flush(); it=(int(m.group(1)),int(m.group(2))); rows={}; continue
    m=re.match(r'stats: +(\d+) +(\d+) +(\d+) \(.*\| (\d+) (\d+) (\d+)$',line)
    if m:
        b,e,sv,b0,b1,b2=map(int,m.groups()); rows[b]=(e,sv,b0,b1,b2)
flush()
