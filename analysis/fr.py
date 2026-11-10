import re,sys
# parse stats output: per iteration, evals/survivors by budget, fr_surv
txt=open(sys.argv[1]).read().split('\n')
it=None; E={}; S={}; F=None; res=[]
def flush():
    if it is None or not E: return
    tot=sum(E.values())
    if tot < float(sys.argv[2]) if len(sys.argv)>2 else 1e7: return
    T={}
    for b in S:
        below=sum(v for bb,v in E.items() if bb<b)
        T[b]=below/S[b] if S[b] else 0
    parts={b:(F[b]*T.get(b,0) if F else 0) for b in S}
    lo=max(parts.values()) if parts else 0; hi=sum(parts.values())
    res.append((it,tot,lo/tot,hi/tot,{b:F[b] for b in sorted(S) if F and F[b]}))
for line in txt:
    m=re.match(r'stats: iteration (\d+) \(best (\d+)\)',line)
    if m:
        flush(); it=(int(m.group(1)),int(m.group(2))); E={}; S={}; F=None; continue
    m=re.match(r'stats: +(\d+) +(\d+) +(\d+) \(',line)
    if m:
        b,e,s=map(int,m.groups()); E[b]=e; S[b]=s; continue
    m=re.match(r'stats: forced-route CB would prune survivors by budget: \[(.*)\]',line)
    if m: F=[int(x) for x in m.group(1).split(',')]
flush()
for r in res: print('iter %s evals %.3g  FR saving est %.1f%%..%.1f%%  fr_surv %s'%(r[0],r[1],100*r[2],100*r[3],r[4]))
