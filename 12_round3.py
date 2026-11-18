"""Step 12: third review round: immune-mediated hepatitis phenotype vs docetaxel (exploratory), adagrasib+ICI case details, crude NSCLC ratio."""
import pandas as pd, numpy as np, os, math, glob, json
from math import comb
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R=os.path.join(B,"results"); W=os.path.expanduser("~/fw"); RS=os.path.join(W,"rescan"); DC=os.path.join(W,"doce")
IMH={"immune-mediated hepatitis","immune-mediated hepatic disorder","autoimmune hepatitis","immune-mediated cholangitis"}
c=pd.read_csv(os.path.join(R,"cohort_reports_v2_flags.csv"),dtype={"primaryid":str},keep_default_na=False)
for k in ["ici_exposed","hepato"]: c[k]=c[k].astype(str)=="True"
kept=set(open(os.path.join(W,"kept_primaryids.txt")).read().split())
ici_any={}; 
for fn in glob.glob(os.path.join(RS,"*_ici.tsv")):
    for line in open(fn):
        p,k,role=line.rstrip("\n").split("\t")
        if p in kept: ici_any.setdefault(p,[]).append(f"{k}:{role}")
doce_ps=set()
for fn in glob.glob(os.path.join(RS,"*_doce.tsv")):
    for line in open(fn):
        p,role=line.rstrip("\n").split("\t")
        if p in kept and role=="PS": doce_ps.add(p)
ind={}
for fn in glob.glob(os.path.join(DC,"*_indi.tsv")):
    for line in open(fn):
        p,pt=line.rstrip("\n").split("\t")
        if p in doce_ps: ind.setdefault(p,set()).add(pt.lower())
def nsclc(s):
    il="; ".join(sorted(s))
    if "small cell lung" in il and "non-small" not in il and "non small" not in il: return False
    if any(k in il for k in ["colorectal","colon","rectal","rectum"]): return False
    return any(k in il for k in ["lung","nsclc","bronch"])
want=set(c.primaryid)|doce_ps; pts={}
with open(os.path.join(W,"reac_kept.tsv")) as f:
    for line in f:
        p,_,pt=line.rstrip("\n").partition("\t")
        if p in want: pts.setdefault(p,set()).add(pt)
imh=lambda p: any(x.lower() in IMH for x in pts.get(p,[]))
def fisher(a,b,c_,d):
    n1=a+b; n2=c_+d; k=a+c_; N=n1+n2
    def pr(x): return comb(n1,x)*comb(n2,k-x)/comb(N,k)
    p0=pr(a); return min(1.0,sum(pr(x) for x in range(max(0,k-n2),min(k,n1)+1) if pr(x)<=p0*(1+1e-7)))
def orci(a,b,c_,d):
    if 0 in (a,b,c_,d): a,b,c_,d=a+.5,b+.5,c_+.5,d+.5
    r=(a/b)/(c_/d); se=math.sqrt(1/a+1/b+1/c_+1/d); return [round(r,2),round(r*math.exp(-1.96*se),2),round(r*math.exp(1.96*se),2)]
out={}
c["imh"]=c.primaryid.map(imh)
rows=[]
for lab,df in [("sotorasib all",c[c.target_drug=="sotorasib"]),("adagrasib all",c[c.target_drug=="adagrasib"]),("sotorasib NSCLC",c[(c.target_drug=="sotorasib")&(c.tumor=="NSCLC")]),("adagrasib NSCLC",c[(c.target_drug=="adagrasib")&(c.tumor=="NSCLC")])]:
    e=df[df.ici_exposed]; u=df[~df.ici_exposed]; a,b_,c_,d=int(e.imh.sum()),int((~e.imh).sum()),int(u.imh.sum()),int((~u.imh).sum())
    rows.append(dict(group=lab,ici=f"{a}/{len(e)}",noici=f"{c_}/{len(u)}",OR=orci(a,b_,c_,d),fisher_p=round(fisher(a,b_,c_,d),4)))
dn=[p for p in doce_ps if nsclc(ind.get(p,set()))]; da=list(doce_ps)
for lab,ps in [("docetaxel all",da),("docetaxel NSCLC",dn)]:
    e=[p for p in ps if p in ici_any]; u=[p for p in ps if p not in ici_any]; a=sum(imh(p) for p in e); c_=sum(imh(p) for p in u)
    rows.append(dict(group=lab,ici=f"{a}/{len(e)}",noici=f"{c_}/{len(u)}",OR=orci(a,len(e)-a,c_,len(u)-c_),fisher_p=round(fisher(a,len(e)-a,c_,len(u)-c_),4)))
out["imh"]=rows
# ratio of ORs soto NSCLC vs doce NSCLC (Haldane), crude
def lor(a,n1,c_,n0):
    b,d=n1-a,n0-c_
    if 0 in (a,b,c_,d): a,b,c_,d=a+.5,b+.5,c_+.5,d+.5
    return math.log((a/b)/(c_/d)), 1/a+1/b+1/c_+1/d
def ratio(r1,r2):
    l=r1[0]-r2[0]; se=math.sqrt(r1[1]+r2[1]); z=l/se; return [round(math.exp(l),2),round(math.exp(l-1.96*se),2),round(math.exp(l+1.96*se),2),round(math.erfc(abs(z)/math.sqrt(2)),3)]
def parse(s): a,n=s.split("/"); return int(a),int(n)
g={r["group"]:r for r in rows}
for k in ["sotorasib NSCLC","adagrasib NSCLC"]:
    a,n1=parse(g[k]["ici"]); c_,n0=parse(g[k]["noici"]); a2,m1=parse(g["docetaxel NSCLC"]["ici"]); c2,m0=parse(g["docetaxel NSCLC"]["noici"])
    out[f"imh_ratio_{k}_vs_docetaxel"]=ratio(lor(a,n1,c_,n0),lor(a2,m1,c2,m0))
# crude NSCLC hepatotoxicity ratio soto/doce and ada/doce
dh=set()
hep=set(l.strip().lower() for l in open(os.path.join(B,"scripts","hepato_pts.txt")) if l.strip())
for p in dn:
    if any(x.lower() in hep for x in pts.get(p,[])): dh.add(p)
e=[p for p in dn if p in ici_any]; u=[p for p in dn if p not in ici_any]
rd=lor(sum(p in dh for p in e),len(e),sum(p in dh for p in u),len(u))
for k in ["sotorasib","adagrasib"]:
    df=c[(c.target_drug==k)&(c.tumor=="NSCLC")]; ee=df[df.ici_exposed]; uu=df[~df.ici_exposed]
    out[f"crude_hep_ratio_{k}_vs_docetaxel_NSCLC"]=ratio(lor(int(ee.hepato.sum()),len(ee),int(uu.hepato.sum()),len(uu)),rd)
# adagrasib + ICI case details
t=c[(c.target_drug=="adagrasib")&c.ici_exposed]
det=[]
for r in t.itertuples():
    det.append(dict(pid=r.primaryid,hep=bool(r.hepato),ici=";".join(ici_any.get(r.primaryid,[])),interval=r.ici_to_kras_days,tto=r.tto_days,country=r.country,rept=r.rept_cod,quarter=r.quarter,tumor=r.tumor,LT=r.outc_LT,HO=r.outc_HO,DE=r.outc_DE,hep_pts=[x for x in pts.get(r.primaryid,[]) if x.lower() in hep]))
out["ada_ici_cases"]=det
s=c[(c.target_drug=="sotorasib")&c.ici_exposed]
out["soto_ici_hep_roles"]=pd.Series([";".join(sorted(set(x.split(":")[1] for x in ici_any.get(p,[])))) for p in s[s.hepato].primaryid]).value_counts().to_dict()
out["soto_ici_intervals_hep"]=s[s.hepato].ici_to_kras_days.replace("",np.nan).dropna().tolist()
json.dump(out,open(os.path.join(R,"round3_analyses.json"),"w"),indent=1,default=str); print(json.dumps(out,indent=1,default=str))
