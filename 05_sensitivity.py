"""Step 5: sensitivity analyses of within-cohort comparisons (and full-database ROR for the narrow definition)."""
import pandas as pd, numpy as np, os, math, json
B = os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R = os.path.join(B,"results"); M = os.path.join(B,"manuscript")
W = os.environ.get("FAERS_W", os.path.join(B,"work")); X = os.environ.get("FAERS_X", os.path.join(B,"extracted"))
c = pd.read_csv(os.path.join(R,"cohort_reports.csv"), keep_default_na=False)
for k in ["ici_exposed","hepato"]: c[k] = c[k].astype(str)=="True"
demo = pd.read_csv(os.path.join(X,"demo.csv"), dtype=str, keep_default_na=False).drop_duplicates("primaryid").set_index("primaryid")
c["primaryid"] = c.primaryid.astype(str)
c["occp"] = c.primaryid.map(demo.occp_cod).fillna(""); c["quarter"] = c.primaryid.map(demo.quarter).fillna("")
reac = pd.read_csv(os.path.join(X,"reac.csv"), dtype=str, keep_default_na=False)
NARROW = {s.lower() for s in ["Hepatotoxicity","Drug-induced liver injury","Liver injury","Hepatitis","Hepatitis acute","Hepatitis fulminant","Hepatitis cholestatic",
  "Immune-mediated hepatitis","Immune-mediated hepatic disorder","Autoimmune hepatitis","Hepatic failure","Acute hepatic failure","Acute on chronic liver failure",
  "Hepatocellular injury","Cholestatic liver injury","Mixed liver injury","Hepatic cytolysis","Hepatic necrosis","Jaundice","Immune-mediated cholangitis"]}
c["hep_narrow"] = c.primaryid.isin(reac[reac.pt.str.strip().str.lower().isin(NARROW)].primaryid)
def ror(a,b,cc,d):
    a,b,cc,d = [x+0.5 if 0 in (a,b,cc,d) else x for x in (a,b,cc,d)]
    r=(a/b)/(cc/d); se=math.sqrt(1/a+1/b+1/cc+1/d); return r, r*math.exp(-1.96*se), r*math.exp(1.96*se)
def within(s, out):
    rows=[]
    def comp(label, sub, mask):
        a=int(sub[mask][out].sum()); n1=int(mask.sum()); cc=int(sub[~mask][out].sum()); n0=int((~mask).sum())
        if n1==0 or n0==0: rows.append(dict(comparison=label, exposed=f"{a}/{n1}", reference=f"{cc}/{n0}", ROR=np.nan, lo=np.nan, hi=np.nan)); return
        r,l,h=ror(a,n1-a,cc,n0-cc); rows.append(dict(comparison=label, exposed=f"{a}/{n1}", reference=f"{cc}/{n0}", ROR=r, lo=l, hi=h))
    comp("ICI vs no ICI, either drug", s, s.ici_exposed)
    for d in ["sotorasib","adagrasib"]:
        t=s[s.target_drug==d]; comp(f"ICI vs no ICI, {d}", t, t.ici_exposed)
    for g,lab in [(False,"no ICI"),(True,"with ICI")]:
        t=s[s.ici_exposed==g]; comp(f"Sotorasib vs adagrasib, {lab}", t, t.target_drug=="sotorasib")
    return rows
HP = {"MD","PH","HP","OT","RN"}
sets = {"Main analysis": (c,"hepato"),
        "S1 NSCLC indication only": (c[c.tumor=="NSCLC"],"hepato"),
        "S2 Health-professional reporters only": (c[c.occp.isin(HP)],"hepato"),
        "S3 Clinical liver injury terms only": (c,"hep_narrow"),
        "S4 Reports received 2023 Q1 or later": (c[c.quarter>="2023Q1"],"hepato")}
out=[]
for name,(s,o) in sets.items():
    for r in within(s,o): out.append(dict(analysis=name, n_reports=len(s), **r))
T = pd.DataFrame(out); T.to_csv(os.path.join(R,"sensitivity_within_cohort.csv"), index=False)
# full-database ROR for narrow definition
kept_hep_narrow = None
rk = os.path.join(W,"reac_kept.tsv")
if os.path.exists(rk):
    ids=set()
    with open(rk) as f:
        for line in f:
            p,_,pt=line.rstrip("\n").partition("\t")
            if pt.lower() in NARROW: ids.add(p)
    N=len(open(os.path.join(W,"kept_primaryids.txt")).read().split())
    Nh=len(ids); res=[]
    for d in ["sotorasib","adagrasib"]:
        for g in [None,False,True]:
            s=c[c.target_drug==d] if g is None else c[(c.target_drug==d)&(c.ici_exposed==g)]
            a=int(s.hep_narrow.sum()); n1=len(s); cc=Nh-int(c.hep_narrow.sum()); d0=(N-len(c))-cc
            r,l,h=ror(a,n1-a,cc,d0); res.append(dict(drug=d, ici=("all" if g is None else ("ICI+" if g else "ICI-")), n=n1, hep=a, pct=round(100*a/n1,1), ROR=r, lo=l, hi=h))
    pd.DataFrame(res).to_csv(os.path.join(R,"sensitivity_narrow_vs_database.csv"), index=False); print(pd.DataFrame(res).round(2).to_string())
pd.set_option("display.width",200); print(T.round(2).to_string())
print(c.groupby("occp").size().to_dict())
