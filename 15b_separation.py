"""Step 15b: identify covariates with zero-event levels (separation) in NSCLC models and refit with Firth's penalized likelihood."""
import os, json, math
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); S=os.path.join(B,"scripts")
exec(open(os.path.join(S,"11_round2.py")).read().split("# ---- A crude")[0])
import pandas as pd, numpy as np
NS="unknown indication|K-ras gene mutation|Neoplasm malignant|Illness"
def recode(d):
    d=d.copy(); m=(d.tumor=="Other")&d.indication.str.contains(NS,case=False); d.loc[m,"tumor"]="Unknown"; return d
cp=prep(recode(c)); dp=prep(recode(doc)); o={}
def sep(df):
    r={}
    for k in COV:
        if df[k].sum()==0 or df[k].sum()==len(df): continue
        e1=int(df[df[k]==1].hepato.sum()); n1=int(df[k].sum())
        if e1==0 or e1==n1: r[k]=f"{e1}/{n1}"
    return r
for lab,df in [("adagrasib_nsclc",cp[(cp.target_drug=="adagrasib")&(cp.tumor=="NSCLC")]),("docetaxel_nsclc",dp[dp.tumor=="NSCLC"]),("sotorasib_nsclc",cp[(cp.target_drug=="sotorasib")&(cp.tumor=="NSCLC")])]:
    o[lab]=dict(zero_event_levels=sep(df),ml=est(fit(df,["ici"]+COV),["ici"]),firth=est(fit(df,["ici"]+COV,True),["ici"]))
for d in ["sotorasib","adagrasib"]:
    pool=pd.concat([cp[(cp.target_drug==d)&(cp.tumor=="NSCLC")],dp[dp.tumor=="NSCLC"]]).copy(); pool["kras"]=(pool.target_drug==d).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
    cols=["kras","ici","kras_x_ici"]+COV
    m1=fit(pool,cols,True); m0=fit(pool,[x for x in cols if x!="kras_x_ici"],True)
    o[f"pooled_{d}_firth"]=est(m1,["kras_x_ici"])+[lrt(m1,m0)]
    o[f"pooled_{d}_zero_event_levels"]=sep(pool)
json.dump(o,open(os.path.join(R,"editorial_separation.json"),"w"),indent=1); print(json.dumps(o,indent=1))
