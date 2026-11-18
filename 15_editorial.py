"""Step 15 (editorial round): study flow counts and model diagnostics (events per parameter, convergence, separation checks, c statistic).
Reuses step 11 data loading/model functions and step 14 recoding."""
import os, json, math
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); S=os.path.join(B,"scripts")
exec(open(os.path.join(S,"11_round2.py")).read().split("# ---- A crude")[0])
import pandas as pd, numpy as np
o={}
# ---- flow
W_=W
o["raw_reports"]=int(open(os.path.join(W_,"n_raw.txt")).read())
cl=sum(1 for _ in open(os.path.join(W_,"case_latest.tsv")))
o["case_latest_lines"]=cl
o["deleted_caseids_listed"]=sum(1 for _ in open(os.path.join(W_,"deleted_caseids.txt")))
o["kept"]=len(kept)
M=os.path.join(W,"v2","m"); drug=pd.read_csv(os.path.join(M,"drug.csv"),dtype=str,keep_default_na=False)
dk=drug[drug.primaryid.isin(kept)&(drug.target_drug!="")]
o["kras_any_role_reports"]=int(dk.primaryid.nunique())
o["kras_any_role_by_drug"]=dk.groupby("target_drug").primaryid.nunique().to_dict()
ps=dk[dk.role_cod.str.upper()=="PS"]; nps=ps.groupby("primaryid").target_drug.nunique()
o["kras_ps_reports"]=int(ps.primaryid.nunique()); o["both_ps_excluded"]=int((nps>1).sum())
o["not_ps_excluded"]=o["kras_any_role_reports"]-o["kras_ps_reports"]
o["cohort_ps"]=len(c)
# reports in raw KRAS extraction before dedup
o["kras_raw_reports_before_dedup"]=int(drug[drug.target_drug!=""].primaryid.nunique())
# ---- diagnostics
NS="unknown indication|K-ras gene mutation|Neoplasm malignant|Illness"
def recode(d):
    d=d.copy(); m=(d.tumor=="Other")&d.indication.str.contains(NS,case=False); d.loc[m,"tumor"]="Unknown"; return d
cp=prep(recode(c)); dp=prep(recode(doc))
def diag(df,cols,label,firth=False):
    cols=[k for k in cols if df[k].sum()>0 and df[k].sum()<len(df)]
    X=np.column_stack([np.ones(len(df))]+[df[k].values.astype(float) for k in cols]); y=df.hepato.values.astype(float); b=np.zeros(X.shape[1]); it=0
    for it in range(1,201):
        p=1/(1+np.exp(-(X@b))); w=p*(1-p); H=X.T@(X*w[:,None]); Hi=np.linalg.inv(H+1e-10*np.eye(len(b)))
        g=X.T@(y-p+(np.einsum("ij,jk,ik->i",X*np.sqrt(w)[:,None],Hi,X*np.sqrt(w)[:,None])*(0.5-p) if firth else 0))
        st=np.clip(Hi@g,-2,2); b=b+st
        if np.max(np.abs(st))<1e-9: break
    p=1/(1+np.exp(-(X@b))); se=np.sqrt(np.diag(np.linalg.inv(X.T@(X*(p*(1-p))[:,None]))))
    ev=int(y.sum()); ne=int(len(y)-ev); k=len(cols)
    # c statistic
    order=np.argsort(p); ranks=np.empty(len(p)); ranks[order]=np.arange(1,len(p)+1)
    auc=(ranks[y==1].sum()-ev*(ev+1)/2)/(ev*ne)
    # condition number of standardized design
    Xs=X[:,1:]; Xs=(Xs-Xs.mean(0))/Xs.std(0); cond=float(np.linalg.cond(Xs)) if k>0 else 1.0
    return dict(model=label,n=len(y),events=ev,nonevents=ne,parameters=k,epv=round(min(ev,ne)/k,1),iterations=it,converged=bool(it<200),max_abs_coef=round(float(np.max(np.abs(b[1:]))),2),max_se=round(float(np.max(se[1:])),2),c_statistic=round(float(auc),3),condition_number=round(cond,1))
cp["ada"]=(cp.target_drug=="adagrasib").astype(int); cp["ada_x_ici"]=cp.ada*cp.ici
D=[]
D.append(diag(cp,["ada","ici","ada_x_ici"]+COV+IND,"Interaction, all indications, adjusted"))
D.append(diag(cp,["ada","ici","ada_x_ici"]+COV+IND,"Interaction, all indications, Firth",True))
for d in ["sotorasib","adagrasib"]:
    s=cp[(cp.target_drug==d)&(cp.tumor=="NSCLC")]; D.append(diag(s,["ici"]+COV,f"{d.capitalize()}, NSCLC, adjusted"))
D.append(diag(dp[dp.tumor=="NSCLC"],["ici"]+COV,"Docetaxel, NSCLC, adjusted"))
for d in ["sotorasib","adagrasib"]:
    pool=pd.concat([cp[(cp.target_drug==d)&(cp.tumor=="NSCLC")],dp[dp.tumor=="NSCLC"]]).copy(); pool["kras"]=(pool.target_drug==d).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
    D.append(diag(pool,["kras","ici","kras_x_ici"]+COV,f"Pooled NSCLC, {d} vs docetaxel"))
o["diagnostics"]=D
# events in ICI-exposed cells
o["ici_cells"]={d:dict(ici_hep=int(cp[(cp.target_drug==d)&cp.ici_exposed].hepato.sum()),ici_nohep=int((cp[(cp.target_drug==d)&cp.ici_exposed].hepato==False).sum()),nsclc_ici_hep=int(cp[(cp.target_drug==d)&cp.ici_exposed&(cp.tumor=="NSCLC")].hepato.sum())) for d in ["sotorasib","adagrasib"]}
# unspecified indication by drug
o["unspecified_indication"]={d:dict(n=int(((cp.target_drug==d)&(cp.tumor=="Unknown")).sum()),pct=round(100*((cp.target_drug==d)&(cp.tumor=="Unknown")).sum()/(cp.target_drug==d).sum(),1)) for d in ["sotorasib","adagrasib"]}
json.dump(o,open(os.path.join(R,"editorial_round.json"),"w"),indent=1)
print(json.dumps(o,indent=1))
