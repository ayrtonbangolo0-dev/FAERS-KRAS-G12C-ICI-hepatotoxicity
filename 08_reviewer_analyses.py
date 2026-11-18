"""Step 8: analyses requested in review: interaction, adjusted model, Omega DDI, docetaxel reference group,
ICI role/trial splits, switch cases, interval groups, pooled time to onset, AEMS continuity."""
import pandas as pd, numpy as np, os, math, glob, json
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R=os.path.join(B,"results"); X=os.environ["FAERS_X"]; W=os.environ["FAERS_W"]; RS=os.path.expanduser("~/fw/rescan")
out={}
def pnorm(z): return 0.5*math.erfc(abs(z)/math.sqrt(2))*2
def ror(a,b,c,d):
    a,b,c,d=[x+0.5 if 0 in (a,b,c,d) else x for x in (a,b,c,d)]; r=(a/b)/(c/d); se=math.sqrt(1/a+1/b+1/c+1/d); return r,r*math.exp(-1.96*se),r*math.exp(1.96*se),se
c=pd.read_csv(os.path.join(R,"cohort_reports.csv"),dtype={"primaryid":str},keep_default_na=False)
for k in ["ici_exposed","hepato","outc_DE","outc_LT","outc_HO","serious","severe_hep"]: c[k]=c[k].astype(str)=="True"
demo=pd.read_csv(os.path.join(X,"demo.csv"),dtype=str,keep_default_na=False).drop_duplicates("primaryid").set_index("primaryid")
c["occp"]=c.primaryid.map(demo.occp_cod).fillna(""); c["quarter"]=c.primaryid.map(demo.quarter).fillna(""); c["year"]=c.quarter.str[:4]
drug=pd.read_csv(os.path.join(X,"drug.csv"),dtype=str,keep_default_na=False)
# ---- 1 ratio of RORs (interaction) crude
s=c[c.target_drug=="sotorasib"]; a=c[c.target_drug=="adagrasib"]
def cnt(d): return int(d[d.ici_exposed].hepato.sum()),int((d.ici_exposed&~d.hepato).sum()),int(d[~d.ici_exposed].hepato.sum()),int((~d.ici_exposed&~d.hepato).sum())
rs=ror(*cnt(s)); ra=ror(*cnt(a)); lr=math.log(ra[0]/rs[0]); se=math.sqrt(rs[3]**2+ra[3]**2)
out["ratio_of_RORs_ada_over_soto"]=dict(ratio=round(math.exp(lr),2),lo=round(math.exp(lr-1.96*se),2),hi=round(math.exp(lr+1.96*se),2),z=round(lr/se,2),p=round(pnorm(lr/se),3),soto=[round(x,2) for x in rs[:3]],ada=[round(x,2) for x in ra[:3]])
# ---- 2 logistic regression (IRLS)
def logit(df, cols):
    Xm=np.column_stack([np.ones(len(df))]+[df[k].values.astype(float) for k in cols]); y=df.hepato.values.astype(float); b=np.zeros(Xm.shape[1])
    for _ in range(50):
        eta=Xm@b; p=1/(1+np.exp(-eta)); Wt=p*(1-p); H=Xm.T@(Xm*Wt[:,None]); g=Xm.T@(y-p); step=np.linalg.solve(H+1e-8*np.eye(len(b)),g); b+=np.clip(step,-5,5)
        if np.max(np.abs(step))<1e-9: break
    cov=np.linalg.inv(H); se=np.sqrt(np.diag(cov))
    return pd.DataFrame({"term":["const"]+cols,"OR":np.exp(b),"lo":np.exp(b-1.96*se),"hi":np.exp(b+1.96*se),"p":[pnorm(x) for x in b/se]}), b, cov, ["const"]+cols
d=c.copy(); d["soto"]=(d.target_drug=="sotorasib").astype(int); d["ici"]=d.ici_exposed.astype(int); d["soto_x_ici"]=d.soto*d.ici; d["ada_ici"]=(1-d.soto)*d.ici; d["soto_ici"]=d.soto*d.ici
crude,_,_,_=logit(d,["soto","ici","soto_x_ici"])
# covariates
d["us"]=(d.country=="US").astype(int); d["cty_miss"]=(d.country=="").astype(int)
d["rept_per"]=(d.rept_cod=="PER").astype(int); d["rept_dir"]=(d.rept_cod=="DIR").astype(int)
d["consumer"]=(d.occp=="CN").astype(int); d["occp_miss"]=(d.occp=="").astype(int)
for y in ["2022","2023","2024","2025","2026"]: d[f"y{y}"]=(d.year==y).astype(int)
d["ind_crc"]=(d.tumor=="CRC").astype(int); d["ind_other"]=d.tumor.isin(["Pancreatic","Other"]).astype(int); d["ind_unk"]=(d.tumor=="Unknown").astype(int)
age=pd.to_numeric(d.age_yr,errors="coerce"); d["age65"]=(age>=65).astype(int); d["age_miss"]=age.isna().astype(int)
d["female"]=(d.sex=="F").astype(int); d["sex_miss"]=(~d.sex.isin(["F","M"])).astype(int)
cov=["us","rept_per","rept_dir","consumer","occp_miss","y2023","y2024","y2025","y2026","ind_crc","ind_other","ind_unk","age65","age_miss","female","sex_miss"]
adj,b,V,names=logit(d,["soto","ici","soto_x_ici"]+cov)
# ICI effect within sotorasib = ici + soto_x_ici
i1,i2=names.index("ici"),names.index("soto_x_ici"); est=b[i1]+b[i2]; se_=math.sqrt(V[i1,i1]+V[i2,i2]+2*V[i1,i2])
out["adjusted"]=dict(ICI_in_adagrasib=[round(math.exp(b[i1]),2),round(math.exp(b[i1]-1.96*math.sqrt(V[i1,i1])),2),round(math.exp(b[i1]+1.96*math.sqrt(V[i1,i1])),2)],
  ICI_in_sotorasib=[round(math.exp(est),2),round(math.exp(est-1.96*se_),2),round(math.exp(est+1.96*se_),2)],
  interaction_soto_x_ici=[round(adj.OR[i2],2),round(adj.lo[i2],2),round(adj.hi[i2],2),round(adj.p[i2],3)],
  soto_vs_ada_noICI=[round(adj.OR[names.index('soto')],2),round(adj.lo[names.index('soto')],2),round(adj.hi[names.index('soto')],2)])
s1=names.index("soto"); est=b[s1]+b[i2]; se_=math.sqrt(V[s1,s1]+V[i2,i2]+2*V[s1,i2]); out["adjusted"]["soto_vs_ada_withICI"]=[round(math.exp(est),2),round(math.exp(est-1.96*se_),2),round(math.exp(est+1.96*se_),2)]
out["crude_logit_interaction"]=[round(crude.OR[3],2),round(crude.lo[3],2),round(crude.hi[3],2),round(crude.p[3],3)]
adj.to_csv(os.path.join(R,"adjusted_logistic_model.csv"),index=False)
# ---- 3 database-level: Omega and docetaxel reference group
kept=set(open(os.path.join(W,"kept_primaryids.txt")).read().split())
hep=set(l.strip().lower() for l in open(os.path.join(B,"scripts","hepato_pts.txt")) if l.strip())
heprep=set()
with open(os.path.join(W,"reac_kept.tsv")) as f:
    for line in f:
        p,_,pt=line.rstrip("\n").partition("\t")
        if pt.lower() in hep: heprep.add(p)
ici_any=set(); ici_role={}
for fn in glob.glob(os.path.join(RS,"*_ici.tsv")):
    for line in open(fn):
        p,k,role=line.rstrip("\n").split("\t"); 
        if p in kept: ici_any.add(p); ici_role.setdefault(p,set()).add(role)
doce_ps=set(); doce_any=set()
for fn in glob.glob(os.path.join(RS,"*_doce.tsv")):
    for line in open(fn):
        p,role=line.rstrip("\n").split("\t")
        if p in kept:
            doce_any.add(p)
            if role=="PS": doce_ps.add(p)
soto_any=set(drug[drug.target_drug=="sotorasib"].primaryid)&kept; ada_any=set(drug[drug.target_drug=="adagrasib"].primaryid)&kept
N=len(kept)
def omega(Dset,Iset):
    both=Dset&Iset; donly=Dset-Iset; ionly=Iset-Dset
    n11=len(both); n10=len(donly); n01=len(ionly); n00=N-n11-n10-n01
    e11=len(both&heprep); e10=len(donly&heprep); e01=len(ionly&heprep); e00=len(heprep)-e11-e10-e01
    f00=e00/n00; f10=e10/n10; f01=e01/n01
    o=lambda f: f/(1-f)
    g11=1-1/(max(o(f00),o(f10))+max(o(f00),o(f01))-o(f00)+1)
    E=g11*n11; om=math.log2((e11+0.5)/(E+0.5)); lo=om-3.3*(e11+0.5)**-0.5-2*(e11+0.5)**-1.5
    return dict(n_both=n11,hep_both=e11,pct_both=round(100*e11/n11,1),pct_drug_only=round(100*f10,1),pct_ici_only=round(100*f01,1),pct_neither=round(100*f00,2),expected=round(E,1),omega=round(om,2),omega025=round(lo,2))
out["omega_sotorasib_x_ICI"]=omega(soto_any,ici_any); out["omega_adagrasib_x_ICI"]=omega(ada_any,ici_any); out["omega_docetaxel_x_ICI"]=omega(doce_any,ici_any)
# docetaxel (primary suspect) cohort: ICI vs no ICI
dp=doce_ps; a1=len(dp&ici_any&heprep); n1=len(dp&ici_any); a0=len((dp-ici_any)&heprep); n0=len(dp-ici_any)
r=ror(a1,n1-a1,a0,n0-a0); out["docetaxel_PS_ICI_vs_noICI"]=dict(ici=f"{a1}/{n1} ({100*a1/n1:.1f}%)",no_ici=f"{a0}/{n0} ({100*a0/n0:.1f}%)",ROR=[round(x,2) for x in r[:3]])
# ---- 4 ICI role and trial-code splits
c["ici_roles"]=c.primaryid.map(lambda p: ",".join(sorted(ici_role.get(p,[]))))
c["ici_suspect"]=c.ici_roles.str.contains("PS|SS|I"); c["ici_concomitant_only"]=c.ici_exposed & ~c.ici_suspect
codes=drug[(drug.target_drug!="")&drug.drugname.str.upper().str.contains("AMG|MRTX")&(drug.role_cod.str.upper()=="PS")].primaryid
c["trial_coded"]=c.primaryid.isin(set(codes))
rp={}
for fn in glob.glob(os.path.join(RS,"*_rpsr.tsv")):
    for line in open(fn):
        p,s_=line.rstrip("\n").split("\t"); rp.setdefault(p,set()).add(s_.strip())
c["rpsr"]=c.primaryid.map(lambda p: ",".join(sorted(rp.get(p,[]))))
c["study_report"]=c.trial_coded | c.rpsr.str.contains("SDY")
def within(df,mask,label):
    a_=int(df[mask].hepato.sum()); n1=int(mask.sum()); c_=int(df[~mask].hepato.sum()); n0=int((~mask).sum()); r=ror(a_,n1-a_,c_,n0-c_)
    return dict(comparison=label,exposed=f"{a_}/{n1}",reference=f"{c_}/{n0}",ROR=round(r[0],2),lo=round(r[1],2),hi=round(r[2],2))
rows=[]
for lab,sub in [("ICI role: concomitant only vs no ICI", c[~c.ici_suspect]), ("ICI role: suspect vs no ICI", c[c.ici_suspect | ~c.ici_exposed]),
                ("Excluding study/trial-coded reports", c[~c.study_report])]:
    for dname in ["sotorasib","adagrasib"]:
        t=sub[sub.target_drug==dname]; rows.append(dict(analysis=lab,drug=dname,**within(t,t.ici_exposed,"ICI vs no ICI")))
    t=sub[~sub.ici_exposed]; rows.append(dict(analysis=lab,drug="both",**within(t,t.target_drug=="sotorasib","Sotorasib vs adagrasib, no ICI")))
out["ici_role_counts"]={str(k):int(v) for k,v in c[c.ici_exposed].groupby(["target_drug","ici_suspect"]).size().items()}
out["study_report_counts"]={str(k):int(v) for k,v in c.groupby(["target_drug","study_report","ici_exposed"]).size().items()}
# ---- 5 switch cases: adagrasib PS reports also listing sotorasib
sw=set(drug[(drug.target_drug=="sotorasib")].primaryid)&set(c[c.target_drug=="adagrasib"].primaryid)
c["switch"]=c.primaryid.isin(sw)
out["switch_cases"]=dict(n=int(c.switch.sum()),hepato=int(c[c.switch].hepato.sum()),ici=int(c[c.switch].ici_exposed.sum()),ici_hep=int(c[c.switch&c.ici_exposed].hepato.sum()))
sub=c[~c.switch]
for dname in ["adagrasib"]:
    t=sub[sub.target_drug==dname]; rows.append(dict(analysis="Excluding switch cases",drug=dname,**within(t,t.ici_exposed,"ICI vs no ICI")))
t=sub[~sub.ici_exposed]; rows.append(dict(analysis="Excluding switch cases",drug="both",**within(t,t.target_drug=="sotorasib","Sotorasib vs adagrasib, no ICI")))
pd.DataFrame(rows).to_csv(os.path.join(R,"reviewer_sensitivity.csv"),index=False)
# ---- 6 interval groups
iv=pd.to_numeric(c.ici_to_kras_days,errors="coerce"); c["ivg"]=pd.cut(iv,[-1e9,-1,30,90,1e9],labels=["overlapping dates","0-30 d","31-90 d",">90 d"])
out["interval_groups"]={str(k):f"{int(g.hepato.sum())}/{len(g)}" for k,g in c[iv.notna()].groupby("ivg",observed=False)}
# ---- 7 time to onset all rows and country %
t_=pd.to_numeric(c.tto_days,errors="coerce"); c["tto"]=t_
tto={}
for dname in ["sotorasib","adagrasib"]:
    for lab,m in [("all",c.target_drug==dname),("no ICI",(c.target_drug==dname)&~c.ici_exposed),("ICI",(c.target_drug==dname)&c.ici_exposed)]:
        x=c[m&c.hepato].tto.dropna(); tto[f"{dname} {lab}"]=f"n={len(x)} median {x.median():.0f} IQR {x.quantile(.25):.0f}-{x.quantile(.75):.0f}"
out["tto"]=tto
out["country_all"]={dname:dict(US=round(100*(g.country=="US").mean(),1),other=round(100*((g.country!="US")&(g.country!="")).mean(),1)) for dname,g in c.groupby("target_drug")}
# ---- 8 AEMS continuity
dm=pd.read_csv(os.path.join(W,"all_demo.tsv"),sep="\t",header=None,names=["primaryid","caseid","fda_dt","q"],dtype=str)
pre=set(dm[dm.q<"2026Q1"].caseid); post=dm[dm.q>="2026Q1"]
out["aems"]=dict(reports_2026=int(len(post)),pct_2026_caseids_seen_before=round(100*post.caseid.isin(pre).mean(),1),
  primaryid_is_caseid_plus_version_2026=round(100*np.mean([str(a).startswith(str(b)) for a,b in zip(post.primaryid,post.caseid)]),1),
  primaryid_is_caseid_plus_version_pre=round(100*np.mean([str(a).startswith(str(b)) for a,b in zip(dm[dm.q<"2026Q1"].primaryid,dm[dm.q<"2026Q1"].caseid)]),1),
  deleted_files={os.path.basename(f):sum(1 for _ in open(f)) for f in sorted(glob.glob(os.path.join(B,"work","q","2026*_deleted.txt")))},
  reports_per_quarter=dm.groupby("q").size().to_dict())
c.to_csv(os.path.join(R,"cohort_reports_v2_flags.csv"),index=False)
json.dump(out,open(os.path.join(R,"reviewer_analyses.json"),"w"),indent=1,default=str)
print(json.dumps(out,indent=1,default=str)); print(pd.DataFrame(rows).to_string()); print(adj.round(3).to_string())
