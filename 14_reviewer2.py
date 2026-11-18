"""Step 14: Reviewer 2 analyses.
A any-role and suspect/interacting KRAS G12C inhibitor cohorts (ICI vs no ICI; interaction; NSCLC adjusted; pooled vs docetaxel)
B complete-case NSCLC models and missingness profile
C 'K-ras gene mutation' indication PT: counts and sensitivity analyses
D Firth implementation validation (closed form for saturated 2x2x2 model; finite-difference gradient of penalized log-likelihood)
Reuses data loading and model functions from step 11 (everything before its section A)."""
import os, json, math
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); S=os.path.join(B,"scripts")
src=open(os.path.join(S,"11_round2.py")).read().split("# ---- A crude")[0]
exec(src)   # defines c (PS cohort), doc (docetaxel PS), prep, fit, est, lrt, COV, IND, ror, cnt, kept, ici_any, heprep
import pandas as pd, numpy as np
def cnt(d): return int(d[d.ici_exposed].hepato.sum()),int((d.ici_exposed&~d.hepato).sum()),int(d[~d.ici_exposed].hepato.sum()),int((~d.ici_exposed&~d.hepato).sum())
o={}
M=os.path.join(W,"v2","m")
rd=lambda t: pd.read_csv(os.path.join(M,t+".csv"),dtype=str,keep_default_na=False)
drug,demo,indi=rd("drug"),rd("demo"),rd("indi")
drug=drug[drug.primaryid.isin(kept)]
kd=drug[drug.target_drug!=""].copy()
rank={"PS":0,"SS":1,"I":2,"C":3}; kd["rk"]=kd.role_cod.str.upper().map(rank).fillna(4)
mr=kd.groupby(["primaryid","target_drug"]).rk.min().reset_index(); top=mr.groupby("primaryid").rk.transform("min"); ties=mr[mr.rk==top].groupby("primaryid").size()
kd=kd[kd.primaryid.isin(ties[ties==1].index)]   # drug with the best role is assigned; reports with both drugs in the same best role excluded (as in the PS cohort)
best=kd.sort_values(["primaryid","rk","drug_seq"]).drop_duplicates("primaryid")
co=best[["primaryid","target_drug","drug_seq","rk"]].copy()
co["role"]=co.rk.map({0:"PS",1:"SS",2:"I",3:"C",4:"other"})
co["ici_exposed"]=co.primaryid.isin(ici_any); co["hepato"]=co.primaryid.isin(heprep)
dm=demo.drop_duplicates("primaryid").set_index("primaryid")
co["sex"]=co.primaryid.map(dm.sex); co["country"]=co.primaryid.map(dm.occr_country); co["rept_cod"]=co.primaryid.map(dm.rept_cod)
co["occp"]=co.primaryid.map(dm.occp_cod).fillna(""); co["year"]=co.primaryid.map(dm.quarter).str[:4]
co["age_yr"]=[agey(dm.at[p,"age"],dm.at[p,"age_cod"]) if p in dm.index else np.nan for p in co.primaryid]
ind=indi.merge(co[["primaryid","drug_seq"]],left_on=["primaryid","indi_drug_seq"],right_on=["primaryid","drug_seq"])
co["indication"]=co.primaryid.map(ind.groupby("primaryid").indi_pt.apply(lambda s:"; ".join(sorted(set(s))))).fillna("")
il=co.indication.str.lower()
co["tumor"]=np.select([il.str.contains("colorectal|colon|rectal|rectum|caecum|cecum|sigmoid"),il.str.contains("lung|nsclc|bronch"),il.str.contains("pancrea"),il==""],["CRC","NSCLC","Pancreatic","Unknown"],"Other")
# check: PS subset reproduces published cohort
ps=co[co.role=="PS"]
o["check_PS"]={d:dict(n=int((ps.target_drug==d).sum()),hep=int(ps[ps.target_drug==d].hepato.sum()),ici=int(ps[ps.target_drug==d].ici_exposed.sum()),nsclc=int(((ps.target_drug==d)&(ps.tumor=="NSCLC")).sum())) for d in ["sotorasib","adagrasib"]}
o["check_PS_match_published_ids"]=len(set(ps.primaryid)^set(c.primaryid))
o["role_counts"]={d:co[co.target_drug==d].role.value_counts().to_dict() for d in ["sotorasib","adagrasib"]}
o["role_by_ici"]={d:{r:dict(n=int(((co.target_drug==d)&(co.role==r)&co.ici_exposed).sum()),hep=int(co[(co.target_drug==d)&(co.role==r)&co.ici_exposed].hepato.sum())) for r in ["PS","SS","I","C"]} for d in ["sotorasib","adagrasib"]}
NS="unknown indication|K-ras gene mutation|Neoplasm malignant|Illness"
def recode(d):
    d=d.copy(); m=(d.tumor=="Other")&d.indication.str.contains(NS,case=False); d.loc[m,"tumor"]="Unknown"; return d
co=recode(co); c=recode(c); doc=recode(doc)
o["indication_coding"]=dict(kras_PS=c.tumor.value_counts().to_dict(),docetaxel=doc.tumor.value_counts().to_dict())
cp_=prep(co); dp=prep(doc)
def r3(x): return [round(v,2) for v in x[:3]]
def cohort_block(df,label):
    res={}
    for d in ["sotorasib","adagrasib"]:
        s=df[df.target_drug==d]; a_,b_,c_,d_=cnt(s)
        res[d]=dict(n=len(s),ici=f"{a_}/{a_+b_}",noici=f"{c_}/{c_+d_}",crude=r3(ror(a_,b_,c_,d_)))
        m=fit(s,["ici"]+COV+IND); res[d]["adj_all"]=est(m,["ici"])
        sn=s[s.tumor=="NSCLC"]; a2,b2,c2,d2=cnt(sn)
        res[d]["nsclc"]=dict(n=len(sn),ici=f"{a2}/{a2+b2}",noici=f"{c2}/{c2+d2}",crude=r3(ror(a2,b2,c2,d2)),adj=est(fit(sn,["ici"]+COV),["ici"]))
        pool=pd.concat([s,dp]); pool=pool[pool.tumor=="NSCLC"].copy(); pool["kras"]=(pool.target_drug==d).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
        cols=["kras","ici","kras_x_ici"]+COV; m1=fit(pool,cols); m0=fit(pool,[x for x in cols if x!="kras_x_ici"])
        res[d]["ratio_vs_docetaxel_nsclc"]=est(m1,["kras_x_ici"])+[lrt(m1,m0)]
    a_,b_,c_,d_=cnt(df); res["both_crude"]=r3(ror(a_,b_,c_,d_)); res["both_ici"]=f"{a_}/{a_+b_}"; res["both_noici"]=f"{c_}/{c_+d_}"
    df=df.copy(); df["ada"]=(df.target_drug=="adagrasib").astype(int); df["ada_x_ici"]=df.ada*df.ici
    for lab,cols in [("crude",["ada","ici","ada_x_ici"]),("adjusted",["ada","ici","ada_x_ici"]+COV+IND)]:
        for fi in [False,True]:
            m1=fit(df,cols,fi); m0=fit(df,[x for x in cols if x!="ada_x_ici"],fi)
            res[f"interaction_{lab}_{'firth' if fi else 'ml'}"]=est(m1,["ada_x_ici"])+[lrt(m1,m0)]
    return res
o["cohort_PS"]=cohort_block(cp_[cp_.role=="PS"],"PS")
o["cohort_suspect_interacting"]=cohort_block(cp_[cp_.role.isin(["PS","SS","I"])],"SI")
o["cohort_any_role"]=cohort_block(cp_,"any")
o["cohort_nonPS_only"]={d:dict(zip(["ici_hep","ici_n","noici_hep","noici_n"],[int(x) for x in (lambda t:(t[0],t[0]+t[1],t[2],t[2]+t[3]))(cnt(cp_[(cp_.role!="PS")&(cp_.target_drug==d)]))])) for d in ["sotorasib","adagrasib"]}
pc0=prep(c); pc0["ada"]=(pc0.target_drug=="adagrasib").astype(int); pc0["ada_x_ici"]=pc0.ada*pc0.ici
cols=["ada","ici","ada_x_ici"]+COV+IND; rc={}
for fi in [False,True]:
    m1=fit(pc0,cols,fi); m0=fit(pc0,[x for x in cols if x!="ada_x_ici"],fi); k="firth" if fi else "ml"
    rc[k]=dict(ratio=est(m1,["ada_x_ici"])+[lrt(m1,m0)],ici_soto=est(m1,["ici"]),ici_ada=est(m1,["ici","ada_x_ici"]),soto_vs_ada_noici=est(m1,["ada"]),soto_vs_ada_ici=est(m1,["ada","ada_x_ici"]))
rc["docetaxel_all_adjusted"]=est(fit(dp,["ici"]+COV+IND),["ici"])
for d in ["sotorasib","adagrasib"]:
    pool=pd.concat([pc0[pc0.target_drug==d],dp]).copy(); pool["kras"]=(pool.target_drug==d).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
    rc[f"pooled_all_{d}_vs_docetaxel"]=est(fit(pool,["kras","ici","kras_x_ici"]+COV+IND),["kras_x_ici"])
o["recoded_published_full_cohort"]=rc
# ---- B complete case (PS cohort, published variables)
pc=prep(c)
def complete(df): return df[pd.to_numeric(df.age_yr,errors="coerce").notna()&df.sex.isin(["F","M"])&(df.country!="")&(df.occp!="")]
COVcc=["us","rept_per","rept_dir","consumer","y2023","y2024","y2025","y2026","age65","female"]
cc={}
for d in ["sotorasib","adagrasib"]:
    s=pc[(pc.target_drug==d)&(pc.tumor=="NSCLC")]; s2=complete(s); a_,b_,c_,d_=cnt(s2)
    cc[d]=dict(n_nsclc=len(s),n_complete=len(s2),ici=f"{a_}/{a_+b_}",noici=f"{c_}/{c_+d_}",crude=r3(ror(a_,b_,c_,d_)),adj=est(fit(s2,["ici"]+COVcc),["ici"]))
    pool=pd.concat([s2,complete(dp[dp.tumor=="NSCLC"])]).copy(); pool["kras"]=(pool.target_drug==d).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
    cols=["kras","ici","kras_x_ici"]+COVcc; m1=fit(pool,cols); m0=fit(pool,[x for x in cols if x!="kras_x_ici"])
    cc[d]["ratio_vs_docetaxel"]=est(m1,["kras_x_ici"])+[lrt(m1,m0)]
dn=complete(dp[dp.tumor=="NSCLC"]); a_,b_,c_,d_=cnt(dn)
cc["docetaxel"]=dict(n_nsclc=int((dp.tumor=="NSCLC").sum()),n_complete=len(dn),ici=f"{a_}/{a_+b_}",noici=f"{c_}/{c_+d_}",adj=est(fit(dn,["ici"]+COVcc),["ici"]))
# full cohort interaction, complete case incl. known indication
fc=complete(pc[pc.tumor!="Unknown"]).copy(); fc["ada"]=(fc.target_drug=="adagrasib").astype(int); fc["ada_x_ici"]=fc.ada*fc.ici
cols=["ada","ici","ada_x_ici"]+COVcc+["ind_crc","ind_other"]; m1=fit(fc,cols); m0=fit(fc,[x for x in cols if x!="ada_x_ici"])
cc["interaction_full_complete"]=dict(n=len(fc),ici_events={d:int(fc[(fc.target_drug==d)&fc.ici_exposed].hepato.sum()) for d in ["sotorasib","adagrasib"]},ratio=est(m1,["ada_x_ici"])+[lrt(m1,m0)],ici_soto=est(m1,["ici"]),ici_ada=est(m1,["ici","ada_x_ici"]))
o["complete_case"]=cc
# missingness profile
mp={}
for d,df in [("sotorasib",pc[pc.target_drug=="sotorasib"]),("adagrasib",pc[pc.target_drug=="adagrasib"]),("docetaxel_nsclc",dp[dp.tumor=="NSCLC"])]:
    row={}
    for var,miss in [("age",pd.to_numeric(df.age_yr,errors="coerce").isna()),("sex",~df.sex.isin(["F","M"])),("country",df.country==""),("reporter",df.occp==""),("indication",df.tumor=="Unknown")]:
        row[var]=dict(n_missing=int(miss.sum()),pct_missing=round(100*miss.mean(),1),hep_missing=round(100*df[miss].hepato.mean(),1) if miss.any() else None,hep_known=round(100*df[~miss].hepato.mean(),1),ici_missing=round(100*df[miss].ici_exposed.mean(),1) if miss.any() else None,ici_known=round(100*df[~miss].ici_exposed.mean(),1))
    row["any_missing_age_sex_country_reporter"]=round(100*(1-len(complete(df))/len(df)),1)
    row["by_report_type_age_missing_pct"]=df.groupby("rept_cod").apply(lambda g: round(100*pd.to_numeric(g.age_yr,errors="coerce").isna().mean(),1)).to_dict()
    mp[d]=row
o["missingness"]=mp
# ---- C K-ras gene mutation PT
kr=pc.indication.str.lower().str.contains("k-ras|kras")
kronly=kr&(pc.tumor=="Unknown")&~pc.indication.str.lower().str.contains("unknown indication|neoplasm malignant|illness")
o["kras_pt"]=dict(any=int(kr.sum()),by_drug=pc[kr].target_drug.value_counts().to_dict(),as_only_tumor_info=int(kronly.sum()),only_by_drug=pc[kronly].target_drug.value_counts().to_dict(),
    with_nsclc=int((kr&(pc.tumor=="NSCLC")).sum()),with_crc=int((kr&(pc.tumor=="CRC")).sum()),
    only_hep=int(pc[kronly].hepato.sum()),only_ici=int(pc[kronly].ici_exposed.sum()),only_ici_hep=int(pc[kronly&pc.ici_exposed].hepato.sum()),
    only_country_top=pc[kronly].country.replace("","NR").value_counts().head(5).to_dict(),unknown_total=int((pc.tumor=="Unknown").sum()))
ks={}
for lab,df in [("exclude_kras_only",pc[~kronly]),("reclassify_kras_only_as_nsclc",pc.assign(tumor=np.where(kronly,"NSCLC",pc.tumor)))]:
    df=df.copy(); df["ind_unk"]=(df.tumor=="Unknown").astype(int); r={}
    for d in ["sotorasib","adagrasib"]:
        s=df[df.target_drug==d]; a_,b_,c_,d_=cnt(s); r[d]=dict(n=len(s),ici=f"{a_}/{a_+b_}",noici=f"{c_}/{c_+d_}",crude=r3(ror(a_,b_,c_,d_)))
        sn=s[s.tumor=="NSCLC"]; a2,b2,c2,d2=cnt(sn); r[d]["nsclc"]=dict(n=len(sn),ici=f"{a2}/{a2+b2}",noici=f"{c2}/{c2+d2}",adj=est(fit(sn,["ici"]+COV),["ici"]))
        pool=pd.concat([sn,dp[dp.tumor=="NSCLC"]]).copy(); pool["kras"]=(pool.target_drug==d).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
        cols=["kras","ici","kras_x_ici"]+COV; m1=fit(pool,cols); m0=fit(pool,[x for x in cols if x!="kras_x_ici"]); r[d]["ratio_vs_docetaxel"]=est(m1,["kras_x_ici"])+[lrt(m1,m0)]
    df["ada"]=(df.target_drug=="adagrasib").astype(int); df["ada_x_ici"]=df.ada*df.ici
    cols=["ada","ici","ada_x_ici"]+COV+IND; m1=fit(df,cols); m0=fit(df,[x for x in cols if x!="ada_x_ici"]); r["interaction_adjusted"]=est(m1,["ada_x_ici"])+[lrt(m1,m0)]
    ks[lab]=r
o["kras_sensitivity"]=ks
# ---- D Firth validation
pc["ada"]=(pc.target_drug=="adagrasib").astype(int); pc["ada_x_ici"]=pc.ada*pc.ici
m=fit(pc,["ada","ici","ada_x_ici"],True)
def cell(d,i): s=pc[(pc.target_drug==d)&(pc.ici_exposed==i)]; return s.hepato.sum()+0.5,(~s.hepato).sum()+0.5
lo=lambda d,i: math.log(cell(d,i)[0]/cell(d,i)[1])
closed=math.exp((lo("adagrasib",True)-lo("adagrasib",False))-(lo("sotorasib",True)-lo("sotorasib",False)))
fv=dict(saturated_firth_interaction=round(math.exp(m["b"][m["names"].index("ada_x_ici")]),6),closed_form_plus_half=round(closed,6))
# finite-difference gradient of penalized log-likelihood at the adjusted Firth solution
cols=["ada","ici","ada_x_ici"]+COV+IND; mf=fit(pc,cols,True); cols2=mf["names"][1:]
X=np.column_stack([np.ones(len(pc))]+[pc[k].values.astype(float) for k in cols2]); y=pc.hepato.values.astype(float)
def pll(b):
    p=1/(1+np.exp(-(X@b))); w=p*(1-p); H=X.T@(X*w[:,None]); return float(np.sum(y*np.log(p)+(1-y)*np.log(1-p)))+0.5*np.linalg.slogdet(H)[1]
g=[]; h=1e-6
for j in range(len(mf["b"])):
    e=np.zeros(len(mf["b"])); e[j]=h; g.append((pll(mf["b"]+e)-pll(mf["b"]-e))/(2*h))
fv["max_abs_fd_gradient_at_firth_solution"]=float(np.max(np.abs(g)))
fv["n_parameters"]=len(mf["b"])
# perturbation check: penalized log-likelihood lower at random nearby points
rng=np.random.default_rng(1); base=pll(mf["b"]); fv["all_perturbations_lower"]=bool(all(pll(mf["b"]+rng.normal(0,0.01,len(mf["b"])))<base for _ in range(200)))
# ML check: 2x2 single-covariate closed form
s=pc[pc.target_drug=="sotorasib"]; mm=fit(s,["ici"]); a_,b_,c_,d_=cnt(s)
fv["ml_single_binary_vs_closed"]=[round(math.exp(mm["b"][1]),6),round((a_*d_)/(b_*c_),6)]
o["firth_validation"]=fv

# ---- E docetaxel reference in matching roles (any role / suspect-interacting), from step 10b
DA=os.path.join(W,"doce_any"); role={}; seqd={}
for fn in glob.glob(os.path.join(DA,"*_role.tsv")):
    for line in open(fn):
        p,k,sq=line.rstrip("\n").split("\t")
        if p in kept: role[p]=int(k)
dmA={}; indA={}
for fn in glob.glob(os.path.join(DA,"*_demo.tsv")):
    for line in open(fn):
        f=line.rstrip("\n").split("\t")
        if f[0] in role: dmA[f[0]]=f
for fn in glob.glob(os.path.join(DA,"*_indi.tsv")):
    for line in open(fn):
        p,pt=line.rstrip("\n").split("\t")
        if p in role: indA.setdefault(p,set()).add(pt)
rows=[]
for p,k in role.items():
    f=dmA.get(p,[p,"","","","","","",""]); il_="; ".join(sorted(indA.get(p,[]))).lower()
    if ("small cell lung" in il_ and "non-small" not in il_ and "non small" not in il_): tum="Other"
    elif any(x in il_ for x in ["colorectal","colon","rectal","rectum"]): tum="CRC"
    elif any(x in il_ for x in ["lung","nsclc","bronch"]): tum="NSCLC"
    elif il_=="": tum="Unknown"
    else: tum="Other"
    rows.append(dict(primaryid=p,target_drug="docetaxel",rk=k,ici_exposed=p in ici_any,hepato=p in heprep,age_yr=agey(f[1],f[2]),sex=f[3],country=f[4],rept_cod=f[5],occp=f[6],year=f[7][:4],tumor=tum,indication=il_))
dA=recode(pd.DataFrame(rows)); dAp=prep(dA)
o["docetaxel_any_check_PS_n"]=int((dA.rk==0).sum())
ref={}
for lab,rk in [("suspect_interacting",[0,1,2]),("any_role",[0,1,2,3])]:
    dd=dAp[dAp.rk.isin(rk)&(dAp.tumor=="NSCLC")]; a_,b_,c_,d_=cnt(dd)
    r={"docetaxel_nsclc":dict(n=len(dd),ici=f"{a_}/{a_+b_}",noici=f"{c_}/{c_+d_}",crude=r3(ror(a_,b_,c_,d_)),adj=est(fit(dd,["ici"]+COV),["ici"]))}
    kc=cp_[cp_.rk.isin(rk)]
    for d in ["sotorasib","adagrasib"]:
        pool=pd.concat([kc[(kc.target_drug==d)&(kc.tumor=="NSCLC")],dd]).copy(); pool["kras"]=(pool.target_drug==d).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
        cols=["kras","ici","kras_x_ici"]+COV; m1=fit(pool,cols); m0=fit(pool,[x for x in cols if x!="kras_x_ici"])
        r[f"{d}_ratio_vs_docetaxel_same_role"]=est(m1,["kras_x_ici"])+[lrt(m1,m0)]
    ref[lab]=r
o["docetaxel_matching_role"]=ref
# exact sotorasib-vs-adagrasib adjusted ORs, recoded indication
pc0["soto"]=(pc0.target_drug=="sotorasib").astype(int); pc0["soto_x_ici"]=pc0.soto*pc0.ici
m1=fit(pc0,["soto","ici","soto_x_ici"]+COV+IND)
o["recoded_soto_vs_ada"]=dict(noici=est(m1,["soto"]),ici=est(m1,["soto","soto_x_ici"]))
json.dump(o,open(os.path.join(R,"reviewer2_analyses.json"),"w"),indent=1,default=str)

co.to_csv(os.path.join(R,"cohort_any_role.csv"),index=False)
print(json.dumps(o,indent=1,default=str))
