"""Step 11: second review round: docetaxel reference (NSCLC, adjusted, pooled), primary-suspect Omega,
Firth penalized logistic check, development-code reports profile and duplicate check, MRTX849 count."""
import pandas as pd, numpy as np, os, math, glob, json
from collections import Counter
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R=os.path.join(B,"results"); W=os.path.expanduser("~/fw"); RS=os.path.join(W,"rescan"); DC=os.path.join(W,"doce")
out={}
def pnorm(z): return math.erfc(abs(z)/math.sqrt(2))
def ror(a,b,c,d):
    a,b,c,d=[x+0.5 if 0 in (a,b,c,d) else x for x in (a,b,c,d)]; r=(a/b)/(c/d); se=math.sqrt(1/a+1/b+1/c+1/d); return r,r*math.exp(-1.96*se),r*math.exp(1.96*se),se
c=pd.read_csv(os.path.join(R,"cohort_reports_v2_flags.csv"),dtype={"primaryid":str},keep_default_na=False)
for k in ["ici_exposed","hepato","trial_coded","study_report","switch","ici_suspect"]: c[k]=c[k].astype(str)=="True"
kept=set(open(os.path.join(W,"kept_primaryids.txt")).read().split())
hep=set(l.strip().lower() for l in open(os.path.join(B,"scripts","hepato_pts.txt")) if l.strip())
cohort=set(c.primaryid); heprep=set(); pts={}
with open(os.path.join(W,"reac_kept.tsv")) as f:
    for line in f:
        p,_,pt=line.rstrip("\n").partition("\t")
        if pt.lower() in hep: heprep.add(p)
        if p in cohort: pts.setdefault(p,set()).add(pt)
ici_any=set()
for fn in glob.glob(os.path.join(RS,"*_ici.tsv")):
    for line in open(fn):
        p=line.split("\t")[0]
        if p in kept: ici_any.add(p)
doce_any=set(); doce_ps=set()
for fn in glob.glob(os.path.join(RS,"*_doce.tsv")):
    for line in open(fn):
        p,role=line.rstrip("\n").split("\t")
        if p in kept:
            doce_any.add(p)
            if role=="PS": doce_ps.add(p)
# docetaxel covariates
dm={}
for fn in glob.glob(os.path.join(DC,"*_demo.tsv")):
    for line in open(fn):
        f=line.rstrip("\n").split("\t")
        if f[0] in doce_ps: dm[f[0]]=f
ind={}
for fn in glob.glob(os.path.join(DC,"*_indi.tsv")):
    for line in open(fn):
        p,pt=line.rstrip("\n").split("\t")
        if p in doce_ps: ind.setdefault(p,set()).add(pt)
def agey(a,u):
    try: a=float(a)
    except: return np.nan
    u=u.upper(); return a if u=="YR" else a/12 if u=="MON" else a/52 if u=="WK" else a/365 if u=="DY" else a*10 if u=="DEC" else np.nan
rows=[]
for p in doce_ps:
    f=dm.get(p,[p,"","","","","","",""])
    il="; ".join(sorted(ind.get(p,[]))).lower()
    if ("small cell lung" in il and "non-small" not in il and "non small" not in il): tum="Other"
    elif any(k in il for k in ["colorectal","colon","rectal","rectum"]): tum="CRC"
    elif any(k in il for k in ["lung","nsclc","bronch"]): tum="NSCLC"
    elif il=="": tum="Unknown"
    else: tum="Other"
    rows.append(dict(primaryid=p,target_drug="docetaxel",ici_exposed=p in ici_any,hepato=p in heprep,age_yr=agey(f[1],f[2]),sex=f[3],country=f[4],rept_cod=f[5],occp=f[6],year=f[7][:4],tumor=tum,indication=il))
doc=pd.DataFrame(rows)
out["docetaxel_ps_n"]=len(doc); out["docetaxel_indications"]={k:int(v) for k,v in doc.tumor.value_counts().items()}
top=Counter(); [top.update(s.split("; ")) for s in doc.indication if s]; out["docetaxel_top_indication_pts"]=top.most_common(12)
def prep(d):
    d=d.copy()
    d["ici"]=d.ici_exposed.astype(int); d["us"]=(d.country=="US").astype(int)
    d["rept_per"]=(d.rept_cod=="PER").astype(int); d["rept_dir"]=(d.rept_cod=="DIR").astype(int)
    d["consumer"]=(d.occp=="CN").astype(int); d["occp_miss"]=(d.occp=="").astype(int)
    for y in ["2023","2024","2025","2026"]: d[f"y{y}"]=(d.year.astype(str)==y).astype(int)
    age=pd.to_numeric(d.age_yr,errors="coerce"); d["age65"]=(age>=65).astype(int); d["age_miss"]=age.isna().astype(int)
    d["female"]=(d.sex=="F").astype(int); d["sex_miss"]=(~d.sex.isin(["F","M"])).astype(int)
    d["ind_crc"]=(d.tumor=="CRC").astype(int); d["ind_other"]=d.tumor.isin(["Pancreatic","Other"]).astype(int); d["ind_unk"]=(d.tumor=="Unknown").astype(int)
    return d
COV=["us","rept_per","rept_dir","consumer","occp_miss","y2023","y2024","y2025","y2026","age65","age_miss","female","sex_miss"]
IND=["ind_crc","ind_other","ind_unk"]
def fit(df, cols, firth=False):
    cols=[k for k in cols if df[k].sum()>0 and df[k].sum()<len(df)]
    X=np.column_stack([np.ones(len(df))]+[df[k].values.astype(float) for k in cols]); y=df.hepato.values.astype(float); b=np.zeros(X.shape[1])
    for it in range(200):
        p=1/(1+np.exp(-(X@b))); w=p*(1-p); H=X.T@(X*w[:,None]); Hi=np.linalg.inv(H+1e-10*np.eye(len(b)))
        if firth:
            h=np.einsum("ij,jk,ik->i",X*np.sqrt(w)[:,None],Hi,X*np.sqrt(w)[:,None]); g=X.T@(y-p+h*(0.5-p))
        else: g=X.T@(y-p)
        step=np.clip(Hi@g,-2,2); b=b+step
        if np.max(np.abs(step))<1e-9: break
    p=1/(1+np.exp(-(X@b))); w=p*(1-p); H=X.T@(X*w[:,None]); V=np.linalg.inv(H)
    ll=float(np.sum(y*np.log(p)+(1-y)*np.log(1-p)))+(0.5*np.linalg.slogdet(H)[1] if firth else 0)
    return dict(b=b,V=V,names=["const"]+cols,ll=ll)
def est(m, terms):
    idx=[m["names"].index(t) for t in terms]; e=sum(m["b"][i] for i in idx); v=sum(m["V"][i,j] for i in idx for j in idx); s=math.sqrt(v)
    return [round(math.exp(e),2),round(math.exp(e-1.96*s),2),round(math.exp(e+1.96*s),2),round(pnorm(e/s),4)]
def lrt(m1,m0): 
    from math import erfc, sqrt
    x=2*(m1["ll"]-m0["ll"]); return round(erfc(sqrt(max(x,0)/2)),4)  # chi2 1 df
# ---- A crude ratio of ICI RORs sotorasib / docetaxel
def cnt(d): return int(d[d.ici_exposed].hepato.sum()),int((d.ici_exposed&~d.hepato).sum()),int(d[~d.ici_exposed].hepato.sum()),int((~d.ici_exposed&~d.hepato).sum())
s=c[c.target_drug=="sotorasib"]; a=c[c.target_drug=="adagrasib"]
rs=ror(*cnt(s)); rd=ror(*cnt(doc)); ra=ror(*cnt(a))
for lab,r1 in [("soto_over_doce",rs),("ada_over_doce",ra)]:
    lr=math.log(r1[0]/rd[0]); se=math.sqrt(r1[3]**2+rd[3]**2); out["crude_ratio_"+lab]=[round(math.exp(lr),2),round(math.exp(lr-1.96*se),2),round(math.exp(lr+1.96*se),2),round(pnorm(lr/se),3)]
# ---- B docetaxel and sotorasib NSCLC, adjusted ICI OR
cp=prep(c); dp=prep(doc)
res={}
for lab,df in [("docetaxel_NSCLC",dp[dp.tumor=="NSCLC"]),("sotorasib_NSCLC",cp[(cp.target_drug=="sotorasib")&(cp.tumor=="NSCLC")]),("adagrasib_NSCLC",cp[(cp.target_drug=="adagrasib")&(cp.tumor=="NSCLC")]),("docetaxel_all",dp)]:
    r=ror(*cnt(df)); m=fit(df,["ici"]+COV+([] if "NSCLC" in lab else IND))
    res[lab]=dict(n=len(df),ici=f"{int(df[df.ici_exposed].hepato.sum())}/{int(df.ici_exposed.sum())}",noici=f"{int(df[~df.ici_exposed].hepato.sum())}/{int((~df.ici_exposed).sum())}",crude=[round(x,2) for x in r[:3]],adjusted=est(m,["ici"]))
out["ici_or_by_drug"]=res
# ---- C pooled models sotorasib vs docetaxel (and adagrasib vs docetaxel)
for other in ["sotorasib","adagrasib"]:
    for scope in ["NSCLC","all"]:
        pool=pd.concat([cp[cp.target_drug==other],dp]); pool=pool[pool.tumor=="NSCLC"] if scope=="NSCLC" else pool
        pool["kras"]=(pool.target_drug==other).astype(int); pool["kras_x_ici"]=pool.kras*pool.ici
        cols=["kras","ici","kras_x_ici"]+COV+([] if scope=="NSCLC" else IND)
        m1=fit(pool,cols); m0=fit(pool,[x for x in cols if x!="kras_x_ici"])
        out[f"pooled_{other}_vs_docetaxel_{scope}"]=dict(n=len(pool),interaction=est(m1,["kras_x_ici"]),lrt_p=lrt(m1,m0),ici_in_docetaxel=est(m1,["ici"]),ici_in_kras=est(m1,["ici","kras_x_ici"]))
# ---- D primary-suspect Omega (drug as PS; reports with drug in other roles removed)
N_all=len(kept)
def omega_ps(D_ps, D_any):
    excl=D_any-D_ps; N=N_all-len(excl)
    both=D_ps&ici_any; donly=D_ps-ici_any; ionly=ici_any-D_any
    n11,n10,n01=len(both),len(donly),len(ionly); n00=N-n11-n10-n01
    e11,e10,e01=len(both&heprep),len(donly&heprep),len(ionly&heprep); e00=len(heprep-excl)-e11-e10-e01
    o=lambda f:f/(1-f); f00,f10,f01=e00/n00,e10/n10,e01/n01
    g=1-1/(max(o(f00),o(f10))+max(o(f00),o(f01))-o(f00)+1); E=g*n11; om=math.log2((e11+0.5)/(E+0.5)); lo=om-3.3*(e11+0.5)**-0.5-2*(e11+0.5)**-1.5
    return dict(n_both=n11,hep_both=e11,pct_both=round(100*e11/n11,1),pct_drug_only=round(100*f10,1),pct_ici_only=round(100*f01,1),pct_neither=round(100*f00,2),expected=round(E,1),omega=round(om,2),omega025=round(lo,2))
drug=pd.read_csv(os.path.join(W,"v2","m","drug.csv"),dtype=str,keep_default_na=False)
soto_any=set(drug[drug.target_drug=="sotorasib"].primaryid)&kept; ada_any=set(drug[drug.target_drug=="adagrasib"].primaryid)&kept
out["omega_ps_sotorasib"]=omega_ps(set(s.primaryid),soto_any); out["omega_ps_adagrasib"]=omega_ps(set(a.primaryid),ada_any); out["omega_ps_docetaxel"]=omega_ps(doce_ps,doce_any)
# ---- E Firth check of cohort models
cp["soto"]=(cp.target_drug=="sotorasib").astype(int); cp["soto_x_ici"]=cp.soto*cp.ici
fr={}
for lab,cols in [("crude",["soto","ici","soto_x_ici"]),("adjusted",["soto","ici","soto_x_ici"]+COV[:5]+["y2023","y2024","y2025","y2026"]+IND+COV[9:])]:
    for firth in [False,True]:
        m1=fit(cp,cols,firth); m0=fit(cp,[x for x in cols if x!="soto_x_ici"],firth)
        fr[f"{lab}_{'firth' if firth else 'ml'}"]=dict(ici_in_adagrasib=est(m1,["ici"]),ici_in_sotorasib=est(m1,["ici","soto_x_ici"]),interaction=est(m1,["soto_x_ici"]),lrt_p=lrt(m1,m0),soto_vs_ada_noici=est(m1,["soto"]),soto_vs_ada_ici=est(m1,["soto","soto_x_ici"]))
out["firth"]=fr
out["covariate_cells_ada_ici"]={k:int(cp[(cp.target_drug=="adagrasib")&cp.ici_exposed][k].sum()) for k in COV+IND}
# ---- F development-code reports
t=c[c.trial_coded]; b_=c[(c.target_drug=="sotorasib")&~c.trial_coded]
def prof(d): return dict(n=len(d),hepato=int(d.hepato.sum()),rept=d.rept_cod.value_counts().to_dict(),country=d.country.replace("","NR").value_counts().head(8).to_dict(),year=d.year.astype(str).value_counts().sort_index().to_dict(),occp=d.occp.replace("","NR").value_counts().to_dict(),ici=int(d.ici_exposed.sum()),tumor=d.tumor.value_counts().to_dict(),
    top_pts=Counter(pt for p in d.primaryid for pt in pts.get(p,[])).most_common(10),top_hep_pts=Counter(pt for p in d.primaryid for pt in pts.get(p,[]) if pt.lower() in hep).most_common(8))
out["coded_profile"]=prof(t); out["coded_by_drug"]=t.target_drug.value_counts().to_dict()
out["brand_soto_profile"]=prof(b_)
# duplicate candidates: same sex, country, event date, age within 1 yr, >=1 shared PT
bk={}
for r in b_.itertuples():
    if r.event_dt and r.sex in ("F","M") and r.country: bk.setdefault((r.sex,r.country,str(r.event_dt)),[]).append(r)
dups=[]
for r in t.itertuples():
    if not (r.event_dt and r.sex in ("F","M") and r.country): continue
    for q in bk.get((r.sex,r.country,str(r.event_dt)),[]):
        a1,a2=pd.to_numeric(r.age_yr,errors="coerce"),pd.to_numeric(q.age_yr,errors="coerce")
        agem= (np.isnan(a1) or np.isnan(a2) or abs(a1-a2)<=1)
        shared=pts.get(r.primaryid,set())&pts.get(q.primaryid,set())
        if agem and shared: dups.append(dict(coded=r.primaryid,brand=q.primaryid,coded_caseid=r.caseid,brand_caseid=q.caseid,hep_coded=r.hepato,hep_brand=q.hepato,n_shared=len(shared),age=[r.age_yr,q.age_yr]))
out["coded_with_complete_keys"]=int(sum(1 for r in t.itertuples() if r.event_dt and r.sex in ("F","M") and r.country))
out["duplicate_candidates"]=dict(pairs=len(dups),coded_reports=len(set(d["coded"] for d in dups)),coded_hepato=len(set(d["coded"] for d in dups if d["hep_coded"])),examples=dups[:5])
pd.DataFrame(dups).to_csv(os.path.join(R,"coded_duplicate_candidates.csv"),index=False)
# loose check: any sotorasib brand report with same sex+country+age (+-1) and >=2 shared PTs, ignoring date
bk2={}
for r in b_.itertuples():
    if r.sex in ("F","M") and r.country and str(r.age_yr)!="": bk2.setdefault((r.sex,r.country),[]).append(r)
loose=set()
for r in t.itertuples():
    if not (r.sex in ("F","M") and r.country and str(r.age_yr)!=""): continue
    for q in bk2.get((r.sex,r.country),[]):
        if abs(float(r.age_yr)-float(q.age_yr))<=1 and len(pts.get(r.primaryid,set())&pts.get(q.primaryid,set()))>=2: loose.add(r.primaryid)
out["loose_duplicate_coded_reports"]=len(loose)
# ---- G 2023+ row and MRTX count
out["coded_received_2023_on"]=int((t.year.astype(int)>=2023).sum())
codes=[]
for fn in glob.glob(os.path.join(RS,"*_codes.tsv")):
    for line in open(fn):
        f=line.rstrip("\n").split("\t")
        if len(f)>=3 and "MRTX" in (" ".join(f[2:])).upper(): codes.append(f[0])
out["mrtx_rows_db"]=len(codes); out["mrtx_reports_kept"]=len(set(codes)&kept); out["mrtx_in_adagrasib_cohort"]=len(set(codes)&set(a.primaryid))
json.dump(out,open(os.path.join(R,"round2_analyses.json"),"w"),indent=1,default=str)
print(json.dumps(out,indent=1,default=str))
