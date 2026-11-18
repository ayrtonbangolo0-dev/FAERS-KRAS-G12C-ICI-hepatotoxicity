"""Step 3: KRAS G12C inhibitor hepatotoxicity by co-reported ICI exposure.
Primary cohort: deduplicated reports with sotorasib/adagrasib as PRIMARY SUSPECT (role PS)."""
import pandas as pd, numpy as np, os, math, json
B = os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); X, W, R = [os.path.join(B, d) for d in ("extracted","work","results")]; R = os.environ.get("FAERS_R", R)
X = os.environ.get("FAERS_X", X); W = os.environ.get("FAERS_W", W)
os.makedirs(R, exist_ok=True)
rd = lambda t: pd.read_csv(os.path.join(X, t+".csv"), dtype=str, keep_default_na=False)
drug, demo, reac, outc, indi, ther = [rd(t) for t in ("drug","demo","reac","outc","indi","ther")]
kept = set(open(os.path.join(W,"kept_primaryids.txt")).read().split())
hep = {l.strip().lower() for l in open(os.path.join(B,"scripts","hepato_pts.txt")) if l.strip()}
N = len(kept); N_hep = int(open(os.path.join(W,"n_hepato_reports.txt")).read())
ptc = pd.read_csv(os.path.join(W,"pt_counts.tsv"), sep="\t", header=None, names=["pt","n"], quoting=3, keep_default_na=False)
log = {"raw_reports": int(open(os.path.join(W,"n_raw.txt")).read()), "dedup_reports_database": N, "database_hepato_reports": N_hep}
drug = drug[drug.primaryid.isin(kept)]
ps = drug[(drug.target_drug!="") & (drug.role_cod.str.upper()=="PS")]
both = ps.groupby("primaryid").target_drug.nunique(); ps = ps[ps.primaryid.isin(both[both==1].index)]
coh = ps.drop_duplicates("primaryid")[["primaryid","caseid","target_drug","drug_seq"]].copy()
ici = drug[drug.ici_agent!=""].groupby("primaryid").ici_agent.apply(lambda s: ";".join(sorted(set(s))))
coh["ici"] = coh.primaryid.map(ici).fillna(""); coh["ici_exposed"] = coh.ici!=""
r = reac[reac.primaryid.isin(coh.primaryid)].copy(); r["is_hep"] = r.pt.str.strip().str.lower().isin(hep)
coh["hepato"] = coh.primaryid.isin(r[r.is_hep].primaryid)
coh["severe_hep"] = coh.primaryid.isin(r[r.pt.str.lower().isin({"hepatic failure","acute hepatic failure","hepatitis fulminant","hepatic necrosis","acute on chronic liver failure"})].primaryid)
o = outc[outc.primaryid.isin(coh.primaryid)].groupby("primaryid").outc_cod.apply(lambda s: set(s.str.upper()))
for c in ["DE","LT","HO","DS"]: coh["outc_"+c] = coh.primaryid.map(lambda p: c in o.get(p, set()))
coh["serious"] = coh[["outc_DE","outc_LT","outc_HO","outc_DS"]].any(axis=1)
d = demo.drop_duplicates("primaryid").set_index("primaryid")
coh["sex"] = coh.primaryid.map(d.sex); coh["country"] = coh.primaryid.map(d.occr_country); coh["event_dt"] = coh.primaryid.map(d.event_dt)
coh["rept_cod"] = coh.primaryid.map(d.rept_cod)
def age(p):
    try:
        a, u = float(d.at[p,"age"]), d.at[p,"age_cod"].upper()
        return a if u=="YR" else a/12 if u=="MON" else a/52 if u=="WK" else a/365 if u=="DY" else a*10 if u=="DEC" else np.nan
    except Exception: return np.nan
coh["age_yr"] = coh.primaryid.map(age)
ind = indi.merge(coh[["primaryid","drug_seq"]], left_on=["primaryid","indi_drug_seq"], right_on=["primaryid","drug_seq"])
ind_pt = ind.groupby("primaryid").indi_pt.apply(lambda s: "; ".join(sorted(set(s))))
coh["indication"] = coh.primaryid.map(ind_pt).fillna("")
il = coh.indication.str.lower()
coh["tumor"] = np.select([il.str.contains("colorectal|colon|rectal|rectum|caecum|cecum|sigmoid"), il.str.contains("lung|nsclc|bronch"), il.str.contains("pancrea"), il==""],["CRC","NSCLC","Pancreatic","Unknown"], "Other")
def dt(s):
    s = str(s).strip()
    try: return pd.Timestamp(s[:4]+"-"+s[4:6]+"-"+s[6:8]) if len(s)==8 else None
    except Exception: return None
th = ther.merge(coh[["primaryid","drug_seq"]], left_on=["primaryid","dsg_drug_seq"], right_on=["primaryid","drug_seq"])
start = th.assign(sd=th.start_dt.map(dt)).dropna(subset=["sd"]).groupby("primaryid").sd.min()
ev = coh.set_index("primaryid").event_dt.map(dt)
tto = (ev - coh.set_index("primaryid").index.map(start).to_series(index=coh.primaryid.values)).dt.days
coh["tto_days"] = coh.primaryid.map(tto)
coh.loc[(coh.tto_days<0)|(coh.tto_days>3650),"tto_days"] = np.nan
# ICI -> KRAS interval when both start/end dates exist
icirows = drug[(drug.ici_agent!="") & drug.primaryid.isin(coh.primaryid)][["primaryid","drug_seq"]]
it = ther.merge(icirows, left_on=["primaryid","dsg_drug_seq"], right_on=["primaryid","drug_seq"])
it["last"] = [dt(e) or dt(s) for s, e in zip(it.start_dt, it.end_dt)]
ici_last = it.dropna(subset=["last"]).groupby("primaryid")["last"].max()
coh["ici_to_kras_days"] = coh.primaryid.map(lambda p: (start[p]-ici_last[p]).days if p in start.index and p in ici_last.index else np.nan)
coh.to_csv(os.path.join(R,"cohort_reports.csv"), index=False)
# ---- disproportionality vs whole database
def stats(a, b, c, d_):
    a,b,c,d_ = [x+0.5 if 0 in (a,b,c,d_) else x for x in (a,b,c,d_)]
    ror = (a/b)/(c/d_); se = math.sqrt(1/a+1/b+1/c+1/d_)
    prr = (a/(a+b))/(c/(c+d_)); sep = math.sqrt(1/a-1/(a+b)+1/c-1/(c+d_))
    n = a+b+c+d_; E = (a+b)*(a+c)/n; ic = math.log2((a+0.5)/(E+0.5))
    chi = n*(abs(a*d_-b*c)-n/2)**2/((a+b)*(c+d_)*(a+c)*(b+d_))
    return dict(ROR=ror, ROR_lo=ror*math.exp(-1.96*se), ROR_hi=ror*math.exp(1.96*se), PRR=prr, PRR_lo=prr*math.exp(-1.96*sep), PRR_hi=prr*math.exp(1.96*sep), chi2=chi, IC=ic, IC025=ic-3.3*(a+0.5)**-0.5-2*(a+0.5)**-1.5)
rows = []
for drugname in ["sotorasib","adagrasib","any"]:
    for grp in ["all","ICI+","ICI-"]:
        s = coh if drugname=="any" else coh[coh.target_drug==drugname]
        s = s if grp=="all" else s[s.ici_exposed==(grp=="ICI+")]
        a = int(s.hepato.sum()); n1 = len(s); b = n1-a
        c = N_hep - int(coh.hepato.sum()) if True else 0; d_ = (N - len(coh)) - c
        rows.append(dict(drug=drugname, ici_group=grp, n_reports=n1, n_hepato=a, pct_hepato=round(100*a/n1,1) if n1 else np.nan, **(stats(a,b,c,d_) if n1 else {})))
disp = pd.DataFrame(rows); disp.to_csv(os.path.join(R,"table2_disproportionality_vs_database.csv"), index=False)
# ---- within-cohort comparisons (reporting odds of hepatotoxicity)
def within(s, mask, label):
    a=int(s[mask].hepato.sum()); b=int(mask.sum())-a; c=int(s[~mask].hepato.sum()); d_=int((~mask).sum())-c
    st = stats(a,b,c,d_); return dict(comparison=label, exposed_hep=a, exposed_n=a+b, unexposed_hep=c, unexposed_n=c+d_, ROR=st["ROR"], lo=st["ROR_lo"], hi=st["ROR_hi"])
W2 = [within(coh, coh.ici_exposed, "ICI+ vs ICI- (both drugs)")]
for dname in ["sotorasib","adagrasib"]:
    s = coh[coh.target_drug==dname]
    if len(s) and s.ici_exposed.any(): W2.append(within(s, s.ici_exposed, f"ICI+ vs ICI- ({dname})"))
for grp in [True, False]:
    s = coh[coh.ici_exposed==grp]
    if (s.target_drug=="sotorasib").any() and (s.target_drug=="adagrasib").any():
        W2.append(within(s, s.target_drug=="sotorasib", f"sotorasib vs adagrasib ({'ICI+' if grp else 'ICI-'})"))
pd.DataFrame(W2).to_csv(os.path.join(R,"table3_within_cohort_comparisons.csv"), index=False)
# ---- table 1 characteristics
t1 = []
for (dn, g), s in coh.groupby(["target_drug","ici_exposed"]):
    t1.append(dict(drug=dn, ici_exposed=g, n=len(s), female_pct=round(100*(s.sex=="F").mean(),1), age_median=s.age_yr.median(), age_known=s.age_yr.notna().sum(),
      NSCLC=(s.tumor=="NSCLC").sum(), CRC=(s.tumor=="CRC").sum(), pancreatic=(s.tumor=="Pancreatic").sum(), other=(s.tumor=="Other").sum(), unknown=(s.tumor=="Unknown").sum(),
      US_pct=round(100*(s.country=="US").mean(),1), hepato=s.hepato.sum(), severe_hep=s.severe_hep.sum(), death=s.outc_DE.sum(), hosp=s.outc_HO.sum()))
pd.DataFrame(t1).to_csv(os.path.join(R,"table1_characteristics.csv"), index=False)
# ---- time to onset among hepatotoxicity reports
tt = coh[coh.hepato & coh.tto_days.notna()].groupby(["target_drug","ici_exposed"]).tto_days.describe()
tt.to_csv(os.path.join(R,"table4_time_to_onset_hepato.csv"))
# ---- CRC subgroup
crc = coh[coh.tumor=="CRC"]
crc.groupby(["target_drug","ici_exposed"]).agg(n=("primaryid","size"), hepato=("hepato","sum")).to_csv(os.path.join(R,"table5_crc_subgroup.csv"))
# ---- top hepatic PTs
r2 = r[r.is_hep & r.primaryid.isin(coh.primaryid)].merge(coh[["primaryid","target_drug","ici_exposed"]])
r2.groupby(["target_drug","ici_exposed","pt"]).primaryid.nunique().rename("n").reset_index().sort_values("n", ascending=False).to_csv(os.path.join(R,"table6_hepatic_pts.csv"), index=False)
log.update(cohort_PS=len(coh), sotorasib=int((coh.target_drug=="sotorasib").sum()), adagrasib=int((coh.target_drug=="adagrasib").sum()),
           ici_exposed=int(coh.ici_exposed.sum()), crc=int((coh.tumor=="CRC").sum()), with_tto=int(coh.tto_days.notna().sum()),
           with_ici_interval=int(coh.ici_to_kras_days.notna().sum()))
json.dump(log, open(os.path.join(R,"flow_counts.json"),"w"), indent=2)
print(json.dumps(log, indent=2)); print(disp.round(2).to_string()); print(pd.DataFrame(W2).round(2).to_string())
print(pd.read_csv(os.path.join(R,"table1_characteristics.csv")).to_string()); print(tt.round(1).to_string())
print(pd.read_csv(os.path.join(R,"table5_crc_subgroup.csv")).to_string())
