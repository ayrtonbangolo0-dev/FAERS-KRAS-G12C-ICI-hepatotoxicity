"""Step 13: journal Reviewer 1: treatment intervals vs event date, adagrasib case table, unrounded Omega, PT list, Table 5 denominators,
interaction ORs expressed as adagrasib/sotorasib."""
import pandas as pd, numpy as np, os, math, glob, json
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R=os.path.join(B,"results"); W=os.path.expanduser("~/fw"); M=os.path.join(W,"v2","m"); RS=os.path.join(W,"rescan")
out={}
c=pd.read_csv(os.path.join(R,"cohort_reports_v2_flags.csv"),dtype={"primaryid":str},keep_default_na=False)
for k in ["ici_exposed","hepato","outc_LT","outc_HO","outc_DE"]: c[k]=c[k].astype(str)=="True"
drug=pd.read_csv(os.path.join(M,"drug.csv"),dtype=str,keep_default_na=False)
ther=pd.read_csv(os.path.join(M,"ther.csv"),dtype=str,keep_default_na=False)
reac=pd.read_csv(os.path.join(M,"reac.csv"),dtype=str,keep_default_na=False)
hep=[l.strip() for l in open(os.path.join(B,"scripts","hepato_pts.txt")) if l.strip()]; hepl={h.lower() for h in hep}
NARROW=["Hepatotoxicity","Drug-induced liver injury","Liver injury","Hepatitis","Hepatitis acute","Hepatitis fulminant","Hepatitis cholestatic","Immune-mediated hepatitis","Immune-mediated hepatic disorder","Autoimmune hepatitis","Hepatic failure","Acute hepatic failure","Acute on chronic liver failure","Hepatocellular injury","Cholestatic liver injury","Mixed liver injury","Hepatic cytolysis","Hepatic necrosis","Jaundice","Immune-mediated cholangitis"]
nl={x.lower() for x in NARROW}
def dt(s):
    s=str(s).strip()
    try: return pd.Timestamp(s[:4]+"-"+s[4:6]+"-"+s[6:8]) if len(s)==8 else None
    except Exception: return None
cohort=set(c.primaryid)
d=drug[drug.primaryid.isin(cohort)]
th=ther[ther.primaryid.isin(cohort)]
def intervals(p,seqs):
    t=th[(th.primaryid==p)&th.dsg_drug_seq.isin(seqs)]
    return [(dt(a),dt(b)) for a,b in zip(t.start_dt,t.end_dt)]
rows=[]
for r in c[c.ici_exposed].itertuples():
    p=r.primaryid; dd=d[d.primaryid==p]
    kseq=list(dd[dd.target_drug==r.target_drug].drug_seq); iseq=dd[dd.ici_agent!=""]
    ki=intervals(p,kseq); ii=intervals(p,list(iseq.drug_seq))
    ks=min([a for a,b in ki if a is not None],default=None); ke=max([b for a,b in ki if b is not None],default=None)
    is_=min([a for a,b in ii if a is not None],default=None); ie=max([b for a,b in ii if b is not None],default=None)
    ev=dt(r.event_dt)
    # classification against event date
    if ks is None or ev is None or (is_ is None and ie is None): cls="insufficient dates"
    elif ie is not None and ie < ks: cls="sequential: ICI ended before KRAS G12C inhibitor start"
    elif is_ is not None and is_ > ev: cls="ICI started after the event"
    elif is_ is not None and is_ <= ev and ((ie is not None and ie >= ks) or ie is None) and ks <= ev:
        cls="documented overlap before the event" if (ie is not None) else "ICI ongoing (no stop date) at KRAS G12C inhibitor start"
    elif is_ is None and ie is not None and ie >= ks: cls="ICI stop on/after KRAS start, ICI start missing"
    else: cls="other"
    roles=";".join(f"{a}:{b}" for a,b in zip(iseq.ici_agent,iseq.role_cod))
    pts=sorted(set(reac[reac.primaryid==p].pt))
    rows.append(dict(primaryid=p,drug=r.target_drug,hepato=r.hepato,ici=roles,ici_start=is_,ici_end=ie,kras_start=ks,event=ev,
        ici_end_to_kras_start=(ks-ie).days if (ks is not None and ie is not None) else None,
        kras_start_to_event=(ev-ks).days if (ks is not None and ev is not None) else None,cls=cls,
        hep_pts=[x for x in pts if x.lower() in hepl],narrow=any(x.lower() in nl for x in pts),n_pts=len(pts),
        country=r.country,rept=r.rept_cod,year=str(r.quarter)[:4],LT=r.outc_LT,HO=r.outc_HO,DE=r.outc_DE,tumor=r.tumor))
iv=pd.DataFrame(rows); iv.to_csv(os.path.join(R,"ici_interval_vs_event.csv"),index=False)
out["interval_classes_all"]={" | ".join(k):int(v) for k,v in iv.groupby(["drug","cls"]).size().items()}
out["interval_classes_hep"]={" | ".join(k):int(v) for k,v in iv[iv.hepato].groupby(["drug","cls"]).size().items()}
out["ada_hep_cases"]=iv[(iv.drug=="adagrasib")&iv.hepato].to_dict("records")
out["soto_hep_n_multi_pt"]=int((iv[(iv.drug=="sotorasib")&iv.hepato].n_pts>1).sum())
# Table 5 denominators and multiple events
for dname in ["sotorasib","adagrasib"]:
    for lab,m in [("all",c.target_drug==dname),("noici",(c.target_drug==dname)&~c.ici_exposed),("ici",(c.target_drug==dname)&c.ici_exposed)]:
        h=c[m&c.hepato]; t=pd.to_numeric(h.tto_days,errors="coerce")
        out[f"tto_{dname}_{lab}"]=dict(evaluable=int(t.notna().sum()),hep=len(h),within60=int((t<=60).sum()))
hepc=c[c.hepato]; npt=reac[reac.primaryid.isin(hepc.primaryid)].groupby("primaryid").pt.nunique()
out["hep_reports_with_other_reactions_pct"]=round(100*float((npt>1).mean()),1)
# unrounded Omega (re-implemented from step 8/11)
kept=set(open(os.path.join(W,"kept_primaryids.txt")).read().split())
heprep=set()
with open(os.path.join(W,"reac_kept.tsv")) as f:
    for line in f:
        p,_,pt=line.rstrip("\n").partition("\t")
        if pt.lower() in hepl: heprep.add(p)
ici_any=set()
for fn in glob.glob(os.path.join(RS,"*_ici.tsv")):
    for line in open(fn):
        p=line.split("\t")[0]
        if p in kept: ici_any.add(p)
N=len(kept)
def omega(D,excl=set()):
    Nn=N-len(excl); both=D&ici_any; donly=D-ici_any; ionly=ici_any-D-excl
    n11,n10,n01=len(both),len(donly),len(ionly); n00=Nn-n11-n10-n01
    e11,e10,e01=len(both&heprep),len(donly&heprep),len(ionly&heprep); e00=len(heprep-excl)-e11-e10-e01
    o=lambda f:f/(1-f); f00,f10,f01=e00/n00,e10/n10,e01/n01
    g=1-1/(max(o(f00),o(f10))+max(o(f00),o(f01))-o(f00)+1); E=g*n11; om=math.log2((e11+0.5)/(E+0.5)); lo=om-3.3*(e11+0.5)**-0.5-2*(e11+0.5)**-1.5
    return dict(obs=e11,n=n11,E=round(E,3),omega=round(om,4),omega025=round(lo,4))
ada_any=set(drug[drug.target_drug=="adagrasib"].primaryid)&kept; ada_ps=set(c[c.target_drug=="adagrasib"].primaryid)
soto_any=set(drug[drug.target_drug=="sotorasib"].primaryid)&kept; soto_ps=set(c[c.target_drug=="sotorasib"].primaryid)
out["omega_ada_any"]=omega(ada_any); out["omega_ada_ps"]=omega(ada_ps,ada_any-ada_ps)
out["omega_soto_any"]=omega(soto_any); out["omega_soto_ps"]=omega(soto_ps,soto_any-soto_ps)
# interaction expressed as adagrasib/sotorasib
adj=pd.read_csv(os.path.join(R,"adjusted_logistic_model.csv")); r=adj[adj.term=="soto_x_ici"].iloc[0]
out["adjusted_interaction_ada_over_soto"]=[round(1/r.OR,2),round(1/r.hi,2),round(1/r.lo,2),round(r.p,4)]
out["firth_interaction_ada_over_soto_from_rounded"]=[round(1/0.17,2),round(1/0.54,2),round(1/0.06,2)]
json.dump(out,open(os.path.join(R,"reviewer1_journal.json"),"w"),indent=1,default=str)
print(json.dumps(out,indent=1,default=str))
