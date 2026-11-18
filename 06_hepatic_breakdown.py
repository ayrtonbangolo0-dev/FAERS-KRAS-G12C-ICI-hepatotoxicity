"""Step 6: breakdown of hepatotoxicity by preferred term and clinical category, with PT-level ROR vs database."""
import pandas as pd, numpy as np, os, math
B = os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R = os.path.join(B,"results")
W = os.environ.get("FAERS_W", os.path.join(B,"work")); X = os.environ.get("FAERS_X", os.path.join(B,"extracted"))
c = pd.read_csv(os.path.join(R,"cohort_reports.csv"), dtype={"primaryid":str}, keep_default_na=False)
c["ici_exposed"] = c.ici_exposed.astype(str)=="True"
reac = pd.read_csv(os.path.join(X,"reac.csv"), dtype=str, keep_default_na=False); reac["ptl"]=reac.pt.str.strip().str.lower()
hep = [l.strip() for l in open(os.path.join(B,"scripts","hepato_pts.txt")) if l.strip()]
CAT = {"Liver failure":["Hepatic failure","Acute hepatic failure","Hepatitis fulminant","Hepatic necrosis","Acute on chronic liver failure"],
 "Drug-induced liver injury / hepatotoxicity":["Drug-induced liver injury","Hepatotoxicity","Liver injury","Hepatocellular injury","Mixed liver injury","Hepatic cytolysis"],
 "Hepatitis (any)":["Hepatitis","Hepatitis acute","Hepatitis cholestatic"],
 "Immune-mediated or autoimmune hepatitis":["Immune-mediated hepatitis","Immune-mediated hepatic disorder","Autoimmune hepatitis","Immune-mediated cholangitis"],
 "Cholestasis or jaundice":["Cholestasis","Cholestatic liver injury","Jaundice"],
 "Liver test abnormalities only terms":["Hepatic function abnormal","Liver function test abnormal","Liver function test increased","Hepatic enzyme increased","Hepatic enzyme abnormal",
   "Transaminases increased","Transaminases abnormal","Hypertransaminasaemia","Alanine aminotransferase increased","Aspartate aminotransferase increased",
   "Gamma-glutamyltransferase increased","Blood alkaline phosphatase increased","Blood bilirubin increased","Hyperbilirubinaemia","Liver disorder"]}
assert sorted(sum(CAT.values(),[]))==sorted(hep), set(hep)^set(sum(CAT.values(),[]))
r = reac[reac.primaryid.isin(c.primaryid)]
# database counts (unique reports per PT, deduplicated)
ptc = pd.read_csv(os.path.join(W,"pt_counts.tsv"), sep="\t", header=None, names=["pt","n"], quoting=3, keep_default_na=False); ptc["ptl"]=ptc.pt.str.lower()
dbn = ptc.groupby("ptl").n.sum().to_dict(); N = len(open(os.path.join(W,"kept_primaryids.txt")).read().split())
# category database counts need report-level union: compute from reac_kept
catdb = {k:set() for k in CAT}; lookup = {p.lower():k for k,v in CAT.items() for p in v}
with open(os.path.join(W,"reac_kept.tsv")) as f:
    for line in f:
        p,_,pt = line.rstrip("\n").partition("\t"); k = lookup.get(pt.lower())
        if k: catdb[k].add(p)
catdb = {k:len(v) for k,v in catdb.items()}
groups = {"Sotorasib, no ICI":(c.target_drug=="sotorasib")&~c.ici_exposed, "Sotorasib + ICI":(c.target_drug=="sotorasib")&c.ici_exposed,
          "Adagrasib, no ICI":(c.target_drug=="adagrasib")&~c.ici_exposed, "Adagrasib + ICI":(c.target_drug=="adagrasib")&c.ici_exposed,
          "Sotorasib, all":c.target_drug=="sotorasib", "Adagrasib, all":c.target_drug=="adagrasib"}
def ror(a,n1,cdb,ncoh_with):  # cdb = database reports with event; subtract cohort reports with event
    cc = cdb - ncoh_with; d = (N-len(c)) - cc; b = n1-a
    if a==0: return np.nan,np.nan,np.nan
    rr=(a/b)/(cc/d); se=math.sqrt(1/a+1/b+1/cc+1/d); return rr, rr*math.exp(-1.96*se), rr*math.exp(1.96*se)
rows=[]
units = [("category",k,[p.lower() for p in v]) for k,v in CAT.items()] + [("PT",p,[p.lower()]) for p in hep]
for kind,name,pts in units:
    ids_all = set(r[r.ptl.isin(pts)].primaryid); dbcount = catdb[name] if kind=="category" else dbn.get(pts[0],0)
    row = dict(level=kind, term=name, database_reports=dbcount)
    for g,m in groups.items():
        s = c[m]; a = int(s.primaryid.isin(ids_all).sum()); n1=len(s)
        rr,l,h = ror(a,n1,dbcount,int(c.primaryid.isin(ids_all).sum()))
        row[f"{g} n (%)"] = f"{a} ({100*a/n1:.1f})"; row[f"{g} ROR (95%CI)"] = "-" if np.isnan(rr) else f"{rr:.2f} ({l:.2f}-{h:.2f})"
    rows.append(row)
T = pd.DataFrame(rows); T.to_csv(os.path.join(R,"hepatic_breakdown.csv"), index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",40)
cols = ["term"]+[f"{g} n (%)" for g in list(groups)[:4]]
print(T[T.level=="category"][cols+["Sotorasib, all ROR (95%CI)","Adagrasib, all ROR (95%CI)","Sotorasib + ICI ROR (95%CI)","Adagrasib + ICI ROR (95%CI)"]].to_string(index=False))
P = T[T.level=="PT"].copy(); P["tot"]=P["Sotorasib, all n (%)"].str.split().str[0].astype(int)+P["Adagrasib, all n (%)"].str.split().str[0].astype(int)
print(P.sort_values("tot",ascending=False)[cols+["Sotorasib, all ROR (95%CI)","Adagrasib, all ROR (95%CI)"]].head(20).to_string(index=False))
