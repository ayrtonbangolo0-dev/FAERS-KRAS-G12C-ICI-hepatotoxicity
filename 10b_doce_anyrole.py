"""Step 10b (Reviewer 2): covariates for docetaxel reports in ANY role (best role per report), one quarter per call."""
import zipfile, glob, os, re, io, sys
DL=os.path.expanduser("~/mnt/Downloads"); OUT=os.path.expanduser("~/fw/doce_any"); os.makedirs(OUT,exist_ok=True)
RK={"PS":0,"SS":1,"I":2,"C":3}
def rows(zf, prefix):
    n=[x for x in zf.namelist() if re.search(r'(^|/)'+prefix+r'\d\dQ\d\.txt$', x, re.I)][0]
    f=io.TextIOWrapper(zf.open(n), encoding="latin-1", newline="")
    h=[x.strip().lower() for x in f.readline().rstrip("\r\n").split("$")]
    for line in f:
        p=line.rstrip("\r\n").split("$"); p+=[""]*(len(h)-len(p)); yield dict(zip(h,p))
q=sys.argv[1].upper(); z=[x for x in glob.glob(os.path.join(DL,"faers_ascii_*.zip")) if q.lower() in x.lower()][0]
zf=zipfile.ZipFile(z); best={}
for r in rows(zf,"DRUG"):
    t=(r.get("drugname","")+" "+r.get("prod_ai","")).upper()
    if "DOCETAXEL" in t or "TAXOTERE" in t:
        k=(RK.get(r.get("role_cod","").upper(),4),r.get("drug_seq",""))
        p=r["primaryid"]
        if p not in best or k<best[p]: best[p]=k
with open(f"{OUT}/{q}_role.tsv","w") as o:
    for p,(k,s) in best.items(): o.write(f"{p}\t{k}\t{s}\n")
with open(f"{OUT}/{q}_indi.tsv","w") as o:
    for r in rows(zf,"INDI"):
        p=r["primaryid"]
        if p in best and r.get("indi_drug_seq","")==best[p][1]: o.write(f'{p}\t{r.get("indi_pt","")}\n')
with open(f"{OUT}/{q}_demo.tsv","w") as o:
    for r in rows(zf,"DEMO"):
        p=r["primaryid"]
        if p in best: o.write("\t".join([p,r.get("age",""),r.get("age_cod",""),r.get("sex",""),r.get("occr_country",""),r.get("rept_cod",""),r.get("occp_cod",""),q])+"\n")
open(f"{OUT}/{q}.done","w").write("ok"); print("done",q,len(best))
