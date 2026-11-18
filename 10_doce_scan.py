"""Step 10: covariates for docetaxel primary-suspect reports (one quarter per call)."""
import zipfile, glob, os, re, io, sys
DL=os.path.expanduser("~/mnt/Downloads"); OUT=os.path.expanduser("~/fw/doce")
def rows(zf, prefix):
    n=[x for x in zf.namelist() if re.search(r'(^|/)'+prefix+r'\d\dQ\d\.txt$', x, re.I)][0]
    f=io.TextIOWrapper(zf.open(n), encoding="latin-1", newline="")
    h=[x.strip().lower() for x in f.readline().rstrip("\r\n").split("$")]
    for line in f:
        p=line.rstrip("\r\n").split("$"); p+=[""]*(len(h)-len(p)); yield dict(zip(h,p))
q=sys.argv[1].upper(); z=[x for x in glob.glob(os.path.join(DL,"faers_ascii_*.zip")) if q.lower() in x.lower()][0]
zf=zipfile.ZipFile(z); ps={}
for r in rows(zf,"DRUG"):
    if r.get("role_cod","").upper()!="PS": continue
    t=(r.get("drugname","")+" "+r.get("prod_ai","")).upper()
    if "DOCETAXEL" in t or "TAXOTERE" in t: ps.setdefault(r["primaryid"],set()).add(r.get("drug_seq",""))
with open(f"{OUT}/{q}_indi.tsv","w") as o:
    for r in rows(zf,"INDI"):
        p=r["primaryid"]
        if p in ps and r.get("indi_drug_seq","") in ps[p]: o.write(f'{p}\t{r.get("indi_pt","")}\n')
with open(f"{OUT}/{q}_demo.tsv","w") as o:
    for r in rows(zf,"DEMO"):
        p=r["primaryid"]
        if p in ps: o.write("\t".join([p,r.get("age",""),r.get("age_cod",""),r.get("sex",""),r.get("occr_country",""),r.get("rept_cod",""),r.get("occp_cod",""),q])+"\n")
open(f"{OUT}/{q}.done","w").write("ok"); print("done",q,len(ps))
