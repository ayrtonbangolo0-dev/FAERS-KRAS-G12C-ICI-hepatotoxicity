"""Step 7: database-wide flags needed for interaction and reference-group analyses (one quarter per call)."""
import zipfile, glob, os, re, io, sys
DL=os.path.expanduser("~/mnt/Downloads"); OUT=os.path.expanduser("~/fw/rescan")
ICI={"pembrolizumab":["PEMBROLIZUMAB","KEYTRUDA"],"nivolumab":["NIVOLUMAB","OPDIVO","OPDUALAG"],"atezolizumab":["ATEZOLIZUMAB","TECENTRIQ"],
 "durvalumab":["DURVALUMAB","IMFINZI"],"cemiplimab":["CEMIPLIMAB","LIBTAYO"],"avelumab":["AVELUMAB","BAVENCIO"],"ipilimumab":["IPILIMUMAB","YERVOY"],
 "tremelimumab":["TREMELIMUMAB","IMJUDO"],"dostarlimab":["DOSTARLIMAB","JEMPERLI"],"tislelizumab":["TISLELIZUMAB","TEVIMBRA"],"toripalimab":["TORIPALIMAB","LOQTORZI"],
 "retifanlimab":["RETIFANLIMAB","ZYNYZ"],"sintilimab":["SINTILIMAB"],"camrelizumab":["CAMRELIZUMAB"],"cosibelimab":["COSIBELIMAB","UNLOXCYT"],"serplulimab":["SERPLULIMAB"],"penpulimab":["PENPULIMAB"]}
CODES=re.compile(r"AMG[\s\-]?510|MRTX[\s\-]?849|FULZERASIB|GARSORASIB|GLECIRASIB|SOSIMERASIB|DIVARASIB|OLOMORASIB|D-1553|GFH925|JDQ443|OPNURASIB")
def cls(t):
    for k,v in ICI.items():
        if any(n in t for n in v): return k
def rows(zf, prefix):
    n=[x for x in zf.namelist() if re.search(r'(^|/)'+prefix+r'\d\dQ\d\.txt$', x, re.I)][0]
    f=io.TextIOWrapper(zf.open(n), encoding="latin-1", newline="")
    h=[x.strip().lower() for x in f.readline().rstrip("\r\n").split("$")]
    for line in f:
        p=line.rstrip("\r\n").split("$"); p+=[""]*(len(h)-len(p)); yield dict(zip(h,p))
q=sys.argv[1].upper(); z=[x for x in glob.glob(os.path.join(DL,"faers_ascii_*.zip")) if q.lower() in x.lower()][0]
zf=zipfile.ZipFile(z)
target=set(l.split(",")[1] for l in open(os.path.expanduser("~/fw/x/drug.csv")).read().splitlines()[1:])
ici=open(f"{OUT}/{q}_ici.tsv","w"); doc=open(f"{OUT}/{q}_doce.tsv","w"); cod=open(f"{OUT}/{q}_codes.tsv","w")
for r in rows(zf,"DRUG"):
    t=(r.get("drugname","")+" "+r.get("prod_ai","")).upper(); role=r.get("role_cod","").upper()
    k=cls(t)
    if k: ici.write(f'{r["primaryid"]}\t{k}\t{role}\n')
    if "DOCETAXEL" in t or "TAXOTERE" in t: doc.write(f'{r["primaryid"]}\t{role}\n')
    if CODES.search(t): cod.write(f'{r["primaryid"]}\t{role}\t{r.get("drugname","")}\t{r.get("prod_ai","")}\n')
rp=open(f"{OUT}/{q}_rpsr.tsv","w")
for r in rows(zf,"RPSR"):
    if r["primaryid"] in target: rp.write(f'{r["primaryid"]}\t{r.get("rpsr_cod","")}\n')
for f in (ici,doc,cod,rp): f.close()
open(f"{OUT}/{q}.done","w").write("ok"); print("done",q)
