"""Step 1: stream 22 FAERS ASCII quarters (2021Q2-2026Q3), extract KRAS G12C inhibitor
reports + background tables needed for deduplication and disproportionality."""
import zipfile, glob, os, re, io, csv, sys, time
DL = os.path.expanduser("~/mnt/Downloads")
OUT = os.path.join(DL, "FAERS_KRAS_G12C")
W, X = os.path.join(OUT, "work"), os.path.join(OUT, "extracted")
TARGET = {"sotorasib": ["SOTORASIB", "LUMAKRAS", "LUMYKRAS", "AMG-510", "AMG 510", "AMG510"], "adagrasib": ["ADAGRASIB", "KRAZATI", "MRTX849", "MRTX 849", "MRTX-849"]}
ICI = {"pembrolizumab":["PEMBROLIZUMAB","KEYTRUDA"],"nivolumab":["NIVOLUMAB","OPDIVO","OPDUALAG"],
 "atezolizumab":["ATEZOLIZUMAB","TECENTRIQ"],"durvalumab":["DURVALUMAB","IMFINZI"],
 "cemiplimab":["CEMIPLIMAB","LIBTAYO"],"avelumab":["AVELUMAB","BAVENCIO"],"ipilimumab":["IPILIMUMAB","YERVOY"],
 "tremelimumab":["TREMELIMUMAB","IMJUDO"],"dostarlimab":["DOSTARLIMAB","JEMPERLI"],"tislelizumab":["TISLELIZUMAB","TEVIMBRA"],
 "toripalimab":["TORIPALIMAB","LOQTORZI"],"retifanlimab":["RETIFANLIMAB","ZYNYZ"],"sintilimab":["SINTILIMAB"],
 "camrelizumab":["CAMRELIZUMAB"],"cosibelimab":["COSIBELIMAB","UNLOXCYT"],"serplulimab":["SERPLULIMAB"],"penpulimab":["PENPULIMAB"]}
def classify(text, d):
    t = text.upper()
    for k, names in d.items():
        if any(n in t for n in names): return k
    return None
def rows(zf, prefix):
    name = [n for n in zf.namelist() if re.search(r'(^|/)'+prefix+r'\d\dQ\d\.txt$', n, re.I)]
    if not name: raise SystemExit(f"missing {prefix} in {zf.filename}")
    with zf.open(name[0]) as fh:
        f = io.TextIOWrapper(fh, encoding="latin-1", newline="")
        hdr = [h.strip().lower() for h in f.readline().rstrip("\r\n").split("$")]
        for line in f:
            p = line.rstrip("\r\n").split("$")
            if len(p) < len(hdr): p += [""]*(len(hdr)-len(p))
            yield dict(zip(hdr, p))
zips = sorted(glob.glob(os.path.join(DL, "faers_ascii_*.zip")), key=lambda s: s.lower())
zips = [z for z in zips if re.search(r'(2021q[234]|202[2-6]q\d)', z, re.I)]
if len(sys.argv)>1: zips=[z for z in zips if sys.argv[1].lower() in z.lower()]
print(len(zips), "quarters", flush=True)
QQ=sys.argv[1].upper() if len(sys.argv)>1 else "ALL"; W=os.environ.get("OUT_W",os.path.join(W,"q")); X=os.environ.get("OUT_X",os.path.join(X,"q"))
demo_all = open(os.path.join(W, QQ+"_demo.tsv"), "w")      # primaryid caseid fda_dt quarter
reac_all = open(os.path.join(W, QQ+"_reac.tsv"), "w")      # primaryid pt
deleted = open(os.path.join(W, QQ+"_deleted.txt"), "w")
outs = {}
def out(tbl, fields):
    if tbl not in outs:
        fh = open(os.path.join(X, QQ+"_"+tbl + ".csv"), "w", newline="")
        wr = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore"); wr.writeheader(); outs[tbl] = (fh, wr, fields)
    return outs[tbl][1]
for z in zips:
    q = re.search(r'(\d{4}q\d)', z, re.I).group(1).upper(); t0 = time.time()
    zf = zipfile.ZipFile(z)
    for n in zf.namelist():
        if re.search(r'delete.*\.txt$', n, re.I):
            for line in io.TextIOWrapper(zf.open(n), encoding="latin-1"):
                c = line.strip().split("$")[0]
                if c.isdigit(): deleted.write(c + "\n")
    # DRUG pass: find target reports (any role) and keep all drug rows for them
    tgt = {}
    for r in rows(zf, "DRUG"):
        k = classify(r.get("drugname","") + " " + r.get("prod_ai",""), TARGET)
        if k: tgt.setdefault(r["primaryid"], set()).add((k, r.get("role_cod","")))
    ids = set(tgt)
    dfields = ["quarter","primaryid","caseid","drug_seq","role_cod","drugname","prod_ai","route","dose_vbm","dechal","rechal","target_drug","ici_agent"]
    for r in rows(zf, "DRUG"):
        if r["primaryid"] in ids:
            txt = r.get("drugname","") + " " + r.get("prod_ai","")
            r.update(quarter=q, target_drug=classify(txt, TARGET) or "", ici_agent=classify(txt, ICI) or "")
            out("drug", dfields).writerow(r)
    for r in rows(zf, "DEMO"):
        if not os.environ.get("SKIP_ALL"): demo_all.write(f'{r["primaryid"]}\t{r["caseid"]}\t{r.get("fda_dt","")}\t{q}\n')
        if r["primaryid"] in ids:
            r["quarter"] = q
            out("demo", ["quarter","primaryid","caseid","caseversion","i_f_code","event_dt","init_fda_dt","fda_dt","rept_cod","mfr_sndr","age","age_cod","age_grp","sex","wt","wt_cod","rept_dt","occp_cod","reporter_country","occr_country"]).writerow(r)
    for r in rows(zf, "REAC"):
        if not os.environ.get("SKIP_ALL"): reac_all.write(f'{r["primaryid"]}\t{r["pt"].strip()}\n')
        if r["primaryid"] in ids:
            r["quarter"] = q; out("reac", ["quarter","primaryid","caseid","pt","drug_rec_act"]).writerow(r)
    for tbl, pre, f in [("outc","OUTC",["quarter","primaryid","caseid","outc_cod"]),
                        ("indi","INDI",["quarter","primaryid","caseid","indi_drug_seq","indi_pt"]),
                        ("ther","THER",["quarter","primaryid","caseid","dsg_drug_seq","start_dt","end_dt","dur","dur_cod"])]:
        for r in rows(zf, pre):
            if r["primaryid"] in ids:
                if "outc_code" in r and "outc_cod" not in r: r["outc_cod"] = r["outc_code"]
                r["quarter"] = q; out(tbl, f).writerow(r)
    print(f"{q}: {len(ids)} KRAS G12C reports, {time.time()-t0:.0f}s", flush=True)
for fh, _, _ in outs.values(): fh.close()
demo_all.close(); reac_all.close(); deleted.close()
open(os.path.join(W, QQ+".done"),"w").write("ok")
print("DONE", QQ, flush=True)
