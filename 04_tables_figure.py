"""Step 4: manuscript Tables 1-3 (csv + json) and Figure 1 forest plot (300 dpi TIFF/PNG)."""
import pandas as pd, numpy as np, json, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R = os.path.join(B,"results"); M = os.path.join(B,"manuscript"); os.makedirs(M, exist_ok=True)
c = pd.read_csv(os.path.join(R,"cohort_reports.csv"), keep_default_na=False)
for k in ["ici_exposed","hepato","severe_hep","outc_DE","outc_LT","outc_HO","outc_DS","serious"]: c[k] = c[k].astype(str)=="True"
c["age_yr"] = pd.to_numeric(c.age_yr, errors="coerce"); c["tto_days"] = pd.to_numeric(c.tto_days, errors="coerce")
c.loc[c.tumor=="Other","tumor"] = np.where(c.loc[c.tumor=="Other","indication"].str.contains("unknown indication|K-ras gene mutation|Neoplasm malignant|Illness", case=False), "Unknown", "Other")
groups = [("Sotorasib, no ICI", (c.target_drug=="sotorasib")&~c.ici_exposed), ("Sotorasib + ICI", (c.target_drug=="sotorasib")&c.ici_exposed),
          ("Adagrasib, no ICI", (c.target_drug=="adagrasib")&~c.ici_exposed), ("Adagrasib + ICI", (c.target_drug=="adagrasib")&c.ici_exposed), ("All reports", c.target_drug!="")]
def npct(m, d): return f"{int(m.sum()):,} ({100*m.sum()/len(d):.1f})" if len(d) else "-"
t1 = {}
for g, m in groups:
    s = c[m]; a = s.age_yr.dropna(); col = {}
    col["Reports, n"] = f"{len(s):,}"
    col["Age, median (IQR), yr"] = f"{a.median():.0f} ({a.quantile(.25):.0f}–{a.quantile(.75):.0f})" if len(a) else "-"
    col["  Age reported"] = npct(s.age_yr.notna(), s)
    col["Sex: female"] = npct(s.sex=="F", s); col["Sex: male"] = npct(s.sex=="M", s); col["Sex: not reported"] = npct(~s.sex.isin(["F","M"]), s)
    col["Reporter: health professional"] = npct(s.rept_cod.isin(["EXP","PER"]) & False, s)  # placeholder replaced below
    col["Report type: expedited"] = npct(s.rept_cod=="EXP", s); col["Report type: periodic"] = npct(s.rept_cod=="PER", s); col["Report type: direct"] = npct(s.rept_cod=="DIR", s)
    col.pop("Reporter: health professional")
    col["Country: United States"] = npct(s.country=="US", s); col["Country: other"] = npct((s.country!="US")&(s.country!=""), s); col["Country: not reported"] = npct(s.country=="", s)
    for t in ["NSCLC","CRC","Pancreatic","Other","Unknown"]:
        col[{"NSCLC":"Indication: NSCLC","CRC":"Indication: colorectal cancer","Pancreatic":"Indication: pancreatic cancer","Other":"Indication: other cancer","Unknown":"Indication: not specified"}[t]] = npct(s.tumor==t, s)
    col["Outcome: death"] = npct(s.outc_DE, s); col["Outcome: life-threatening"] = npct(s.outc_LT, s); col["Outcome: hospitalization"] = npct(s.outc_HO, s); col["Any serious outcome"] = npct(s.serious, s)
    t1[g] = col
T1 = pd.DataFrame(t1); T1.index.name = "Characteristic"; T1.to_csv(os.path.join(M,"Table1_characteristics.csv"))
ici_agents = c[c.ici_exposed].ici.str.split(";").explode().value_counts().to_dict()
d = pd.read_csv(os.path.join(R,"table2_disproportionality_vs_database.csv"))
lab = {("sotorasib","all"):"Sotorasib, all",("sotorasib","ICI-"):"Sotorasib, no ICI",("sotorasib","ICI+"):"Sotorasib + ICI",
       ("adagrasib","all"):"Adagrasib, all",("adagrasib","ICI-"):"Adagrasib, no ICI",("adagrasib","ICI+"):"Adagrasib + ICI",
       ("any","all"):"Either drug, all",("any","ICI-"):"Either drug, no ICI",("any","ICI+"):"Either drug + ICI"}
d["Group"] = [lab[(a,b)] for a,b in zip(d.drug,d.ici_group)]
T2 = pd.DataFrame({"Group": d.Group, "Reports, n": d.n_reports.map("{:,}".format), "Hepatotoxicity, n (%)": [f"{a} ({p:.1f})" for a,p in zip(d.n_hepato,d.pct_hepato)],
   "ROR (95%CI)": [f"{r:.2f} ({l:.2f}–{h:.2f})" for r,l,h in zip(d.ROR,d.ROR_lo,d.ROR_hi)], "PRR (95%CI)": [f"{r:.2f} ({l:.2f}–{h:.2f})" for r,l,h in zip(d.PRR,d.PRR_lo,d.PRR_hi)],
   "χ²": d.chi2.map("{:.1f}".format), "IC (IC025)": [f"{i:.2f} ({j:.2f})" for i,j in zip(d.IC,d.IC025)], "Signal": ["Yes" if (l>1 and p>=2 and ch>=4 and j>0) else "No" for l,p,ch,j in zip(d.ROR_lo,d.PRR,d.chi2,d.IC025)]})
T2.to_csv(os.path.join(M,"Table2_disproportionality.csv"), index=False)
w = pd.read_csv(os.path.join(R,"table3_within_cohort_comparisons.csv"))
wl = {"ICI+ vs ICI- (both drugs)":"ICI vs no ICI — either drug","ICI+ vs ICI- (sotorasib)":"ICI vs no ICI — sotorasib","ICI+ vs ICI- (adagrasib)":"ICI vs no ICI — adagrasib",
      "sotorasib vs adagrasib (ICI-)":"Sotorasib vs adagrasib — no ICI","sotorasib vs adagrasib (ICI+)":"Sotorasib vs adagrasib — with ICI"}
w["Comparison"] = w.comparison.map(wl); w = w.set_index("Comparison").loc[list(wl.values())].reset_index()
tto = {}
for g, m in groups[:4]:
    s = c[m & c.hepato].tto_days.dropna(); tto[g] = f"{s.median():.0f} ({s.quantile(.25):.0f}–{s.quantile(.75):.0f}); n={len(s)}" if len(s) else "-"
T3 = pd.DataFrame({"Comparison": w.Comparison, "Exposed: hepatotoxicity/reports": [f"{a}/{n} ({100*a/n:.1f}%)" for a,n in zip(w.exposed_hep,w.exposed_n)],
   "Reference: hepatotoxicity/reports": [f"{a}/{n} ({100*a/n:.1f}%)" for a,n in zip(w.unexposed_hep,w.unexposed_n)], "ROR (95%CI)": [f"{r:.2f} ({l:.2f}–{h:.2f})" for r,l,h in zip(w.ROR,w.lo,w.hi)]})
T3.to_csv(os.path.join(M,"Table3_within_cohort.csv"), index=False)
json.dump({"T1": T1.reset_index().to_dict("split"), "T2": T2.to_dict("split"), "T3": T3.to_dict("split"), "tto": tto, "ici_agents": ici_agents,
           "flow": json.load(open(os.path.join(R,"flow_counts.json")))}, open(os.path.join(M,"tables.json"),"w"), ensure_ascii=False, indent=1)
# ---------- Figure 1 forest plot
wf = w.copy()
_ra = os.path.join(R,"reviewer_analyses.json")
if os.path.exists(_ra):
    _d = json.load(open(_ra))["docetaxel_PS_ICI_vs_noICI"]; import re as _re
    _e = _re.match(r"(\d+)/(\d+)", _d["ici"]); _n = _re.match(r"(\d+)/(\d+)", _d["no_ici"])
    wf = pd.concat([wf, pd.DataFrame([{"Comparison":"ICI vs no ICI — docetaxel, all","exposed_hep":int(_e.group(1)),"exposed_n":int(_e.group(2)),"unexposed_hep":int(_n.group(1)),"unexposed_n":int(_n.group(2)),"ROR":_d["ROR"][0],"lo":_d["ROR"][1],"hi":_d["ROR"][2]}])], ignore_index=True)
    _r2 = os.path.join(R,"round2_analyses.json")
    if os.path.exists(_r2):
        _q = json.load(open(_r2))["ici_or_by_drug"]["docetaxel_NSCLC"]
        _e = _re.match(r"(\d+)/(\d+)", _q["ici"]); _n = _re.match(r"(\d+)/(\d+)", _q["noici"])
        wf = pd.concat([wf, pd.DataFrame([{"Comparison":"ICI vs no ICI — docetaxel, NSCLC","exposed_hep":int(_e.group(1)),"exposed_n":int(_e.group(2)),"unexposed_hep":int(_n.group(1)),"unexposed_n":int(_n.group(2)),"ROR":_q["crude"][0],"lo":_q["crude"][1],"hi":_q["crude"][2]}])], ignore_index=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9})
A = d[d.drug!="any"].copy(); order = ["Sotorasib, all","Sotorasib, no ICI","Sotorasib + ICI","Adagrasib, all","Adagrasib, no ICI","Adagrasib + ICI"]
A = A.set_index("Group").loc[order].reset_index()
fig, axes = plt.subplots(2, 1, figsize=(9.6, 7.8), gridspec_kw={"height_ratios":[6,6],"hspace":0.55})
def panel(ax, labels, est, lo, hi, txt, title, xlim, colors, hdr="n/N  ·  ROR (95%CI)"):
    y = np.arange(len(labels))[::-1]
    ax.errorbar(est, y, xerr=[est-lo, hi-est], fmt="none", ecolor="#444", elinewidth=1.1, capsize=2.5)
    ax.scatter(est, y, marker="s", s=38, c=colors, zorder=3, edgecolors="black", linewidths=0.5)
    ax.axvline(1, color="#888", ls="--", lw=0.9); ax.set_xscale("log"); ax.set_xlim(*xlim)
    ax.set_yticks(y); ax.set_yticklabels(labels); ax.tick_params(axis="y", length=0)
    for yi, t in zip(y, txt): ax.text(1.02, yi, t, transform=ax.get_yaxis_transform(), va="center", fontsize=8.3)
    ax.text(1.02, len(labels)-0.45, hdr, transform=ax.get_yaxis_transform(), fontsize=8.3, fontweight="bold")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=9.5, pad=22, x=-0.45); ax.set_xlabel("Reporting odds ratio (log scale)")
    for s in ["top","right"]: ax.spines[s].set_visible(False)
    ax.set_xticks([t for t in [0.5,1,2,5,10,20,50] if xlim[0]<=t<=xlim[1]]); ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v,_: f"{v:g}")); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_ylim(-0.6, len(labels)-0.4)
def _rel(t):
    t = t.replace(", no ICI", ", no ICI recorded").replace(" + ICI", ", ICI co-reported")
    t = t.replace("ICI vs no ICI", "ICI co-reported vs none recorded").replace("— no ICI", "— no ICI recorded").replace("— with ICI", "— ICI co-reported")
    return t
A["Group"] = [_rel(g) for g in A.Group]; wf["Comparison"] = [_rel(g) for g in wf.Comparison]
colA = ["#2c7fb8" if "Sotorasib" in g else "#d95f0e" for g in A.Group]
panel(axes[0], A.Group, A.ROR.values, A.ROR_lo.values, A.ROR_hi.values, [f"{a}/{n}  ·  {r:.2f} ({l:.2f}–{h:.2f})" for a,n,r,l,h in zip(A.n_hepato,A.n_reports,A.ROR,A.ROR_lo,A.ROR_hi)],
      "A  Hepatotoxicity signal versus the full FAERS/AEMS database", (0.7, 60), colA)
panel(axes[1], wf.Comparison, wf.ROR.values, wf.lo.values, wf.hi.values, [f"{a}/{n} vs {b}/{m}  ·  {r:.2f} ({l:.2f}–{h:.2f})" for a,n,b,m,r,l,h in zip(wf.exposed_hep,wf.exposed_n,wf.unexposed_hep,wf.unexposed_n,wf.ROR,wf.lo,wf.hi)],
      "B  Within-cohort comparisons", (0.3, 60), ["#555"]*len(w)+["#999"]*(len(wf)-len(w)), hdr="First vs second group  ·  ROR (95%CI)")
fig.subplots_adjust(left=0.37, right=0.65, top=0.93, bottom=0.08)
for ext in ["tiff","png"]: fig.savefig(os.path.join(M, f"Figure1_forest_plot.{ext}"), dpi=300, **({"pil_kwargs":{"compression":"tiff_lzw"}} if ext=="tiff" else {}))
print(T1.to_string()); print(T2.to_string()); print(T3.to_string()); print(tto, ici_agents)
