"""Step 9: full sensitivity grid (Table 6) on the final cohort."""
import pandas as pd, numpy as np, os, math
B=os.path.expanduser("~/mnt/Downloads/FAERS_KRAS_G12C"); R=os.path.join(B,"results")
c=pd.read_csv(os.path.join(R,"cohort_reports_v2_flags.csv"),dtype={"primaryid":str},keep_default_na=False)
for k in ["ici_exposed","hepato","ici_suspect","study_report","switch"]: c[k]=c[k].astype(str)=="True"
reac=pd.read_csv(os.path.join(os.environ["FAERS_X"],"reac.csv"),dtype=str,keep_default_na=False)
NARROW={s.lower() for s in ["Hepatotoxicity","Drug-induced liver injury","Liver injury","Hepatitis","Hepatitis acute","Hepatitis fulminant","Hepatitis cholestatic","Immune-mediated hepatitis","Immune-mediated hepatic disorder","Autoimmune hepatitis","Hepatic failure","Acute hepatic failure","Acute on chronic liver failure","Hepatocellular injury","Cholestatic liver injury","Mixed liver injury","Hepatic cytolysis","Hepatic necrosis","Jaundice","Immune-mediated cholangitis"]}
c["hep_narrow"]=c.primaryid.isin(reac[reac.pt.str.strip().str.lower().isin(NARROW)].primaryid)
def ror(a,n1,cc,n0):
    b,d=n1-a,n0-cc; a,b,cc,d=[x+0.5 if 0 in (a,b,cc,d) else x for x in (a,b,cc,d)]; r=(a/b)/(cc/d); se=math.sqrt(1/a+1/b+1/cc+1/d)
    return f"{r:.2f} ({r*math.exp(-1.96*se):.2f}-{r*math.exp(1.96*se):.2f})"
def grid(s,o="hepato"):
    out=[]
    for d in ["sotorasib","adagrasib"]:
        t=s[s.target_drug==d]; m=t.ici_exposed; out.append(ror(int(t[m][o].sum()),int(m.sum()),int(t[~m][o].sum()),int((~m).sum())))
    for g in [False,True]:
        t=s[s.ici_exposed==g]; m=t.target_drug=="sotorasib"; out.append(ror(int(t[m][o].sum()),int(m.sum()),int(t[~m][o].sum()),int((~m).sum())))
    return out
HP={"MD","PH","HP","OT","RN"}
sets=[("Main analysis",c,"hepato"),("NSCLC indication only",c[c.tumor=="NSCLC"],"hepato"),("Health-professional reporters only",c[c.occp.isin(HP)],"hepato"),
 ("Clinical liver injury terms only",c,"hep_narrow"),("Reports received from 2023",c[c.quarter>="2023Q1"],"hepato"),
 ("ICI recorded as concomitant only",c[~c.ici_suspect],"hepato"),("ICI recorded as suspect or interacting",c[c.ici_suspect|~c.ici_exposed],"hepato"),
 ("Excluding likely trial reports",c[~c.study_report],"hepato"),("Excluding possible switch cases",c[~c.switch],"hepato")]
rows=[[name,f"{len(s)}"]+grid(s,o) for name,s,o in sets]
T=pd.DataFrame(rows,columns=["Analysis","Reports","ICI vs no ICI, sotorasib","ICI vs no ICI, adagrasib","Sotorasib vs adagrasib, no ICI","Sotorasib vs adagrasib, with ICI"])
T.to_csv(os.path.join(B,"manuscript","Table6_sensitivity_grid.csv"),index=False); pd.set_option("display.width",250); print(T.to_string(index=False))
