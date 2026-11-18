import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10.5})
BLUE,ORANGE,DARK,GREY="#2E75B6","#D95F02","#404040","#7F7F7F"
A=[("Sotorasib, all","657/3515",9.54,8.77,10.39,BLUE),("Sotorasib, no ICI recorded","628/3436",9.29,8.51,10.13,BLUE),("Sotorasib, ICI co-reported","29/79",24.08,15.24,38.05,BLUE),
   ("Adagrasib, all","44/984",1.94,1.44,2.63,ORANGE),("Adagrasib, no ICI recorded","38/963",1.71,1.23,2.36,ORANGE),("Adagrasib, ICI co-reported†","6/21",16.61,6.44,42.80,ORANGE)]
B=[("ICI co-reported vs none recorded — either drug","35/100 vs 666/4399",3.02,1.98,4.59,DARK),("ICI co-reported vs none recorded — sotorasib","29/79 vs 628/3436",2.59,1.63,4.13,DARK),
   ("ICI co-reported vs none recorded — adagrasib†","6/21 vs 38/963",9.74,3.58,26.49,DARK),("Sotorasib vs adagrasib — no ICI recorded","628/3436 vs 38/963",5.44,3.89,7.62,DARK),
   ("Sotorasib vs adagrasib — ICI co-reported†","29/79 vs 6/21",1.45,0.51,4.15,DARK)]
C=[("Docetaxel, all indications: ICI co-reported vs none, ROR","55/559 vs 381/10943",3.03,2.25,4.07,GREY),("Docetaxel, NSCLC: ICI co-reported vs none, ROR","16/169 vs 81/907",1.07,0.61,1.87,GREY),
   ("Sotorasib, NSCLC: adjusted ICI OR","20/55 vs 409/1820",1.50,0.84,2.69,BLUE),("Adagrasib, NSCLC: adjusted ICI OR†","6/19 vs 30/558",6.59,1.87,23.28,ORANGE),("Docetaxel, NSCLC: adjusted ICI OR","16/169 vs 81/907",1.22,0.68,2.21,GREY),
   ("Ratio to docetaxel, sotorasib (primary suspect)","",1.32,0.59,2.97,BLUE),("Ratio to docetaxel, adagrasib (primary suspect)†","",5.02,1.40,17.94,ORANGE),
   ("Ratio to docetaxel, sotorasib (any role)","",1.90,1.16,3.13,BLUE),("Ratio to docetaxel, adagrasib (any role)","",2.68,1.17,6.11,ORANGE)]
panels=[("A  Prespecified: hepatotoxicity signal against all other FAERS/AEMS reports","n/N  ·  ROR (95%CI)",A,(0.7,60),[1,2,5,10,20,50],"Reporting odds ratio (log scale)",False),
        ("B  Prespecified: within-cohort comparisons","First vs second group  ·  ROR (95%CI)",B,(0.3,60),[0.5,1,2,5,10,20,50],"Reporting odds ratio (log scale)",False),
        ("C  Exploratory: docetaxel reference and adjusted NSCLC models","Events  ·  estimate (95%CI)",C,(0.3,60),[0.5,1,2,5,10,20,50],"Odds ratio or ratio of odds ratios (log scale)",True)]
h=[len(p[2]) for p in panels]
fig,axs=plt.subplots(3,1,figsize=(10.5,14.2),gridspec_kw=dict(height_ratios=[x+1.6 for x in h],hspace=0.62))
for ax,(title,head,rows,lim,ticks,xl,expl) in zip(axs,panels):
    n=len(rows); ys=np.arange(n)[::-1]
    for y,(lab,ev,e,lo,hi,col) in zip(ys,rows):
        ax.plot([lo,hi],[y,y],color="#404040",lw=1.3); ax.plot([lo,lo],[y-.12,y+.12],color="#404040",lw=1.3); ax.plot([hi,hi],[y-.12,y+.12],color="#404040",lw=1.3)
        if expl: ax.plot(e,y,marker="D",ms=8,mfc="white",mec=col,mew=1.8,zorder=3)
        else: ax.plot(e,y,marker="s",ms=8.5,mfc=col,mec="black",mew=0.6,zorder=3)
        ax.text(1.02,y,(ev+"  ·  " if ev else "")+f"{e:.2f} ({lo:.2f}–{hi:.2f})",transform=ax.get_yaxis_transform(),va="center",fontsize=10)
    ax.text(1.02,n-0.35,head,transform=ax.get_yaxis_transform(),fontweight="bold",fontsize=10)
    ax.axvline(1,ls="--",color="grey",lw=1); ax.set_xscale("log"); ax.set_xlim(*lim); ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks]); ax.minorticks_off()
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows]); ax.set_ylim(-0.7,n-0.3+0.4); ax.set_xlabel(xl)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    ax.annotate(title,xy=(0.01,1.0),xycoords=("figure fraction","axes fraction"),xytext=(0,14),textcoords="offset points",fontweight="bold",fontsize=12,ha="left",va="bottom")
fig.text(-0.17,0.045,"Squares: prespecified crude estimates. Open diamonds: exploratory estimates added after the initial analysis.\n† Based on 6 ICI co-reported adagrasib hepatotoxicity reports.",fontsize=9.5)
fig.savefig("f1.png",dpi=300,bbox_inches="tight",pad_inches=0.15); print("ok")

from PIL import Image
Image.open("f1.png").convert("RGB").save("Figure2_forest_plot.png",optimize=True)
