"""
Script Name: make_notebooks.py
Purpose:     Generate the four analysis notebooks (QC, composition, time trends, turnover
             + model) from code defined here, so notebook content is version-controlled as
             plain Python and regenerated reproducibly.
Inputs:      none (reads nothing; writes notebook files)
Outputs:     notebooks/01_qc.ipynb ... 04_turnover_model.ipynb
Notes:       Notebooks import from src/rsv and read data/processed/rsv_features.tsv.
             Run: python3 scripts/make_notebooks.py   then execute with nbconvert.
"""
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
import os

os.makedirs("notebooks", exist_ok=True)
SETUP = """\
import sys, os
sys.path.insert(0, os.path.abspath("../src"))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns
sns.set_theme(style="whitegrid")
pd.set_option("display.width", 120)
from rsv import io, composition
"""


def build(path, cells):
    nb = new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nbf.write(nb, path)
    print("wrote", path)


# ---- 01 QC ----
build("notebooks/01_qc.ipynb", [
    new_markdown_cell("# 01 — Data QC\n\nMerged RSV-A/B GISAID set: counts, sequence length, "
                      "ambiguous bases, collection date, and geography. Justifies the QC filter."),
    new_code_cell(SETUP),
    new_code_cell('seqs = io.read_fasta("../data/full/RSV_full.fasta")\n'
                  'meta = io.load_metadata("../data/full/RSV_full_metadata.tsv")\n'
                  'print("sequences:", len(seqs), "| metadata rows:", len(meta))\n'
                  'meta[["Subtype","year","region"]].describe(include="all")'),
    new_markdown_cell("### Sequence length by subtype"),
    new_code_cell('meta["len"] = meta["Sequence Length"]\n'
                  'fig, ax = plt.subplots(figsize=(7,4))\n'
                  'for s,c in zip("AB",["#2c7fb8","#de2d26"]):\n'
                  '    meta.loc[meta.Subtype==s,"len"].plot.hist(bins=120, alpha=.5, ax=ax, label=s, color=c)\n'
                  'ax.set_xlim(14800,15600); ax.axvline(14900,ls="--",c="k"); ax.axvline(15400,ls="--",c="k")\n'
                  'ax.set_xlabel("sequence length (bp)"); ax.legend(title="subtype"); plt.show()'),
    new_markdown_cell("### Apply the QC filter"),
    new_code_cell('clean, report = io.qc_filter(meta, seqs)\n'
                  'pd.Series(report)'),
    new_markdown_cell("### Sequences per year and region (QC-passing)"),
    new_code_cell('display(clean.groupby(["year"]).size().rename("n").to_frame().T)\n'
                  'clean.groupby(["region","Subtype"]).size().unstack(fill_value=0)'),
])

# ---- 02 composition ----
build("notebooks/02_composition.ipynb", [
    new_markdown_cell("# 02 — Composition by subtype and lineage\n\nGC content and the four GC "
                      "dinucleotide O/E ratios, compared between RSV-A and RSV-B and across lineages."),
    new_code_cell(SETUP),
    new_code_cell('d = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
                  'FEATS = ["gc_content","oe_CG","oe_GC","oe_GG","oe_CC"]\n'
                  'd.groupby("Subtype")[FEATS].agg(["mean","std"]).round(4)'),
    new_markdown_cell("### A vs B: is the difference significant? (Mann–Whitney U)"),
    new_code_cell('from scipy.stats import mannwhitneyu\n'
                  'for f in FEATS:\n'
                  '    a=d[d.Subtype=="A"][f]; b=d[d.Subtype=="B"][f]\n'
                  '    u,p = mannwhitneyu(a,b)\n'
                  '    print(f"{f:10s} A={a.mean():.4f} B={b.mean():.4f}  p={p:.2e}")'),
    new_markdown_cell("### Composition across the most common lineages"),
    new_code_cell('top = (d.groupby(["Subtype","Lineage"]).size().rename("n").reset_index()\n'
                  '        .sort_values("n",ascending=False).groupby("Subtype").head(8))\n'
                  'lin = d.merge(top[["Subtype","Lineage"]], on=["Subtype","Lineage"])\n'
                  'summ = lin.groupby(["Subtype","Lineage"])[FEATS+["gc_content"]].mean()\n'
                  'summ = lin.groupby(["Subtype","Lineage"]).agg(n=("gc_content","size"),\n'
                  '        gc=("gc_content","mean"), oe_CG=("oe_CG","mean")).round(4)\n'
                  'summ.to_csv("../results/composition_by_lineage_summary.tsv", sep="\\t")\n'
                  'summ'),
    new_code_cell('g = sns.catplot(data=lin, x="oe_CG", y="Lineage", col="Subtype", kind="box",\n'
                  '    sharey=False, height=4, aspect=1, color="#6baed6")\n'
                  'g.set_titles("RSV-{col_name}"); plt.show()'),
])

# ---- 03 time trends ----
build("notebooks/03_time_trends.ipynb", [
    new_markdown_cell("# 03 — Composition over time\n\nYearly trends in GC content and CpG O/E "
                      "(2015–2026), including the COVID-19 disruption (2020–2021)."),
    new_code_cell(SETUP),
    new_code_cell('d = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
                  'd = d[(d.year>=2015)&(d.year<=2026)]\n'
                  'yr = d.groupby(["Subtype","year"]).agg(n=("gc_content","size"),\n'
                  '      gc=("gc_content","mean"), oe_CG=("oe_CG","mean")).reset_index()\n'
                  'yr = yr[yr.n>=20]\n'
                  'yr.head()'),
    new_markdown_cell("### Trend correlation (mean composition vs year)"),
    new_code_cell('from scipy.stats import spearmanr\n'
                  'for s in "AB":\n'
                  '    x=yr[yr.Subtype==s]\n'
                  '    for f in ["gc","oe_CG"]:\n'
                  '        r,p = spearmanr(x.year, x[f])\n'
                  '        print(f"{s} {f:6s} rho={r:+.3f} p={p:.3f}")'),
    new_code_cell('fig,ax=plt.subplots(1,2,figsize=(12,4))\n'
                  'for s,c in zip("AB",["#2c7fb8","#de2d26"]):\n'
                  '    x=yr[yr.Subtype==s]\n'
                  '    ax[0].plot(x.year,x.gc,"-o",color=c,label=s)\n'
                  '    ax[1].plot(x.year,x.oe_CG,"-o",color=c,label=s)\n'
                  'ax[0].set_title("mean GC content"); ax[1].set_title("mean CpG O/E")\n'
                  'ax[0].axvspan(2020,2021,color="grey",alpha=.15)\n'
                  'ax[1].axvspan(2020,2021,color="grey",alpha=.15)\n'
                  'ax[0].legend(title="subtype"); plt.show()'),
    new_markdown_cell("Grey band = COVID-19 disruption. Note A and B diverge in CpG O/E after 2021."),
])

# ---- 04 turnover + model ----
build("notebooks/04_turnover_model.ipynb", [
    new_markdown_cell("# 04 — Lineage turnover and prediction\n\nDo expanding lineages differ in "
                      "composition from fading ones? Can composition predict subtype and year?"),
    new_code_cell(SETUP),
    new_code_cell('d = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
                  'FEATS = ["gc_content","freq_CG","freq_GC","freq_GG","freq_CC",\n'
                  '         "oe_CG","oe_GC","oe_GG","oe_CC"]'),
    new_markdown_cell("### Expanding vs fading lineages\n"
                      "Compare each lineage's share in an early window (2017–2020) vs a late "
                      "window (2023–2026); label those gaining/losing >3 percentage points."),
    new_code_cell('def turnover(sub):\n'
                  '    x=d[d.Subtype==sub]\n'
                  '    early=x[x.year.between(2017,2020)].Lineage.value_counts(normalize=True)\n'
                  '    late =x[x.year.between(2023,2026)].Lineage.value_counts(normalize=True)\n'
                  '    chg=pd.DataFrame({"early":early,"late":late}).fillna(0)\n'
                  '    chg["delta"]=chg.late-chg.early\n'
                  '    comp=x.groupby("Lineage")[FEATS].mean()\n'
                  '    chg=chg.join(comp)\n'
                  '    chg["group"]=np.where(chg.delta>0.03,"expanding",np.where(chg.delta<-0.03,"fading","stable"))\n'
                  '    return chg\n'
                  'for s in "AB":\n'
                  '    t=turnover(s)\n'
                  '    g=t[t.group!="stable"].groupby("group")[["oe_CG","gc_content"]].mean().round(4)\n'
                  '    print(f"=== RSV-{s} ==="); print(g, "\\n")'),
    new_markdown_cell("### Predict subtype and year from composition (Random Forest)"),
    new_code_cell('from sklearn.ensemble import RandomForestClassifier\n'
                  'from sklearn.model_selection import cross_val_score\n'
                  'sm = d.sample(min(8000,len(d)), random_state=1)\n'
                  'for tgt in ["Subtype","year"]:\n'
                  '    y=sm[tgt].astype(str)\n'
                  '    clf=RandomForestClassifier(n_estimators=200,n_jobs=-1,random_state=0)\n'
                  '    acc=cross_val_score(clf,sm[FEATS],y,cv=4).mean()\n'
                  '    base=y.value_counts(normalize=True).max()\n'
                  '    print(f"predict {tgt:8s}: acc={acc:.3f}  (majority baseline {base:.3f})")'),
    new_code_cell('clf=RandomForestClassifier(n_estimators=300,n_jobs=-1,random_state=0)\n'
                  'clf.fit(sm[FEATS], sm["Subtype"])\n'
                  'imp=pd.Series(clf.feature_importances_,index=FEATS).sort_values()\n'
                  'imp.plot.barh(color="#6baed6"); plt.title("feature importance (subtype)"); plt.show()'),
])
print("done")

# ---- 05 within vs between lineage, growth, founder, robustness ----
build("notebooks/05_trend_decomposition.ipynb", [
    new_markdown_cell("# 05 — Is the trend evolution or lineage replacement?\n\n"
                      "A rising overall trend can come from genomes changing **inside** lineages, "
                      "or from high-value lineages **replacing** low-value ones. Adding lineage "
                      "fixed effects absorbs the replacement part, leaving the within-lineage slope.\n\n"
                      "Also covers lineage growth rate, founder composition, and robustness checks."),
    new_code_cell(SETUP + "from rsv import trends"),
    new_code_cell('d = pd.read_csv("../data/processed/rsv_features.tsv", sep="\\t")\n'
                  'd = d[(d.year>=2015)&(d.year<=2026)]\n'
                  'F = ["gc_content","oe_CG","oe_GC","oe_GG","oe_CC"]\n'
                  'print(len(d), "sequences")'),
    new_markdown_cell("## Step 1 — decomposition\n"
                      "Opposite signs between `slope_total` and `slope_within` mean the raw "
                      "trend is an artefact of turnover (a Simpson's paradox)."),
    new_code_cell('dec = pd.concat([trends.decompose_trend(s,f).assign(Subtype=sub)\n'
                  '                 for sub,s in d.groupby("Subtype") for f in F], ignore_index=True)\n'
                  'dec[["Subtype","feature","slope_total","slope_within","pct_from_turnover"]].round(6)'),
    new_markdown_cell("## Per-lineage slopes\nDo individual lineages move in the same direction?"),
    new_code_cell('sl = trends.per_lineage_slopes(d, "oe_CG")\n'
                  'print(f"{len(sl)} lineages, {sl.sig.sum()} significant")\n'
                  'sl.groupby("Subtype").agg(n=("slope","size"), median_slope=("slope","median"),\n'
                  '                          n_positive=("slope", lambda x:(x>0).sum()))'),
    new_markdown_cell("## Step 3 — lineage growth rate vs composition\n"
                      "Growth rate = logistic change in a lineage's share of its subtype per year."),
    new_code_cell('growth = trends.lineage_growth_rate(d)\n'
                  'tab, corr = trends.growth_vs_composition(d, growth, F)\n'
                  'from statsmodels.stats.multitest import multipletests\n'
                  'corr["p_adj"] = multipletests(corr["p"], method="fdr_bh")[1]\n'
                  'corr.sort_values("p").round(4)'),
    new_markdown_cell("## Step 4 — founder composition vs later success\n"
                      "Tree-free analogue of the paper's 'earliest viruses' analysis: does a "
                      "lineage's founding composition predict how successful it later becomes?"),
    new_code_cell('fnd = trends.founder_composition(d, F)\n'
                  '_, fres = trends.founder_vs_success(fnd, growth, F)\n'
                  'fres.head(10).round(4)'),
    new_markdown_cell("## Step 5 — robustness"),
    new_code_cell('rob = trends.robustness(d, "oe_CG")\n'
                  'rob[["Subtype","scenario","n","slope_total","p_total","slope_within","p_within"]].round(6)'),
])

# ---- 06 per-gene ----
build("notebooks/06_per_gene.ipynb", [
    new_markdown_cell("# 06 — Composition gene by gene\n\n"
                      "Each genome is aligned to its subtype's annotated reference and cut into "
                      "the 10 RSV genes, so we can see which gene carries the signal. Alignment "
                      "keeps insertions, so the ON1/BA G-gene duplications are retained."),
    new_code_cell(SETUP + "from rsv import trends"),
    new_code_cell('g = pd.read_csv("../data/processed/rsv_gene_features.tsv", sep="\\t")\n'
                  'g = g[(g.year>=2015)&(g.year<=2026)]\n'
                  'print(len(g), "gene-rows;", g["Accession ID"].nunique(), "sequences")\n'
                  'GENES = ["NS1","NS2","N","P","M","SH","G","F","M2","L"]\n'
                  'g.groupby(["Subtype","gene"]).agg(n=("gc_content","size"),\n'
                  '   gc=("gc_content","mean"), oe_CG=("oe_CG","mean")).round(4)'),
    new_markdown_cell("### Which genes carry the within-lineage CpG O/E trend?"),
    new_code_cell('rows=[]\n'
                  'for (sub,gene),s in g.groupby(["Subtype","gene"]):\n'
                  '    if len(s)<500 or s.Lineage.nunique()<3: continue\n'
                  '    r=trends.decompose_trend(s,"oe_CG"); r.insert(0,"gene",gene); r.insert(0,"Subtype",sub)\n'
                  '    rows.append(r)\n'
                  'gd=pd.concat(rows,ignore_index=True)\n'
                  'from statsmodels.stats.multitest import multipletests\n'
                  'gd["p_within_adj"]=multipletests(gd["p_within"],method="fdr_bh")[1]\n'
                  'gd["gene"]=pd.Categorical(gd["gene"],GENES,ordered=True)\n'
                  'gd.sort_values(["Subtype","gene"])[["Subtype","gene","n","slope_total","slope_within","p_within_adj"]].round(6)'),
])
