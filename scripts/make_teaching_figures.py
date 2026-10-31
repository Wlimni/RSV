"""
Script Name: make_teaching_figures.py
Purpose:     Draw the explanatory diagrams used in the textbook guide (docs/). These are
             teaching illustrations, not results: genome map, virion schematic, the O/E
             idea, codon positions, and why a pooled trend can mislead.
Inputs:      results/gene_summary.tsv (for real gene lengths)
Outputs:     docs/figures/t1..t5 *.png
Notes:       Run: .venv/bin/python scripts/make_teaching_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Rectangle, FancyArrowPatch

OUT = "docs/figures"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False})
BLUE, RED, GREY = "#2C7FB8", "#C62824", "#6B7280"


# ---------------------------------------------------------------- T1 genome map
def genome_map():
    genes = [("NS1", 420, "blocks interferon"), ("NS2", 375, "blocks interferon"),
             ("N", 1176, "coats the RNA"), ("P", 726, "polymerase partner"),
             ("M", 771, "matrix, structure"), ("SH", 195, "small, ion channel"),
             ("G", 969, "attaches to cell"), ("F", 1725, "fuses membranes"),
             ("M2", 585, "transcription"), ("L", 6498, "copies the genome")]
    total = sum(g[1] for g in genes)
    fig, ax = plt.subplots(figsize=(11, 2.9))
    x = 0.0
    for name, ln, role in genes:
        w = ln / total
        colour = RED if name in ("NS1", "NS2") else (BLUE if name in ("G", "F") else "#9CA3AF")
        ax.add_patch(FancyBboxPatch((x, .45), w - .003, .26,
                                    boxstyle="round,pad=0,rounding_size=.01",
                                    facecolor=colour, edgecolor="black", linewidth=.7))
        ax.text(x + w / 2, .58, name, ha="center", va="center", fontsize=9,
                color="white", weight="bold")
        ax.text(x + w / 2, .38, f"{ln:,}", ha="center", va="top", fontsize=7, color=GREY)
        if name in ("NS1", "NS2", "G", "F", "L"):
            ax.text(x + w / 2, .78, role, ha="center", va="bottom", fontsize=7,
                    color="black", rotation=20)
        x += w
    ax.annotate("", xy=(1, .30), xytext=(0, .30),
                arrowprops=dict(arrowstyle="->", color="black", lw=1))
    ax.text(0, .22, "3' end", fontsize=8, color=GREY)
    ax.text(1, .22, "5' end", fontsize=8, color=GREY, ha="right")
    ax.text(.5, .09, "about 15,200 nucleotides in total  ·  numbers are gene length in nucleotides",
            ha="center", fontsize=8, color=GREY)
    ax.set_xlim(-.01, 1.01); ax.set_ylim(0, 1.05); ax.axis("off")
    ax.set_title("The RSV genome: ten genes, always in this order",
                 fontsize=11, weight="bold", loc="left")
    fig.tight_layout(); fig.savefig(f"{OUT}/t1_genome_map.png", dpi=200); plt.close(fig)


# ---------------------------------------------------------------- T2 virion
def virion():
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.add_patch(Ellipse((0, 0), 5.6, 4.2, facecolor="#E8EEF4",
                         edgecolor="#1F3A52", linewidth=1.6))
    # spikes
    for ang in np.linspace(0, 2 * np.pi, 26, endpoint=False):
        x, y = 2.8 * np.cos(ang), 2.1 * np.sin(ang)
        xo, yo = 3.25 * np.cos(ang), 2.45 * np.sin(ang)
        c = BLUE if int(ang * 4) % 2 else RED
        ax.plot([x, xo], [y, yo], color=c, lw=2.1, solid_capstyle="round")
    # coiled nucleocapsid
    t = np.linspace(0, 6.4 * np.pi, 600)
    ax.plot(1.45 * np.cos(t) * (1 - t / 70), .95 * np.sin(t) * (1 - t / 70),
            color="#374151", lw=1.5)
    ax.text(0, -.1, "RNA genome\ncoated in N protein", ha="center", va="center",
            fontsize=8, color="#111827",
            bbox=dict(boxstyle="round,pad=.3", fc="white", ec="none", alpha=.85))
    lbl = [(3.6, 2.0, "F protein\nfuses the virus to the cell\n(target of nirsevimab)", BLUE),
           (-3.9, 1.6, "G protein\ngrabs the cell surface\nmost variable gene", RED),
           (3.3, -2.0, "envelope\nborrowed from the\nhost cell membrane", "#1F3A52")]
    for x, y, s, c in lbl:
        ax.text(x, y, s, fontsize=7.5, color=c, ha="center", va="center")
    ax.set_xlim(-5.2, 5.2); ax.set_ylim(-3.1, 3.1); ax.axis("off")
    ax.set_title("What an RSV particle looks like", fontsize=11, weight="bold")
    fig.tight_layout(); fig.savefig(f"{OUT}/t2_virion.png", dpi=200); plt.close(fig)


# ------------------------------------------------------------- T3 O/E explained
def oe_explained():
    fig, ax = plt.subplots(2, 1, figsize=(9, 4.6),
                           gridspec_kw={"height_ratios": [1, 1.25]})
    seq = "ATGCGTTACGGCATCCGTAAGCTACGA"
    a0 = ax[0]
    for i, ch in enumerate(seq):
        a0.text(i, .55, ch, ha="center", va="center", fontsize=11, family="monospace")
    for i in range(len(seq) - 1):
        if seq[i] == "C" and seq[i + 1] == "G":
            a0.add_patch(Rectangle((i - .45, .32), 1.9, .46, facecolor="#FDE68A",
                                   edgecolor="#B45309", linewidth=1, zorder=0))
    a0.text(len(seq) / 2, .05, "every C immediately followed by a G is highlighted — these are the CpG pairs",
            ha="center", fontsize=8, color=GREY)
    a0.set_xlim(-1, len(seq)); a0.set_ylim(0, 1); a0.axis("off")
    a0.set_title("Step 1 — count the pairs", fontsize=10, weight="bold", loc="left")

    a1 = ax[1]
    a1.axis("off")
    a1.set_title("Step 2 — compare with what chance alone would give",
                 fontsize=10, weight="bold", loc="left")
    txt = [
        (r"if 15% of letters are C and 15% are G, then by chance alone", 0.80, 8.5, "black"),
        (r"expected CpG   =   0.15  $\times$  0.15   =   0.0225   (2.25% of positions)", 0.62, 10, "#1F3A52"),
        (r"suppose we actually count   0.005   (0.5% of positions)", 0.44, 9, "black"),
        (r"O/E   =   0.005  $\div$  0.0225   =   0.22", 0.24, 12, "#B45309"),
        ("a value well below 1 means the virus is avoiding CpG", 0.06, 8.5, GREY)]
    for s, y, fs, c in txt:
        a1.text(.5, y, s, ha="center", fontsize=fs, color=c,
                weight="bold" if fs >= 10 else "normal")
    fig.tight_layout(); fig.savefig(f"{OUT}/t3_oe_explained.png", dpi=200); plt.close(fig)


# ------------------------------------------------------- T4 codon positions
def codon_positions():
    fig, ax = plt.subplots(figsize=(9, 3.2))
    codons = [("ATG", "Met"), ("GCG", "Ala"), ("TTA", "Leu"), ("CGT", "Arg"), ("AAC", "Asn")]
    x = 0
    for cod, aa in codons:
        for j, ch in enumerate(cod):
            c = "#1F5C8B" if j == 2 else "#80B1D3"
            ax.add_patch(Rectangle((x, .45), .9, .35, facecolor=c,
                                   edgecolor="white", linewidth=1.4))
            ax.text(x + .45, .625, ch, ha="center", va="center", color="white",
                    fontsize=11, weight="bold", family="monospace")
            x += 1
        ax.text(x - 1.5, .34, aa, ha="center", va="top", fontsize=8, color=GREY)
        x += .55
    ax.add_patch(Rectangle((0, .92), .9, .14, facecolor="#80B1D3", edgecolor="none"))
    ax.text(1.1, .99, "positions 1 and 2  —  changing these usually changes the protein  (GC12)",
            va="center", fontsize=8.5)
    ax.add_patch(Rectangle((0, .14), .9, .14, facecolor="#1F5C8B", edgecolor="none"))
    ax.text(1.1, .21, "position 3  —  changing this usually does NOT change the protein  (GC3)",
            va="center", fontsize=8.5)
    ax.set_xlim(-.3, x); ax.set_ylim(.05, 1.15); ax.axis("off")
    ax.set_title("Every codon is three letters — and the third one is special",
                 fontsize=11, weight="bold", loc="left")
    fig.tight_layout(); fig.savefig(f"{OUT}/t4_codon_positions.png", dpi=200); plt.close(fig)


# ------------------------------------------------------- T5 pooling can mislead
def pooling():
    rng = np.random.default_rng(3)
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 3.6), sharey=True)
    years = np.linspace(2016, 2025, 260)
    lo = 1 - (years - 2016) / 9
    grpA = rng.random(260) < lo
    yA = 0.342 - 0.0004 * (years - 2016) + rng.normal(0, .0006, 260)
    yB = 0.347 - 0.0004 * (years - 2016) + rng.normal(0, .0006, 260)
    y = np.where(grpA, yA, yB)

    ax[0].scatter(years, y, s=7, c="#9CA3AF", alpha=.75)
    m, b = np.polyfit(years, y, 1)
    ax[0].plot(years, m * years + b, color="#B45310", lw=2)
    ax[0].set_title(f"Pooled: looks like it is RISING\n(slope {m:+.5f}/yr)",
                    fontsize=9.5, weight="bold", color="#B45310")
    ax[0].set_ylabel("composition measure")

    for mask, c, nm in [(grpA, BLUE, "lineage 1"), (~grpA, RED, "lineage 2")]:
        ax[1].scatter(years[mask], y[mask], s=7, c=c, alpha=.75, label=nm)
        mm, bb = np.polyfit(years[mask], y[mask], 1)
        ax[1].plot(years[mask], mm * years[mask] + bb, color=c, lw=2)
    ax[1].set_title("Split by lineage: BOTH are falling", fontsize=9.5,
                    weight="bold", color="#1F3A52")
    ax[1].legend(frameon=False, fontsize=8, loc="lower right")
    for a in ax:
        a.set_xlabel("year")
    fig.suptitle("Why a trend must be measured inside each lineage", fontsize=11,
                 weight="bold", x=.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, .93])
    fig.savefig(f"{OUT}/t5_pooling.png", dpi=200); plt.close(fig)


for fn in (genome_map, virion, oe_explained, codon_positions, pooling):
    fn(); print("drew", fn.__name__)
print("teaching figures written to", OUT)
