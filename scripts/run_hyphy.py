"""
Script Name: run_hyphy.py
Purpose:     Run HyPhy SLAC on each gene within each RSV subtype, to ask whether A and B differ
             in the strength of selection on the protein (dN/dS) in the genes whose CpG
             content differs (NS1, NS2) versus the genes where it does not.
Inputs:      data/processed/cds/{A,B}_{gene}.fasta (from build_codon_data.py)
Outputs:     results/hyphy_summary.tsv  one row per subtype x gene
             data/processed/hyphy/       alignments, trees and SLAC JSON (not pushed)
Notes:       Needs mafft (PATH), IQ-TREE 3 and HyPhy 2.5 from the conda env named `hyphy`
             (conda create -n hyphy -c conda-forge -c bioconda hyphy iqtree).
             Alignment is done on the protein and back-translated, so the codon alignment is
             always in frame. Each gene uses a random sample of --n genomes (seeded) so the tree
             stays small; SLAC is a fast counting method, adequate for a gene-wide dN/dS.
             Run: .venv/bin/python scripts/run_hyphy.py [--n 100]  (AB rows use n/2 per subtype)
"""
import os, sys, json, random, subprocess, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np, pandas as pd
from Bio.Seq import Seq
from rsv import genes

ap = argparse.ArgumentParser()
ap.add_argument("--n", type=int, default=100)
ap.add_argument("--env", default="/opt/anaconda3/envs/hyphy/bin")
args = ap.parse_args()
GENES = ["NS1", "NS2", "N", "P", "M", "SH", "F", "M2-1"]
W = "data/processed/hyphy"
os.makedirs(W, exist_ok=True)


def codon_align(seqs, path):
    """Align at protein level with mafft, then put each sequence's codons back under it."""
    prot = {k: str(Seq(v).translate()) for k, v in seqs.items()}
    with open(path + ".prot.fa", "w") as fh:
        fh.write("".join(f">{k}\n{p}\n" for k, p in prot.items()))
    out = subprocess.run(["mafft", "--quiet", "--auto", path + ".prot.fa"],
                         capture_output=True, text=True, check=True).stdout
    aln, name, buf = {}, None, []
    for line in out.split("\n"):
        if line.startswith(">"):
            if name:
                aln[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    aln[name] = "".join(buf)
    with open(path, "w") as fh:
        for k, p in aln.items():
            cod, i = [], 0
            for ch in p:
                if ch == "-":
                    cod.append("---")
                else:
                    cod.append(seqs[k][3 * i:3 * i + 3]); i += 1
            fh.write(f">{k}\n{''.join(cod)}\n")


rows = []
# "AB" puts both subtypes on one tree: its dN/dS is dominated by the deep A-versus-B split
for sub in ["A", "B", "AB"]:
    for gene in GENES:
        d = {}
        for part in sub:
            x = genes.read_fasta(f"data/processed/cds/{part}_{gene}.fasta")
            x = {k: v[:-3] for k, v in x.items() if v[-3:] in ("TAA", "TAG", "TGA")}   # strip stop
            x = {k: v for k, v in x.items() if set(v) <= set("ACGT")}                  # SLAC rejects N
            random.seed(0)
            take = args.n // len(sub)
            d.update({k: x[k] for k in random.sample(sorted(x), min(take, len(x)))})
        keep = sorted(d)
        aln = f"{W}/{sub}_{gene}.aln.fa"
        codon_align({k: d[k] for k in keep}, aln)
        subprocess.run([f"{args.env}/iqtree3", "-s", aln, "-m", "GTR+G", "-fast", "-T", "2",
                        "--prefix", f"{W}/{sub}_{gene}", "-redo", "-quiet"], check=True)
        js = f"{W}/{sub}_{gene}.slac.json"
        subprocess.run([f"{args.env}/hyphy", "slac", "--alignment", aln,
                        "--tree", f"{W}/{sub}_{gene}.treefile", "--branches", "All",
                        "--samples", "0", "--output", js], capture_output=True, text=True, check=True)
        j = json.load(open(js))
        hdr = [h[0] for h in j["MLE"]["headers"]]
        site = pd.DataFrame(j["MLE"]["content"]["0"]["by-branch"]["AVERAGED"], columns=hdr)
        S, NS = site["S"].sum(), site["N"].sum()
        eS, eNS = site["ES"].sum(), site["EN"].sum()
        dS, dN = S / eS, NS / eNS
        # HyPhy's own gene-wide dN/dS (a single rate ratio fitted by ML on the whole tree)
        om = j["fits"]["Global MG94xREV"]["Rate Distributions"]
        omega = next(v for k, v in om.items() if "non-synonymous/synonymous" in k)[0][0]
        rows.append({"Subtype": sub, "gene": gene, "n_seq": len(keep), "codons": len(site),
                     "observed_S": S, "observed_NS": NS, "dS": dS, "dN": dN, "omega": omega,
                     "sites_purifying": int((site["P [dN/dS < 1]"] < 0.1).sum()),
                     "sites_diversifying": int((site["P [dN/dS > 1]"] < 0.1).sum())})
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv("results/hyphy_summary.tsv", sep="\t", index=False)
