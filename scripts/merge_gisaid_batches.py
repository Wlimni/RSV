"""
Script Name: merge_gisaid_batches.py
Purpose:     Merge per-date-range GISAID EpiRSV downloads (RSV A and B) into one FASTA and one metadata table.
Inputs:      Batch folders (default: every folder under Raw/), each holding one FASTA and up to three TSVs
             (patient status, sequencing technology, dates and location) with GISAID default names.
Outputs:     <out_dir>/RSV_full.fasta, <out_dir>/RSV_full_metadata.tsv, <out_dir>/merge_report.txt
Notes:       TSV type is detected from its columns (filenames are not reliable). Join key = Accession ID.
             Run: python3 -I scripts/merge_gisaid_batches.py [--raw Raw] [--out data/full]
"""
import argparse, csv, glob, os

ap = argparse.ArgumentParser()
ap.add_argument("--raw", default="data/Raw")
ap.add_argument("--out", default="data/full")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

meta, fasta = {}, {}
n_dup_fasta = n_dup_meta = 0
n_rows = 0

for b in sorted(glob.glob(os.path.join(args.raw, "*", ""))):
    for f in sorted(glob.glob(os.path.join(b, "*.tsv"))):
        with open(f, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                n_rows += 1
                acc = row["Accession ID"]
                merged = meta.setdefault(acc, {})
                for k, v in row.items():
                    if k and v not in (None, "") and k not in merged:  # first non-empty value wins
                        merged[k] = v
    for f in sorted(glob.glob(os.path.join(b, "*.fasta"))):
        h, s = None, []
        def flush():
            global n_dup_fasta
            if h:
                acc = h.split("|")[1]
                n_dup_fasta += acc in fasta
                fasta[acc] = (h, "".join(s).upper())
        with open(f) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    flush()
                    h, s = line[1:], []
                else:
                    s.append(line)
        flush()

cols = []
for r in meta.values():
    for k in r:
        if k not in cols:
            cols.append(k)

with open(os.path.join(args.out, "RSV_full_metadata.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(meta.values())
with open(os.path.join(args.out, "RSV_full.fasta"), "w") as fh:
    for h, s in fasta.values():
        fh.write(f">{h}\n{s}\n")

sub = {}
for r in meta.values():
    sub[r.get("Subtype", "NA")] = sub.get(r.get("Subtype", "NA"), 0) + 1
report = [
    f"fasta records: {len(fasta)} | duplicate accessions across fasta files: {n_dup_fasta}",
    f"metadata accessions: {len(meta)} (from {n_rows} TSV rows)",
    f"in fasta, no metadata: {len(set(fasta) - set(meta))}",
    f"in metadata, no fasta: {len(set(meta) - set(fasta))}",
    f"subtype counts: {sub}",
    f"columns: {cols}",
]
open(os.path.join(args.out, "merge_report.txt"), "w").write("\n".join(report) + "\n")
print("\n".join(report))
