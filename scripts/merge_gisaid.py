"""
Script Name: merge_gisaid.py
Purpose:     Merge the GISAID downloads of one virus (several TSVs and FASTAs in one folder) into
             one metadata table and one FASTA. Used for NiV and Mpox; the RSV batches have their
             own script (merge_gisaid_batches.py) because they sit in one folder per batch.
Inputs:      data/Raw/<virus>/**/*.tsv and *.fasta, subfolders allowed (one per search) (GISAID default
             names, which are not reliable)
Outputs:     data/processed/<virus>/<virus>_metadata.tsv, <virus>.fasta, merge_report.txt
Notes:       Join key = Accession ID; for each field the first non-empty value wins. FASTA header
             is >name|accession|date, so the accession is the 2nd '|' field. Duplicate accessions
             across files are counted and kept once.
             Run: .venv/bin/python scripts/merge_gisaid.py NiV   (or Mpox)
"""
import csv, glob, os, sys

virus = sys.argv[1]
raw, out = f"data/Raw/{virus}", f"data/processed/{virus}"
os.makedirs(out, exist_ok=True)
csv.field_size_limit(10**8)

meta, n_rows = {}, 0
for f in sorted(glob.glob(f"{raw}/**/*.tsv", recursive=True)):
    with open(f, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n_rows += 1
            m = meta.setdefault(row["Accession ID"], {})
            for k, v in row.items():
                if k and v not in (None, "") and k not in m:
                    m[k] = v

fasta, n_dup = {}, 0
for f in sorted(glob.glob(f"{raw}/**/*.fasta", recursive=True)):
    h, s = None, []
    with open(f) as fh:
        for line in list(fh) + [">END|END|END"]:
            line = line.strip()
            if line.startswith(">"):
                if h:
                    acc = h.split("|")[1]
                    n_dup += acc in fasta
                    fasta[acc] = (h, "".join(s).upper())
                h, s = line[1:], []
            else:
                s.append(line)

cols = []
for r in meta.values():
    for k in r:
        if k not in cols:
            cols.append(k)
with open(f"{out}/{virus}_metadata.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(meta.values())
with open(f"{out}/{virus}.fasta", "w") as fh:
    for acc, (h, s) in fasta.items():
        fh.write(f">{h}\n{s}\n")

with open(f"{out}/merge_report.txt", "w") as fh:
    fh.write(f"tsv rows read: {n_rows}\nunique accessions in metadata: {len(meta)}\n"
             f"unique sequences: {len(fasta)} (duplicates across files: {n_dup})\n"
             f"in metadata but no sequence: {len(set(meta) - set(fasta))}\n"
             f"in fasta but no metadata: {len(set(fasta) - set(meta))}\n")
print(open(f"{out}/merge_report.txt").read())
