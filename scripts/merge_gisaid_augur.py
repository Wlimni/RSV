"""
Script Name: merge_gisaid_augur.py
Purpose:     Merge GISAID "Input for the Augur pipeline" downloads (one .tar per search, each holding a
             metadata TSV and a FASTA) into one metadata table and one FASTA.
Inputs:      data/Raw/<virus>/<search folder>/*.tar
Outputs:     data/processed/<virus>/<virus>_metadata.tsv   (adds a `search` column = folder name)
             data/processed/<virus>/<virus>.fasta          (header = >EPI_ISL accession)
             data/processed/<virus>/merge_report.txt
Notes:       The tar members have absolute paths (/123.metadata.tsv), so files are read straight from
             the archive and nothing is extracted to disk. FASTA headers are strain names; they are
             mapped to the EPI_ISL accession through the metadata `strain` column. A genome that
             appears in several searches (for example two labels for the same host) is kept once.
             Run: .venv/bin/python scripts/merge_gisaid_augur.py Cov_19
"""
import glob, io, os, sys, tarfile
import pandas as pd

virus = sys.argv[1]
raw, out = f"data/Raw/{virus}", f"data/processed/{virus}"
os.makedirs(out, exist_ok=True)

meta, fasta, report = [], {}, []
for tpath in sorted(glob.glob(f"{raw}/*/*.tar")):
    search = os.path.basename(os.path.dirname(tpath))
    with tarfile.open(tpath) as tf:
        names = tf.getnames()
        mname = next(n for n in names if n.endswith("metadata.tsv"))
        fname = next(n for n in names if n.endswith(".fasta"))
        m = pd.read_csv(io.BytesIO(tf.extractfile(mname).read()), sep="\t", dtype=str)
        m["search"] = search
        strain2acc = dict(zip(m.strain, m.gisaid_epi_isl))
        name, buf = None, []
        n_seq = 0
        for line in io.TextIOWrapper(tf.extractfile(fname), encoding="utf-8"):
            line = line.strip()
            if line.startswith(">"):
                if name:
                    fasta[strain2acc.get(name, name)] = "".join(buf).upper(); n_seq += 1
                name, buf = line[1:], []
            else:
                buf.append(line)
        if name:
            fasta[strain2acc.get(name, name)] = "".join(buf).upper(); n_seq += 1
    report.append(f"{search}: {len(m)} metadata rows, {n_seq} sequences")
    meta.append(m)

meta = pd.concat(meta, ignore_index=True)
n_dup = meta.gisaid_epi_isl.duplicated().sum()
meta = meta.drop_duplicates("gisaid_epi_isl")
meta.to_csv(f"{out}/{virus}_metadata.tsv", sep="\t", index=False)
with open(f"{out}/{virus}.fasta", "w") as fh:
    for acc in meta.gisaid_epi_isl:
        if acc in fasta:
            fh.write(f">{acc}\n{fasta[acc]}\n")
report += [f"unique genomes: {len(meta)} (duplicates across searches: {n_dup})",
           f"in metadata but no sequence: {len(set(meta.gisaid_epi_isl) - set(fasta))}",
           f"in fasta but no metadata: {len(set(fasta) - set(meta.gisaid_epi_isl))}"]
open(f"{out}/merge_report.txt", "w").write("\n".join(report) + "\n")
print("\n".join(report))
