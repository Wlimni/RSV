"""
Script Name: build_features.py
Purpose:     Produce the analysis-ready feature table: one row per QC-passing sequence with
             metadata (subtype, lineage, date, region) joined to composition features
             (GC content, 4 dinucleotide freqs, 4 O/E ratios).
Inputs:      data/full/RSV_full.fasta, data/full/RSV_full_metadata.tsv
Outputs:     data/processed/rsv_features.tsv, data/processed/qc_report.txt
Notes:       Composition is computed on the whole genome here; per-gene features are a later
             step (needs alignment/annotation). Called from notebooks/01_data.ipynb; can also
             be run standalone: python3 -I scripts/build_features.py
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd
from rsv import io, composition

FULL = "data/full"
OUT = "data/processed"


def main():
    os.makedirs(OUT, exist_ok=True)

    seqs = io.read_fasta(os.path.join(FULL, "RSV_full.fasta"))
    meta = io.load_metadata(os.path.join(FULL, "RSV_full_metadata.tsv"))
    clean, report = io.qc_filter(meta, seqs)

    # composition features for each QC-passing sequence
    feats = pd.DataFrame(
        [composition.all_features(seqs[a]) for a in clean["Accession ID"]],
        index=clean.index,
    )
    keep = ["Accession ID", "Subtype", "Lineage", "Genotype", "date", "year",
            "region", "country", "seq_len_actual", "ambig_frac"]
    out = pd.concat([clean[keep], feats], axis=1)
    out.to_csv(os.path.join(OUT, "rsv_features.tsv"), sep="\t", index=False)

    lines = ["QC filter counts:"] + [f"  {k}: {v}" for k, v in report.items()]
    lines += [f"\nfeature rows written: {len(out)}",
              f"columns: {list(out.columns)}",
              f"\nby subtype:\n{out.groupby('Subtype').size().to_string()}"]
    open(os.path.join(OUT, "qc_report.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
