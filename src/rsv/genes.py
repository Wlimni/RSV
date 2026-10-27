"""
Script Name: genes.py
Purpose:     Extract individual RSV genes from whole genomes by aligning to an annotated
             reference, so composition can be measured per gene (as the paper did per segment).
Inputs:      A reference GenBank file (CDS coordinates) + reference FASTA, and query sequences.
Outputs:     dict of {gene: {accession: gene_sequence}} and tidy per-gene feature tables.
Notes:       Alignment keeps insertions (no --keeplength), because RSV-A ON1 and RSV-B BA
             genotypes carry well-known duplications in the G gene; dropping them would
             bias exactly the gene we care most about.
"""
import os, re, subprocess, tempfile
from collections import OrderedDict


def read_cds(gb_path):
    """Parse gene name -> (start, end) 1-based inclusive CDS coordinates from a GenBank file."""
    genes = OrderedDict()
    cur = None
    for line in open(gb_path):
        m = re.match(r"^     CDS\s+(?:join\()?(<?\d+)\.\.(>?\d+)", line)
        if m:
            cur = (int(m.group(1).lstrip("<")), int(m.group(2).lstrip(">")))
            continue
        if cur:
            g = re.search(r'/gene="([^"]+)"', line)
            if g:
                genes.setdefault(g.group(1), cur)
                cur = None
    return genes


def read_fasta(path):
    seqs, name, buf = OrderedDict(), None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def align_to_reference(ref_fasta, query_seqs, threads=-1, tmpdir=None):
    """
    Align query sequences against a single reference with mafft --addfragments.
    Returns (aligned_reference_row, {accession: aligned_row}).
    Insertions relative to the reference appear as extra alignment columns.
    """
    with tempfile.TemporaryDirectory(dir=tmpdir) as td:
        qpath = os.path.join(td, "q.fasta")
        with open(qpath, "w") as fh:
            for acc, s in query_seqs.items():
                fh.write(f">{acc}\n{s}\n")
        cmd = ["mafft", "--quiet", "--6merpair", "--thread", str(threads),
               "--addfragments", qpath, ref_fasta]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    rows, name, buf = OrderedDict(), None, []
    for line in out.split("\n"):
        if line.startswith(">"):
            if name:
                rows[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        rows[name] = "".join(buf)
    ref_name = next(iter(rows))
    ref_row = rows.pop(ref_name)
    return ref_row, rows


def ref_pos_to_column(ref_row):
    """Map 1-based reference position -> alignment column index."""
    pos2col, p = {}, 0
    for col, ch in enumerate(ref_row):
        if ch != "-":
            p += 1
            pos2col[p] = col
    return pos2col


def slice_genes(ref_row, rows, cds):
    """
    Cut each aligned sequence into its genes using reference CDS coordinates.
    Gaps are stripped, so each returned gene sequence is the query's own nucleotides
    (including any insertion it carries inside that gene).
    """
    pos2col = ref_pos_to_column(ref_row)
    last = max(pos2col)
    out = {}
    for gene, (s, e) in cds.items():
        s, e = max(1, s), min(last, e)
        c0, c1 = pos2col[s], pos2col[e]
        out[gene] = {acc: row[c0:c1 + 1].replace("-", "").upper()
                     for acc, row in rows.items()}
    return out
