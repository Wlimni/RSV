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
from collections import Counter, OrderedDict


def read_cds(gb_path):
    """
    Parse CDS features into {name: (start, end)}, 1-based inclusive, in file order.

    A gene can carry more than one CDS. RSV's M2 is the case that matters: it holds
    M2-1 (ORF 1) and M2-2 (ORF 2), which overlap by 35 nt. Keying on /gene alone would
    silently keep whichever came first and label it "M2", so where a /gene name is used
    by several CDS, each is named from its /product instead ("M2-1 protein" -> "M2-1").

    A join(...) location is reduced to its outermost coordinates. RSV CDS are contiguous,
    so this is exact here; for a spliced CDS it would span the intron, and the caller
    would need a real feature parser.
    """
    recs, loc, gene, product = [], None, None, None

    def flush():
        if loc is not None:
            recs.append({"loc": loc, "gene": gene, "product": product})

    for line in open(gb_path):
        m = re.match(r"^     CDS\s+(.*)$", line)
        if m:
            flush()
            nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
            loc = (min(nums), max(nums)) if nums else None
            gene = product = None
            continue
        if loc is not None:
            if re.match(r"^     \S", line):      # next feature began
                flush()
                loc = gene = product = None
                continue
            g = re.search(r'/gene="([^"]+)"', line)
            if g:
                gene = g.group(1)
            q = re.search(r'/product="([^"]+)"', line)
            if q:
                product = q.group(1)
    flush()

    counts = Counter(r["gene"] for r in recs if r["gene"])
    out = OrderedDict()
    for r in recs:
        if r["gene"] is None or r["loc"] is None:
            continue
        name = r["gene"]
        if counts[name] > 1 and r["product"]:
            # "M2-1 protein" -> "M2-1"; falls back to the gene name if the product is bare
            name = r["product"].split()[0] or name
        out.setdefault(name, r["loc"])
    return out


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
