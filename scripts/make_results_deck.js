/*
 * Script Name: make_results_deck.js
 * Purpose:     Build the results deck for the supervisor: a one-slide method summary, the six
 *              result figures each with a one-sentence takeaway and a short explanation, and
 *              next steps.
 * Inputs:      results/figures/fig1..fig6 PNGs (fig3 is the O/E-only version, because the full
 *              figure is too tall for a 16:9 slide).
 * Outputs:     docs/RSV_Results_Update.pptx
 * Notes:       Run from the repository root:
 *                NODE_PATH=<dir holding pptxgenjs> node scripts/make_results_deck.js
 *              SKILL_DIR points at the pptx skill, whose apply_theme.js writes the theme colours
 *              (pptxgenjs cannot). Every number on a slide is checked against results/*.tsv by
 *              hand; if the analysis is rerun, re-read the slides before sending.
 */
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const SKILL_DIR = process.env.SKILL_DIR;
if (!SKILL_DIR) throw new Error("set SKILL_DIR to the pptx skill directory (it holds apply_theme.js)");
const { applyTheme } = require(path.join(SKILL_DIR, "scripts", "apply_theme.js"));

const OUT = "docs/RSV_Results_Update.pptx";
const FIG = "results/figures";

// Blue and red are reserved for RSV-A and RSV-B, as in every figure; the slide chrome uses
// neutral ink and teal so it never competes with them.
const THEME = {
  name: "RSV results",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1B2631", lt1: "FFFFFF", dk2: "34495E", lt2: "F2F4F6",
    accent1: "2C7FB8", accent2: "C62824", accent3: "0F8B8D",
    accent4: "D98E04", accent5: "6B7B8C", accent6: "C5CED7",
    hlink: "0F8B8D", folHlink: "6B7B8C",
  },
};

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";                       // 13.33 x 7.5 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "RSV composition: results update";
pres.author = "William Ni Hong Hung";
const C = pres.SchemeColor;

// ---- layouts -----------------------------------------------------------------------------
pres.defineSlideMaster({
  title: "TITLE",
  background: { color: C.text1 },
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: C.background1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.4, w: 12.13, h: 1.05,
        fontFace: "Cambria", fontSize: 24, bold: true, color: C.text1, align: "left",
        valign: "top", margin: 0 }, text: "" } },
    { text: { text: "RSV composition  ·  results update", options: { x: 0.6, y: 7.0, w: 5, h: 0.3,
        fontSize: 10, color: C.accent5, margin: 0 } } },
  ],
  slideNumber: { x: 12.18, y: 7.0, w: 0.55, h: 0.3, fontSize: 10, color: C.accent5, align: "right" },
});

// ---- helpers -----------------------------------------------------------------------------
// PNG width and height sit at fixed byte offsets in the header, which avoids an image library
function pngSize(file) {
  const b = fs.readFileSync(file);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}

// small A/B key in the footer of every figure slide
function legend(slide) {
  slide.addShape(pres.ShapeType.ellipse, { x: 10.55, y: 7.08, w: 0.14, h: 0.14,
    fill: { color: C.accent1 }, line: { color: C.accent1 }, objectName: "key A" });
  slide.addText("RSV-A", { x: 10.75, y: 7.0, w: 0.7, h: 0.3, fontSize: 10, color: C.accent5,
    margin: 0, isTextBox: true, objectName: "key A label" });
  slide.addShape(pres.ShapeType.ellipse, { x: 11.5, y: 7.08, w: 0.14, h: 0.14,
    fill: { color: C.accent2 }, line: { color: C.accent2 }, objectName: "key B" });
  slide.addText("RSV-B", { x: 11.7, y: 7.0, w: 0.7, h: 0.3, fontSize: 10, color: C.accent5,
    margin: 0, isTextBox: true, objectName: "key B label" });
}

const label = (text, extra = {}) => ({ text, options: Object.assign(
  { fontSize: 11, bold: true, color: C.accent3, charSpacing: 2, breakLine: true, paraSpaceAfter: 4 }, extra) });
const para = (text, extra = {}) => ({ text, options: Object.assign(
  { fontSize: 15, color: C.text1, breakLine: true, paraSpaceAfter: 6 }, extra) });
const bullet = (text, last = false) => para(text, { bullet: { indent: 14 }, breakLine: !last });

const BOX = { x: 0.6, y: 1.7, w: 8.0, h: 5.2 };       // where the figure may sit
const RIGHT = 12.73;                                   // right content edge

// figure on the left in a thin frame, explanatory panel on the right
function figureSlide({ section, title, img, panel, notes }) {
  const slide = pres.addSlide({ masterName: "CONTENT", sectionTitle: section });
  slide.addText(title, { placeholder: "title" });

  const px = pngSize(path.join(FIG, img));
  const aspect = px.w / px.h;
  let w = BOX.w, h = w / aspect;
  if (h > BOX.h) { h = BOX.h; w = h * aspect; }
  const y = BOX.y + (BOX.h - h) / 2;

  slide.addShape(pres.ShapeType.rect, { x: BOX.x, y, w, h, fill: { color: C.background1 },
    line: { color: C.accent6, width: 0.75 }, objectName: "figure frame" });
  slide.addImage({ path: path.join(FIG, img), x: BOX.x, y, w, h, altText: title });

  const x0 = BOX.x + w + 0.35, pw = RIGHT - x0;
  slide.addShape(pres.ShapeType.roundRect, { x: x0, y: BOX.y, w: pw, h: BOX.h, rectRadius: 0.08,
    fill: { color: C.background2 }, line: { color: C.background2 }, objectName: "explanation panel" });
  panel(slide, x0, pw);
  legend(slide);
  slide.addNotes(notes);
  return slide;
}

// the common panel: what it shows (bullets) and one comparison line
function standardPanel(shows, vsLabel, vs) {
  return (slide, x0, pw) => {
    const runs = [label("WHAT IT SHOWS")];
    shows.forEach((s, i) => runs.push(bullet(s, false)));
    runs.push(label(vsLabel, { paraSpaceBefore: 10 }));
    runs.push(para(vs, { breakLine: false, paraSpaceAfter: 0 }));
    slide.addText(runs, { x: x0 + 0.25, y: BOX.y + 0.22, w: pw - 0.5, h: BOX.h - 0.44,
      valign: "top", margin: 0, isTextBox: true, objectName: "explanation text" });
  };
}

// ============================ 1. title ======================================================
pres.addSection({ title: "Overview" });
{
  const s = pres.addSlide({ masterName: "TITLE", sectionTitle: "Overview" });
  // the two subtypes as overlapping circles: the one subject every figure is about
  s.addShape(pres.ShapeType.ellipse, { x: 8.7, y: 2.0, w: 3.0, h: 3.0, fill: { color: C.accent1, transparency: 12 },
    line: { color: C.accent1 }, objectName: "RSV-A circle" });
  s.addShape(pres.ShapeType.ellipse, { x: 9.95, y: 2.8, w: 3.0, h: 3.0, fill: { color: C.accent2, transparency: 12 },
    line: { color: C.accent2 }, objectName: "RSV-B circle" });
  s.addText("RSV-A", { x: 9.0, y: 2.7, w: 1.4, h: 0.5, fontSize: 18, bold: true, color: C.background1,
    isTextBox: true, margin: 0, objectName: "RSV-A label" });
  s.addText("RSV-B", { x: 11.15, y: 4.4, w: 1.4, h: 0.5, fontSize: 18, bold: true, color: C.background1,
    align: "right", isTextBox: true, margin: 0, objectName: "RSV-B label" });

  s.addText("Does the influenza GC and CpG finding reproduce in human RSV?", { x: 0.9, y: 1.7, w: 7.7, h: 2.3,
    fontFace: "Cambria", fontSize: 38, bold: true, color: C.background1, valign: "top", margin: 0,
    isTextBox: true, objectName: "deck title" });
  s.addText("Results update for Yongtao Ye", { x: 0.9, y: 4.25, w: 7.7, h: 0.5, fontSize: 20,
    color: C.accent6, margin: 0, isTextBox: true, objectName: "subtitle" });
  s.addText([
    { text: "William Ni Hong Hung  ·  2 November 2026", options: { breakLine: true } },
    { text: "Method adapted from Ye et al. 2026, Nature Microbiology" },
  ], { x: 0.9, y: 5.4, w: 7.7, h: 0.8, fontSize: 14, color: C.accent6, margin: 0, isTextBox: true,
    paraSpaceAfter: 4, objectName: "byline" });
  s.addNotes("The question from the paper: GC content and CpG fall over time in influenza lineages that have " +
    "moved into a new mammalian host. I applied the same measurements to 32,442 human RSV genomes. " +
    "Short version: the time-trend result is negative, and the interesting finding is a CpG pattern in " +
    "NS1 and NS2 that differs between the two subtypes.");
}

// ============================ 2. method (the wordy one) =====================================
{
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Overview" });
  s.addText("We measured composition in 32,442 RSV genomes and fitted the trend inside each lineage, as the influenza paper did",
    { placeholder: "title" });

  // funnel: how many genomes survive each stage
  const chips = [
    ["33,040", "genomes pooled from GISAID, duplicates removed"],
    ["33,038", "pass quality control (2 removed)"],
    ["32,442", "inside the window and the 9 kept lineages"],
    ["9", "lineages: 6 RSV-A and 3 RSV-B"],
  ];
  const gap = 0.5, cw = (12.13 - 3 * gap) / 4;
  chips.forEach(([big, small], i) => {
    const x = 0.6 + i * (cw + gap);
    s.addShape(pres.ShapeType.roundRect, { x, y: 1.6, w: cw, h: 1.0, rectRadius: 0.08,
      fill: { color: C.background2 }, line: { color: C.background2 }, objectName: `funnel box ${i + 1}` });
    s.addText([
      { text: big, options: { fontSize: 26, bold: true, color: C.accent3, breakLine: true } },
      { text: small, options: { fontSize: 12, color: C.text2 } },
    ], { x: x + 0.15, y: 1.7, w: cw - 0.3, h: 0.85, valign: "top", margin: 0, isTextBox: true,
      objectName: `funnel text ${i + 1}` });
    if (i < 3) {
      s.addShape(pres.ShapeType.rightArrow, { x: x + cw + 0.1, y: 1.99, w: 0.3, h: 0.24,
        fill: { color: C.accent5 }, line: { color: C.accent5 }, objectName: `funnel arrow ${i + 1}` });
    }
  });

  // the five steps as a table: what was done, and the rule or number behind it
  const head = (t) => ({ text: t, options: { bold: true, color: C.background1, fill: { color: C.text2 }, fontSize: 13 } });
  const step = (n, name) => ({ text: [
    { text: `${n}   `, options: { bold: true, color: C.accent3 } },
    { text: name, options: { bold: true } }] });
  const rows = [
    [head("Step"), head("What we did"), head("Rule or number")],
    [step(1, "Data"),
     "GISAID whole genomes of RSV-A and RSV-B collected 2001 to 2026, with subtype, lineage, date and region",
     "15,458 RSV-A  ·  17,580 RSV-B"],
    [step(2, "Quality control"),
     "Keep complete genomes with a collection date; drop contamination and poor reads",
     "Length 14,900 to 15,400 nt  ·  at most 1% ambiguous bases"],
    [step(3, "Window and lineages"),
     "Cut 2026 at the end of May (later months are deposition lag); trim sparse lineage ends; keep lineages with at least 5 years and 200 genomes",
     "May holds 92% of 2025's count, June 53%"],
    [step(4, "Genes and measures"),
     "Align to the subtype reference (MAFFT), cut 10 genes (M2-1 only); measure GC, GC3 and GC12, and the O/E ratio of CpG, GpC, GpG, CpC and UpA",
     "O/E = observed divided by expected from base content"],
    [step(5, "Statistics"),
     "Slope of each measure on collection date within each lineage; standard errors clustered on sub-lineage because genomes are related",
     "Subtype contrast judged across lineages (6 A vs 4 B), not genomes"],
  ];
  const tableRows = rows.map((r, ri) => r.map((cell, ci) => {
    if (ri === 0) return cell;
    const o = typeof cell === "string" ? { text: cell, options: {} } : cell;
    o.options = Object.assign({ fill: { color: ri % 2 ? C.background1 : C.background2 } }, o.options || {});
    return o;
  }));
  s.addTable(tableRows, { x: 0.6, y: 2.8, w: 12.13, colW: [2.45, 6.2, 3.48],
    rowH: [0.36, 0.74, 0.74, 0.74, 0.74, 0.74], fontFace: "Calibri", fontSize: 14, color: C.text1,
    valign: "middle", border: { type: "solid", pt: 0.5, color: "C5CED7" }, margin: [0.04, 0.12, 0.04, 0.12],
    objectName: "method table" });
  s.addNotes("QC detail: the two genomes removed were an RSV-B with 2.11% ambiguous bases, and an RSV-B of 15,854 nt " +
    "carrying a 583-nt non-RSV tail at 49% GC against the virus's 34%, which is contamination. That second case is " +
    "why the length filter has an upper limit as well as a lower one.\n\n" +
    "Window: the download was made in October 2026 and sequences take months to be deposited. Compared with the " +
    "same months of 2025, 2026 is at 78 to 93% through May, then 53% in June, 29% in July and 4% in August, so the " +
    "months after May are a non-random subset and are cut.\n\n" +
    "Statistics: genomes inside a lineage share ancestry, so plain regression treats one expanding clade as hundreds " +
    "of independent points. Clustering the standard errors on sub-lineage inflates them about 6.6-fold (median). " +
    "Inference uses a t distribution with clusters minus one degrees of freedom, since a lineage often has only a " +
    "handful of clusters.\n\n" +
    "Codon positions (GC3 and GC12) are measured gene by gene, in frame. An earlier version measured them on the " +
    "concatenated genes, which put every gene downstream of an odd-length gene out of frame.");
}

// ============================ 3 to 8. the six figures =======================================
pres.addSection({ title: "Results" });

figureSlide({
  section: "Results",
  title: "GC content has no consistent trend: 4 of 9 lineages are significant, and they point both ways",
  img: "fig1_lineage_trends.png",
  panel: standardPanel([
    "One panel per lineage: each dot is a genome, the line is the fitted trend.",
    "The slopes are tiny, at most 0.012% GC per year.",
    "p values count related genomes as one cluster; plain regression would call 7 of 9 significant.",
  ], "VERSUS THE PAPER", "All 12 influenza lineages fell (P < 0.01). RSV does not."),
  notes: "Significant lineages after clustering: A.D.2 (up), A.D.4 (down), B.D.4 (up) and B.D.E (up). Five slopes are " +
    "positive and four negative. The largest is A.D.5 at minus 0.0117% per year.\n\n" +
    "Caution on A.D.4: it has only 368 genomes and falls back to clustering by calendar year (10 clusters), which " +
    "makes it more significant than plain regression, not less. Treat that single p value (0.0014) with care.\n\n" +
    "The paper's lineages had all been sampled across a host shift. Human RSV has circulated in people for a very " +
    "long time, so the absence of a directional trend is a reasonable biological result rather than a failed analysis.",
});

figureSlide({
  section: "Results",
  title: "GC3 and GC12 show no reliable contrast, unlike influenza where GC3 fell 5.65x faster",
  img: "fig2_gc3_vs_gc12.png",
  panel: standardPanel([
    "Change per year at all sites, first and second codon sites (GC12) and third sites (GC3).",
    "GC3 rises in 7 of 9 lineages; GC12 is split, 4 down and 5 up.",
    "Only 3 of 9 GC3 and 1 of 9 GC12 slopes are significant.",
  ], "VERSUS THE PAPER", "Both fell, GC3 faster, which points at pressure on the RNA. RSV shows no such signal."),
  notes: "An earlier version of this analysis reported GC3 rising in 6 lineages and GC12 falling in 6, as an " +
    "opposite-directions result. That came partly from a bug: codon positions were measured on the concatenated " +
    "genes, which is out of frame downstream of any gene whose length is not a multiple of 3 (about 4% of genomes). " +
    "Correcting it changed the slopes by up to 4.7e-4 per year and flipped 3 of 9 GC3 signs, so the misframing was " +
    "correlated with collection date.\n\n" +
    "The corrected pattern is itself fragile: refitting from only the seven genes that rarely truncate agrees in " +
    "sign for just 7 of 9 lineages (GC3) and 6 of 9 (GC12). So I report it as unresolved, not as a finding.\n\n" +
    "Magnitude: GC3 moves 2.65x as fast as GC12 by mean absolute slope, against 5.65x in the paper. A ratio of " +
    "signed sums is not meaningful here because the slopes carry both signs.",
});

figureSlide({
  section: "Results",
  title: "CpG O/E shows no consistent trend either: only 1 of 9 lineages is significant",
  img: "fig3_oe_only.png",
  panel: standardPanel([
    "Observed over expected CpG per year, plus three controls, corrected for base content.",
    "CpG rises in 6 lineages and falls in 3; only A.D.5 is significant, and it is falling.",
    "At a finer lineage level the median direction flips (15 of 20 negative), so direction depends on how lineages are drawn.",
  ], "VERSUS THE PAPER", "CpG O/E fell in 9 of 12 influenza lineages."),
  notes: "UpA is shown here for contrast but is not in the paper. GpC, GpG and CpC are the paper's comparison " +
    "dinucleotides.\n\n" +
    "The single significant lineage, A.D.5, is declining, which is the direction the paper reports. One lineage " +
    "out of nine is not a pattern, and it does not survive the finer lineage definition: at level 3 the median CpG " +
    "slope is slightly positive (3 of 9 negative); at level 4 it is negative (15 of 20 negative). When direction " +
    "depends on granularity, granularity is driving it rather than biology.\n\n" +
    "Slopes are about 1e-3 O/E units per year against a CpG O/E of about 0.24.",
});

figureSlide({
  section: "Results",
  title: "CpG is not just GC moving, but it is less tied to GC than in influenza",
  img: "fig4_gc_correlation.png",
  panel: standardPanel([
    "Correlation of each dinucleotide's frequency with GC content, within each lineage.",
    "Near 1 would mean it adds nothing beyond GC. CpG averages 0.43, with a range of 0.25 to 0.62.",
    "UpA is left out: U and A are the bases GC does not count, so its correlation is negative by definition.",
  ], "VERSUS THE PAPER", "Correlations of 0.59 to 0.73. RSV's dinucleotides average 0.35 to 0.51."),
  notes: "Means over lineages: CpC 0.51, GpC 0.49, CpG 0.43, GpG 0.35. The paper reports 0.72, 0.67, 0.73 and 0.59 " +
    "for its four dinucleotides.\n\n" +
    "I originally computed the UpA correlation too, got about minus 0.65 in every lineage, and removed it: since " +
    "A and U are the complement of G and C, the correlation is forced negative and carries no information about selection.",
});

figureSlide({
  section: "Results",
  title: "CpG suppression is uneven along the genome: NS1 and NS2 sit lowest in RSV-A",
  img: "fig5_gene_composition.png",
  panel: standardPanel([
    "Mean GC content (top) and CpG O/E (bottom) for each of the 10 genes, 3' to 5'.",
    "GC content varies by gene but is nearly the same in RSV-A and RSV-B.",
    "RSV-A CpG O/E is 0.016 in NS1 and 0.072 in NS2, against 0.16 to 0.41 in every other gene.",
  ], "WHY IT MATTERS", "NS1 and NS2 are the genes that switch off the interferon alarm."),
  notes: "M2 is M2-1 only. The reference annotates two overlapping reading frames (M2-1 and M2-2, overlapping by 32 to 35 " +
    "nt depending on subtype), and counting both would count that stretch twice.\n\n" +
    "The G gene is the most GC-rich (43% in RSV-A, 40% in RSV-B) and SH is among the lowest.\n\n" +
    "Our search of PubMed and the web found no published comparison of CpG between RSV-A and RSV-B, and no ZAP study " +
    "in RSV or any pneumovirus. Engineered work exists (Liang 2017, Le Nouen 2020, Meng 2014) but not this observation.",
});

// fig6 gets stat callouts instead of bullets: it is the main finding
figureSlide({
  section: "Results",
  title: "RSV-A has almost no CpG in NS1 and NS2, and the gap holds across independent lineages",
  img: "fig6_slide.png",
  panel: (slide, x0, pw) => {
    const stats = [
      ["12.7x", "more CpG (O/E) in RSV-B than RSV-A in NS1; NS2 is 7.2x"],
      ["0.001–0.027  vs  0.163–0.210", "NS1 lineage means, RSV-A vs RSV-B: the ranges do not overlap"],
      ["30.5%  vs  30.3%", "GC content of NS1 in A and B: the base composition is the same"],
    ];
    stats.forEach(([big, small], i) => {
      const y = BOX.y + 0.2 + i * 1.05;
      slide.addShape(pres.ShapeType.roundRect, { x: x0 + 0.2, y, w: pw - 0.4, h: 0.95, rectRadius: 0.06,
        fill: { color: C.background1 }, line: { color: C.accent6, width: 0.75 }, objectName: `stat card ${i + 1}` });
      slide.addText([
        { text: big, options: { fontSize: 20, bold: true, color: C.accent3, breakLine: true } },
        { text: small, options: { fontSize: 14, color: C.text1 } },
      ], { x: x0 + 0.4, y: y + 0.03, w: pw - 0.8, h: 0.89, valign: "middle", margin: 0, isTextBox: true,
        objectName: `stat text ${i + 1}` });
    });
    slide.addText([
      label("HOW TO READ IT"),
      para("Top: CpG O/E, bars are genome means and whiskers span the lineage means. Bottom: GC content, the control.", { paraSpaceAfter: 8 }),
      para("6 of 10 genes separate across lineages (NS1, NS2, P, M, SH, M2-1); the whiskers of N, G, F and L overlap.", { breakLine: false, paraSpaceAfter: 0 }),
    ], { x: x0 + 0.25, y: BOX.y + 3.4, w: pw - 0.5, h: 1.7, valign: "top", margin: 0, isTextBox: true,
      objectName: "how to read" });
  },
  notes: "The slide shows two of the three panels so the type stays readable. The full figure, which adds CpG counts per gene with the share of the gene on each bar, is results/figures/fig6_subtype_contrast.png.\n\n" +
    "Genome-level p values here are below 1e-300 for every gene, which is pseudoreplication: genomes in a lineage " +
    "share ancestry, and RSV-A and RSV-B are two clades, so the per-genome comparison is closer to a single " +
    "evolutionary contrast than to 33,000 samples. The claim therefore rests on lineage means (6 RSV-A against 4 " +
    "RSV-B lineages). A Mann-Whitney test on 6 against 4 bottoms out at p = 0.0095 even with perfect separation, so " +
    "I quote the ranges and not a p value.\n\n" +
    "Genome-level detail: 84.8% of RSV-A genomes carry no CpG anywhere in NS1, and 34.2% none in NS2. The earlier " +
    "claim that CpG differs between subtypes in every gene was an artefact of counting genomes; at lineage level only " +
    "six genes separate.\n\n" +
    "Is it CpG-specific? RSV's CpG O/E is 0.24 in both subtypes while UpA is 0.90 and GpC, GpG and CpC sit near 1. " +
    "UpA is the one other motif with experimental evidence of restricting RNA virus replication, through a separate " +
    "recognition site (Goonawardane et al. 2021). If RSV were under broad pressure against immune-visible motifs, UpA " +
    "would be suppressed too. It is not.\n\n" +
    "This is an association in sequence, not a demonstrated mechanism. Testing ZAP would need laboratory work.",
});

// ============================ 9. conclusions and next steps ================================
pres.addSection({ title: "Next" });
{
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: "Next" });
  s.addText("The trend result is negative; the NS1 and NS2 CpG loss is the finding to pursue", { placeholder: "title" });

  const colW = 5.95, xL = 0.6, xR = 6.78;
  s.addText("WHAT WE FOUND", { x: xL, y: 1.6, w: colW, h: 0.3, fontSize: 12, bold: true, color: C.accent3,
    charSpacing: 2, margin: 0, isTextBox: true, objectName: "found label" });
  s.addText("NEXT STEPS", { x: xR, y: 1.6, w: colW, h: 0.3, fontSize: 12, bold: true, color: C.accent3,
    charSpacing: 2, margin: 0, isTextBox: true, objectName: "next label" });

  const found = [
    ["1 to 4 of 9", "No consistent trend", "Lineages are significant in only 1 to 4 of 9 for GC, GC3, GC12 and CpG, and point both ways."],
    ["12.7x", "NS1 and NS2 in RSV-A", "CpG O/E far below RSV-B, and it holds across independent lineages."],
    ["0.24 vs 0.90", "CpG-specific pressure", "CpG O/E against UpA O/E: UpA is barely suppressed."],
  ];
  found.forEach(([stat, head, body], i) => {
    const y = 2.0 + i * 1.3;
    s.addShape(pres.ShapeType.roundRect, { x: xL, y, w: colW, h: 1.18, rectRadius: 0.06,
      fill: { color: C.background2 }, line: { color: C.background2 }, objectName: `found card ${i + 1}` });
    s.addText(stat, { x: xL + 0.2, y, w: 1.85, h: 1.18, fontSize: 20, bold: true, color: C.accent3,
      valign: "middle", margin: 0, isTextBox: true, objectName: `found stat ${i + 1}` });
    s.addText([
      { text: head, options: { fontSize: 16, bold: true, color: C.text1, breakLine: true } },
      { text: body, options: { fontSize: 14, color: C.text1 } },
    ], { x: xL + 2.1, y, w: colW - 2.3, h: 1.18, valign: "middle", margin: 0, isTextBox: true,
      objectName: `found text ${i + 1}` });
  });

  const next = [
    ["Is the loss at synonymous sites?", "Silent changes at third positions would point to selection on the RNA rather than the protein."],
    ["Is RSV-A NS1 CpG still declining?", "Within lineages over time: still falling, or already at the floor. This is the alternative explanation to rule out."],
    ["Add a paramyxovirus panel?", "Nipah and Hendra (sporadic) against human RSV, for the sustained-versus-sporadic contrast. Needs your steer."],
  ];
  next.forEach(([head, body], i) => {
    const y = 2.0 + i * 1.3;
    s.addShape(pres.ShapeType.roundRect, { x: xR, y, w: colW, h: 1.18, rectRadius: 0.06,
      fill: { color: C.background1 }, line: { color: C.accent6, width: 0.75 }, objectName: `next card ${i + 1}` });
    s.addShape(pres.ShapeType.ellipse, { x: xR + 0.2, y: y + 0.36, w: 0.46, h: 0.46,
      fill: { color: C.accent3 }, line: { color: C.accent3 }, objectName: `next badge ${i + 1}` });
    s.addText(String(i + 1), { x: xR + 0.2, y: y + 0.36, w: 0.46, h: 0.46, fontSize: 16, bold: true,
      color: C.background1, align: "center", valign: "middle", margin: 0, isTextBox: true,
      objectName: `next number ${i + 1}` });
    s.addText([
      { text: head, options: { fontSize: 16, bold: true, color: C.text1, breakLine: true } },
      { text: body, options: { fontSize: 14, color: C.text1 } },
    ], { x: xR + 0.85, y, w: colW - 1.05, h: 1.18, valign: "middle", margin: 0, isTextBox: true,
      objectName: `next text ${i + 1}` });
  });

  s.addShape(pres.ShapeType.roundRect, { x: 0.6, y: 5.98, w: 12.13, h: 0.9, rectRadius: 0.06,
    fill: { color: C.background2 }, line: { color: C.background2 }, objectName: "caveat box" });
  s.addText([
    { text: "CAVEATS  ", options: { bold: true, color: C.accent3, fontSize: 12, charSpacing: 2 } },
    { text: "An association in sequence, not a demonstrated mechanism (consistent with ZAP, untested).  Only 6 RSV-A and 4 RSV-B lineages.  Genes truncated at their start may be off-frame.  8.9% of sequences are exact duplicates, not yet tested for.",
      options: { fontSize: 14, color: C.text1 } },
  ], { x: 0.8, y: 5.98, w: 11.73, h: 0.9, valign: "middle", margin: 0, isTextBox: true, objectName: "caveat text" });
  s.addNotes("Suggested question to put to Yongtao: which of the three next steps first? My own ranking is the " +
    "synonymous-sites test, because it asks whether the NS1/NS2 CpG loss looks like selection on the RNA, which is " +
    "the paper's own logic (GC3 changing faster than GC12). The within-lineage check is the one a reviewer is most " +
    "likely to raise as an alternative explanation. The paramyxovirus panel is the largest piece of work and depends " +
    "on whether the sustained-versus-sporadic question is the goal.\n\n" +
    "Caveats in full: (1) the NS1/NS2 result is an observation about sequence; ZAP would need wet-lab work to " +
    "confirm. (2) The subtype contrast rests on 6 and 4 lineages. (3) A gene truncated at its start is still off-frame " +
    "within its own slice; an in-frame filter and a rebuild (about 25 minutes) would close it. (4) 8.9% of sequences " +
    "are byte-identical duplicates (largest group 12); a dedup sensitivity check has not been run.");
}

(async () => {
  fs.mkdirSync(path.dirname(OUT), { recursive: true });
  await pres.writeFile({ fileName: OUT });
  await applyTheme(OUT, THEME);
  console.log("wrote", OUT);
})();
