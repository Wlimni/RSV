/*
 * Script Name: make_onboarding.js
 * Purpose:     Build the RSV project onboarding document (.docx) for someone joining with
 *              limited background in virology or bioinformatics.
 * Inputs:      none (content is inline)
 * Outputs:     docs/RSV_Onboarding.docx
 * Notes:       Run with: node make_onboarding.js
 */
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  LevelFormat, PageBreak, convertInchesToTwip
} = require('docx');
const fs = require('fs');

const ACCENT = '1F4E79';     // deep blue, used for headings
const ACCENT_A = '2C7FB8';   // RSV-A blue, matches the figures
const ACCENT_B = 'C62824';   // RSV-B red
const GREY = '595959';
const BOXBG = 'EEF3F8';
const WARNBG = 'FDF4E3';

const PAGE_W = 12240, PAGE_H = 15840;          // US Letter in DXA
const MARGIN = convertInchesToTwip(1);
const CONTENT_W = PAGE_W - 2 * MARGIN;         // usable width for tables

// ---------- small builders ----------
const p = (text, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 120, line: 276 },
  alignment: opts.align,
  children: [new TextRun({
    text, size: opts.size ?? 22, bold: opts.bold, italics: opts.italics,
    color: opts.color, font: 'Calibri'
  })]
});

// paragraph from mixed runs: ['plain', ['bold bit', {bold:true}], ...]
const rich = (parts, opts = {}) => new Paragraph({
  spacing: { after: opts.after ?? 120, line: 276 },
  children: parts.map(x => {
    const [t, o] = Array.isArray(x) ? x : [x, {}];
    return new TextRun({ text: t, size: 22, font: 'Calibri', ...o });
  })
});

const h1 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_1,
  spacing: { before: 360, after: 160 },
  children: [new TextRun({ text, size: 32, bold: true, color: ACCENT, font: 'Calibri' })]
});

const h2 = (text) => new Paragraph({
  heading: HeadingLevel.HEADING_2,
  spacing: { before: 260, after: 120 },
  children: [new TextRun({ text, size: 26, bold: true, color: ACCENT, font: 'Calibri' })]
});

const bullet = (text, level = 0) => new Paragraph({
  numbering: { reference: 'bullets', level },
  spacing: { after: 80, line: 276 },
  children: [new TextRun({ text, size: 22, font: 'Calibri' })]
});

const numbered = (text) => new Paragraph({
  numbering: { reference: 'numbers', level: 0 },
  spacing: { after: 80, line: 276 },
  children: [new TextRun({ text, size: 22, font: 'Calibri' })]
});

// shaded single-cell callout
const callout = (title, lines, bg = BOXBG) => new Table({
  columnWidths: [CONTENT_W],
  width: { size: CONTENT_W, type: WidthType.DXA },
  rows: [new TableRow({
    children: [new TableCell({
      width: { size: CONTENT_W, type: WidthType.DXA },
      shading: { type: ShadingType.CLEAR, fill: bg },
      margins: { top: 160, bottom: 160, left: 180, right: 180 },
      children: [
        new Paragraph({
          spacing: { after: 90 },
          children: [new TextRun({ text: title, bold: true, size: 22, color: ACCENT, font: 'Calibri' })]
        }),
        ...lines.map((l, i) => new Paragraph({
          spacing: { after: i === lines.length - 1 ? 0 : 90, line: 276 },
          children: [new TextRun({ text: l, size: 21, font: 'Calibri' })]
        }))
      ]
    })]
  })]
});

// data table with header row
const table = (headers, rows, widths) => {
  const cw = widths.map(w => Math.round(CONTENT_W * w));
  const diff = CONTENT_W - cw.reduce((a, b) => a + b, 0);
  cw[cw.length - 1] += diff;                      // make widths sum exactly
  // numeric-looking columns (anything after the first) are right-aligned unless the
  // table is mostly prose, which we detect by the column being wide.
  const rightAlign = (i) => i > 0 && widths[i] <= 0.34;
  const mkRow = (cells, opts) => new TableRow({
    tableHeader: opts?.header,
    children: cells.map((t, i) => new TableCell({
      width: { size: cw[i], type: WidthType.DXA },
      shading: opts?.header
        ? { type: ShadingType.CLEAR, fill: 'DCE6F1' }
        : (opts?.bg ? { type: ShadingType.CLEAR, fill: opts.bg } : undefined),
      margins: { top: 90, bottom: 90, left: 120, right: 120 },
      children: [new Paragraph({
        alignment: rightAlign(i) ? AlignmentType.RIGHT : undefined,
        spacing: { after: 0, line: 252 },
        children: [new TextRun({ text: t, size: 20, bold: opts?.header, font: 'Calibri' })]
      })]
    }))
  });
  return new Table({
    columnWidths: cw,
    width: { size: CONTENT_W, type: WidthType.DXA },
    rows: [mkRow(headers, { header: true }), ...rows.map(r => mkRow(r.cells ?? r, { bg: r.bg }))]
  });
};

const spacer = (h = 120) => new Paragraph({ spacing: { after: h }, children: [] });
const rule = () => new Paragraph({
  spacing: { before: 60, after: 160 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: 'BFBFBF' } },
  children: []
});

// ---------- document content ----------
const children = [];

// Title block
children.push(new Paragraph({
  spacing: { before: 1600, after: 80 },
  children: [new TextRun({ text: 'Onboarding Guide', size: 52, bold: true, color: ACCENT, font: 'Calibri' })]
}));
children.push(new Paragraph({
  spacing: { after: 300 },
  children: [new TextRun({
    text: 'Genomic composition and host transmission in RSV',
    size: 30, color: GREY, font: 'Calibri'
  })]
}));
children.push(rule());
children.push(p('Written for someone joining this project with limited background in virology or bioinformatics. No prior knowledge of RSV, genome sequences or statistics is assumed. Every technical term is explained the first time it appears, and there is a glossary at the end.', { size: 22, color: GREY }));
children.push(spacer(200));
children.push(rich([['Project owner: ', { bold: true }], 'William Ni Hong Hung, The University of Hong Kong']));
children.push(rich([['Code repository: ', { bold: true }], 'github.com/Wlimni/RSV']));
children.push(rich([['Document date: ', { bold: true }], '28 October 2026']));

children.push(new Paragraph({ children: [new PageBreak()] }));

// Contents
children.push(h1('Contents'));
[
  '1.  The project in one page',
  '2.  Background: what you need to know first',
  '3.  The paper this project is based on',
  '4.  The purpose of this project',
  '5.  What we have built',
  '6.  Results so far',
  '7.  The limitation that shapes everything next',
  '8.  What is probably next',
  '9.  How to get started with the code',
  '10. Glossary'
].forEach(t => children.push(p(t, { after: 90 })));

children.push(new Paragraph({ children: [new PageBreak()] }));

// 1. One page
children.push(h1('1.  The project in one page'));
children.push(p('Viruses are made of genetic material that can be written down as a long string of four letters: A, C, G and T. This project asks a simple question: if you look only at the letter patterns in a virus genome, without knowing anything else about it, can you tell something useful about how that virus spreads?'));
children.push(p('The idea comes from a 2026 study of influenza. That study found that flu strains which managed to establish themselves long-term in mammals had a measurably different letter composition from those that infected a mammal once and then died out. If genome composition really carries that signal, it could help flag a dangerous new virus early.'));
children.push(p('Our task was to test whether the same kind of signal exists in a different virus family, starting with RSV, a common respiratory virus. We built a complete analysis pipeline, applied it to 33,038 RSV genomes, and found one clear and interesting pattern, plus one important limitation that determines what the project should do next.'));
children.push(spacer());
children.push(callout('The three things to remember', [
  '1. We measure letter composition of virus genomes, especially how often a C is followed by a G (called "CpG").',
  '2. Our clearest finding: RSV subtype A has almost completely removed CpG from the two genes that block the human immune alarm system. Subtype B has not.',
  '3. Our main limitation: the original question was about viruses that spread sustainably versus ones that fizzle out. All RSV data is from humans where it already spreads sustainably, so RSV alone cannot answer that question.'
]));

// 2. Background
children.push(h1('2.  Background: what you need to know first'));

children.push(h2('2.1  What RSV is'));
children.push(p('RSV stands for Respiratory Syncytial Virus. It is a very common respiratory infection. Almost every child catches it by the age of two. In most people it causes something like a bad cold, but in infants, the elderly and people with weak immune systems it can cause serious lung infection and hospitalisation. It is one of the leading causes of infant hospitalisation worldwide.'));
children.push(p('RSV comes in two forms, called subtype A and subtype B. They are closely related but genetically distinct, and they circulate at the same time, often taking turns dominating from season to season. Throughout this project we analyse A and B separately, because they turn out to behave differently.'));
children.push(p('Since 2023 there are new ways to prevent severe RSV: a monoclonal antibody called nirsevimab given to infants, and a vaccine for older adults. Both target a virus protein called F, which becomes relevant later.'));

children.push(h2('2.2  What a virus genome is, and how we read it'));
children.push(p('A genome is the full set of genetic instructions for an organism. For RSV it is a single strand of RNA about 15,200 units long. Each unit is one of four chemical letters, written A, C, G and T. So to a computer, an RSV genome is simply a string of about 15,200 characters.'));
children.push(p('That string is not random. Parts of it are genes, each encoding one protein that does a specific job. RSV has ten genes. In the order they appear along the genome, they are:'));
children.push(p('NS1  -  NS2  -  N  -  P  -  M  -  SH  -  G  -  F  -  M2  -  L', { bold: true, align: AlignmentType.CENTER }));
children.push(p('Four of these matter for this document:'));
children.push(bullet('NS1 and NS2 are "interferon antagonists". Interferon is the alarm signal a human cell sends when it detects a virus. These two genes produce proteins that switch that alarm off, so the virus can replicate before the body reacts.'));
children.push(bullet('G is the attachment protein. It sits on the outside of the virus and grabs onto human cells. Because it is exposed, antibodies attack it constantly, so it is the fastest-changing gene.'));
children.push(bullet('F is the fusion protein, also on the outside. It is the target of nirsevimab and the new vaccine.'));
children.push(bullet('L is by far the largest gene, about 6,500 letters, roughly 43% of the whole genome. Because of its size it dominates any measurement taken across the whole genome at once.'));

children.push(h2('2.3  What we mean by "composition"'));
children.push(p('Composition means the statistics of the letters, ignoring what they code for. The two kinds we measure:'));
children.push(spacer(80));
children.push(rich([['GC content. ', { bold: true }], 'The fraction of all letters that are either G or C. For RSV this is about 0.33, meaning roughly a third of the genome is G or C and two thirds is A or T.']));
children.push(rich([['Dinucleotide frequency. ', { bold: true }], 'How often a specific pair of adjacent letters occurs. We track four pairs that involve only G and C:']));
children.push(spacer(80));
children.push(table(
  ['Pair', 'Written as', 'Means'],
  [
    ['C then G', 'CpG', 'a C immediately followed by a G'],
    ['G then C', 'GpC', 'a G immediately followed by a C'],
    ['G then G', 'GpG', 'two Gs in a row'],
    ['C then C', 'CpC', 'two Cs in a row']
  ],
  [0.22, 0.22, 0.56]
));
children.push(spacer(140));
children.push(p('The "p" in CpG just stands for the chemical bond between the two letters. It is there to make clear we mean a C next to a G along one strand, not a C paired with a G across the two strands of a double helix.'));

children.push(h2('2.4  The observed / expected ratio, and why we need it'));
children.push(p('Counting pairs directly is misleading. A genome with lots of Gs and Cs will naturally contain more CpG pairs just by chance. To ask whether a virus is actively avoiding a pair, we have to compare against what chance alone would produce.'));
children.push(p('That comparison is the observed / expected ratio, written O/E:'));
children.push(spacer(80));
children.push(p('O/E  =  frequency of CpG  ÷  ( frequency of C  ×  frequency of G )', { bold: true, align: AlignmentType.CENTER }));
children.push(spacer(80));
children.push(p('A worked example. Suppose 15% of letters are C and 15% are G. If letters were arranged randomly, the chance that any given position starts a CpG pair is 0.15 × 0.15 = 0.0225, so we would expect about 2.25% of positions to be CpG. If we actually count and find only 0.5%, then:'));
children.push(p('O/E  =  0.005 ÷ 0.0225  =  0.22', { bold: true, align: AlignmentType.CENTER }));
children.push(spacer(80));
children.push(p('Interpretation is straightforward:'));
children.push(bullet('O/E near 1.0 means the pair appears about as often as chance predicts.'));
children.push(bullet('O/E below 1.0 means the pair is suppressed, appearing less often than chance. The virus is avoiding it.'));
children.push(bullet('O/E above 1.0 means the pair is enriched.'));

children.push(h2('2.5  Why CpG in particular matters'));
children.push(p('Human cells carry an antiviral protein called ZAP (zinc-finger antiviral protein). ZAP physically binds to viral RNA at CpG sites and marks that RNA for destruction. In effect, each CpG in a virus genome is a handle the immune system can grab.'));
children.push(p('The consequence is that viruses which circulate in humans over long periods tend to evolve fewer CpG sites than chance would give them. They file off the handles. This is why CpG O/E, rather than GC content alone, is the measurement people pay attention to.'));
children.push(spacer());
children.push(callout('Important caution', [
  'The ZAP explanation is a well-supported idea in the field, but in this project it remains an interpretation of patterns in sequence data. We have not demonstrated the mechanism ourselves. Doing so would require laboratory experiments, not sequence analysis.'
], WARNBG));

children.push(h2('2.6  Where the data comes from'));
children.push(p('GISAID is an international database where laboratories deposit virus genome sequences along with information about each sample: the date it was collected, the country, the host, and which branch of the virus family it belongs to. It was built for influenza and expanded to other viruses including RSV.'));
children.push(p('GISAID data comes with terms of use that restrict redistribution. For that reason the project repository contains no GISAID sequences or metadata at all. It holds only our code, summary tables and figures. Anyone reproducing the work downloads the data themselves.'));

children.push(h2('2.7  Lineages'));
children.push(p('As a virus spreads it accumulates small changes, and the population splits into branches of a family tree. Each named branch is a lineage, for example "A.D.1" or "B.D.E.1". GISAID assigns a lineage to every sequence.'));
children.push(p('Lineages rise and fall. A new lineage may spread rapidly, dominate for a few seasons, then be replaced by another. This turnover is central to understanding our results, because it can create apparent trends in the data that are not real changes in the virus. Section 6 explains this.'));

children.push(new Paragraph({ children: [new PageBreak()] }));

// 3. The paper
children.push(h1('3.  The paper this project is based on'));
children.push(rich([
  'Ye et al., 2026, ', ['Nature Microbiology', { italics: true }],
  ': "Genomic features associated with sustained mammalian transmission of avian influenza A viruses" (doi:10.1038/s41564-025-02257-4). A copy is in the project folder.'
]));

children.push(h2('3.1  The problem the paper addresses'));
children.push(p('Influenza A viruses live naturally in wild birds. Occasionally one jumps into a mammal: a pig, a horse, a mink, a human. When that happens, one of two things follows.'));
children.push(spacer(80));
children.push(table(
  ['Outcome', 'What happens', 'How common'],
  [
    ['Sporadic', 'A handful of infections, then it dies out. Weeks to months. Does not establish a chain of transmission.', 'The usual outcome'],
    ['Sustained', 'It establishes in the new host and keeps circulating for years or decades.', 'Rare']
  ],
  [0.17, 0.58, 0.25]
));
children.push(spacer(140));
children.push(p('The distinction is the whole point. A sporadic spillover is a localised scare. A sustained one is how pandemics begin. Human seasonal flu, swine flu, horse flu and dog flu are all sustained lineages that originally came from birds.'));
children.push(p('So the question is: when a new bird flu infects a mammal, can you tell from its genome which of these two futures it is heading for?'));

children.push(h2('3.2  What they did and found'));
children.push(p('They measured GC content and the four GC dinucleotide frequencies across 115,520 influenza genomes, then grouped the viruses into three classes based on their known history: birds only, sporadic mammalian infections, and sustained mammalian lineages.'));
children.push(p('They found that sustained mammalian lineages had lower GC content and fewer GC dinucleotides, and that these lineages could be reliably separated from the others. Critically, the earliest viruses of each sustained lineage already showed reduced GC content, before the lineage had established itself. That is what makes the signal potentially predictive rather than merely a consequence of having adapted.'));
children.push(p('Their conclusion was carefully worded: reduced GC-related content is probably necessary but not sufficient for sustained mammalian transmission, and it should be one input among several in pandemic risk assessment tools.'));
children.push(spacer());
children.push(callout('The point that matters most for us', [
  'The labels "sustained" and "sporadic" do not come from the genome. They come from knowing each virus\'s history: which host it infected and for how long. The genome is the input; the label is what you are trying to predict.',
  'This means any attempt to repeat the analysis needs a dataset containing examples of both classes. Section 7 explains why that becomes a problem for RSV.'
]));

// 4. Purpose
children.push(h1('4.  The purpose of this project'));
children.push(p('The brief, as given by the supervising team, was: study genomic features associated with host transmission for paramyxoviruses, and try RSV first.'));
children.push(p('Paramyxoviruses are a family of viruses that includes RSV, measles, mumps, and the bat-borne Nipah and Hendra viruses. The underlying hope is that the influenza result generalises, and that a similar composition signal exists in this family too. If it does, the method could eventually extend to other viruses, such as mpox.'));
children.push(p('RSV was chosen as the starting point because it has a very large amount of public sequence data, which makes it a good virus to build and test a pipeline on. The phrase "try RSV first" is important: RSV was always intended as a starting point, not as the final answer.'));
children.push(spacer());
children.push(p('Concretely, the project set out to:'));
children.push(numbered('Collect and clean a large set of RSV genomes.'));
children.push(numbered('Measure GC content and GC dinucleotide composition for each genome, following the paper\'s method.'));
children.push(numbered('Look for patterns across subtype, gene, geography, lineage and time.'));
children.push(numbered('Judge whether the approach is worth extending to the wider paramyxovirus family.'));

children.push(new Paragraph({ children: [new PageBreak()] }));

// 5. What we built
children.push(h1('5.  What we have built'));
children.push(p('The pipeline runs in five stages. Each stage is a script, and the notebooks walk through what each one produces.'));

children.push(h2('5.1  Data collection and merging'));
children.push(p('RSV sequences were downloaded from GISAID in eight batches, split by subtype and date range, because the database limits how much can be exported at once. Each batch arrives as one file of sequences plus up to three separate tables of metadata, which have to be matched up by accession number.'));
children.push(p('Merging produced 33,040 genomes: 15,458 of subtype A and 17,582 of subtype B, collected between 2001 and 2026. Six sequences appeared in more than one batch and were de-duplicated.'));

children.push(h2('5.2  Quality control'));
children.push(p('Not every deposited sequence is usable. Some are incomplete; some contain stretches of "N", meaning the sequencing machine could not determine the letter. Since our entire analysis is counting letters, ambiguous letters would distort the result.'));
children.push(p('We keep a sequence only if it passes all three checks:'));
children.push(bullet('Length between 14,900 and 15,400 letters. A complete RSV genome is about 15,200, so this window keeps full genomes and excludes fragments.'));
children.push(bullet('At most 1% ambiguous letters.'));
children.push(bullet('A complete collection date, so it can be placed in time.'));
children.push(p('33,038 of the 33,040 sequences passed. The filter is deliberately gentle because the data had already been restricted to complete, high-coverage genomes at download time.'));

children.push(h2('5.3  Measuring composition'));
children.push(p('For every genome we compute GC content, the frequency of the four dinucleotides, and their four O/E ratios. Positions containing ambiguous letters are skipped rather than guessed.'));

children.push(h2('5.4  Splitting each genome into its ten genes'));
children.push(p('Measuring the whole genome at once hides a lot, because the genes do different jobs and are under different pressures. The paper did the equivalent by analysing each influenza segment separately.'));
children.push(p('To do this we align each genome against a reference genome whose gene boundaries are already annotated (obtained from the public NCBI database, which has no redistribution restrictions), then cut out each gene by position. This produced 330,380 records: 33,038 genomes times 10 genes.'));
children.push(p('One technical detail is worth recording, because getting it wrong would have corrupted the most interesting gene. Some RSV lineages carry a duplication in the G gene, an extra 72 letters inserted. A naive alignment would discard insertions and silently truncate those sequences. Our method keeps them. We verified this: the G gene comes out at 969 letters against a 897-letter reference, which is exactly the reference plus the known 72-letter duplication.'));

children.push(h2('5.5  Analysis and figures'));
children.push(p('The statistical work is deliberately kept simple: we describe patterns rather than build models. The one piece of real statistics is the separation of genuine change from lineage replacement, explained in section 6.4.'));
children.push(p('Figures are produced in R using ggplot2, which gives publication-quality output, and are displayed inside the notebooks so the analysis and the figures cannot drift out of step.'));
children.push(spacer());
children.push(callout('Tools used', [
  'Python for data handling, composition measurement and statistics (pandas, statsmodels).',
  'MAFFT for aligning genomes to the reference.',
  'R with ggplot2 for figures.',
  'Jupyter notebooks as the main documented pipeline.',
  'Git and GitHub for version control, with all GISAID data excluded.'
]));

children.push(new Paragraph({ children: [new PageBreak()] }));

// 6. Results
children.push(h1('6.  Results so far'));
children.push(p('Four findings, in increasing order of how interesting they are.'));

children.push(h2('6.1  RSV suppresses CpG strongly, and specifically'));
children.push(p('Across both subtypes, CpG O/E is about 0.24. RSV carries only about a quarter of the CpG pairs that its own base composition would predict by chance. The other three dinucleotides sit near or above 1.0.'));
children.push(spacer(80));
children.push(table(
  ['Measure', 'RSV-A', 'RSV-B'],
  [
    ['GC content', '0.3344', '0.3364'],
    [{ text: 'CpG O/E' }, '0.236', '0.241'],
    ['GpC O/E', '0.996', '1.086'],
    ['GpG O/E', '1.130', '1.054'],
    ['CpC O/E', '1.064', '1.008']
  ].map(r => Array.isArray(r) ? r.map(x => typeof x === 'object' ? x.text : x) : r),
  [0.4, 0.3, 0.3]
));
children.push(spacer(140));
children.push(p('This matters because it shows the effect is specific to CpG, not a general shortage of G and C letters. If the virus were simply low in G and C, all four ratios would be depressed. Only CpG is.'));

children.push(h2('6.2  Composition varies a lot between genes'));
children.push(p('Once we look gene by gene, the genome is clearly not uniform. The G gene stands out: it has by far the highest GC content of any gene, about 0.43 in subtype A against roughly 0.33 genome-wide, and the highest GpG ratio.'));
children.push(p('This makes G the natural region to focus on if the project wants to examine one region in detail. It is also the most biologically active region, being under constant antibody pressure.'));

children.push(h2('6.3  The clearest finding: RSV-A has eliminated CpG from its immune-evasion genes'));
children.push(p('This is the most striking result in the project.'));
children.push(spacer(80));
children.push(table(
  ['Gene', 'CpG per gene, A', 'CpG per gene, B', '% of A genomes with zero CpG'],
  [
    ['NS1', '0.16', '1.94', '84.8%'],
    ['NS2', '0.76', '5.64', '34.2%'],
    ['N', '6.86', '6.09', '0%'],
    ['L', '41.2', '38.2', '0%']
  ],
  [0.2, 0.26, 0.26, 0.28]
));
children.push(spacer(140));
children.push(p('Read the NS1 row. In subtype A the entire gene contains on average 0.16 CpG pairs, and in roughly 85% of subtype A genomes there is not a single one anywhere in the gene. In subtype B the same gene averages about 1.94 and almost always has some. NS2 shows the same pattern more weakly.'));
children.push(p('The genes N and L are shown for comparison. They behave almost identically in both subtypes, so this is not a general difference between A and B.'));
children.push(p('Crucially, GC content for NS1 is 0.305 in subtype A and 0.303 in subtype B, essentially identical. So the difference cannot be explained by one subtype simply having fewer Gs and Cs in that gene. It is specific avoidance of the CpG pair.'));
children.push(spacer());
children.push(callout('Why this is interesting', [
  'NS1 and NS2 are the two genes whose job is to shut down the cell\'s interferon alarm. CpG is the motif the antiviral protein ZAP uses to find and destroy viral RNA.',
  'So the genes that disable the host\'s defence are precisely the genes that have been most thoroughly stripped of the marker the host uses to detect them. That is a coherent story, and it is the kind of observation worth following up.',
  'It remains an observation. Confirming the mechanism would require restoring CpG sites in NS1 experimentally and measuring whether the virus replicates worse in cells that have ZAP than in cells that do not.'
]));

children.push(h2('6.4  A trend that looked real, and was not'));
children.push(p('This section is included because it is the most useful methodological lesson from the project.'));
children.push(p('Plotting average CpG O/E per year for subtype B shows a clear upward trend over about a decade. The obvious reading is that RSV-B is evolving towards more CpG, which would be surprising and would contradict the direction the influenza paper found.'));
children.push(p('But an average can move for two entirely different reasons:'));
children.push(numbered('The genomes inside each lineage are genuinely changing.'));
children.push(numbered('No lineage changes at all, but a lineage that already had high CpG grows while a lineage with low CpG shrinks. The average rises with no evolution whatsoever.'));
children.push(p('To separate these, we fit the trend twice: once plainly, and once while holding lineage constant. The second fit cannot be affected by replacement, so any trend that survives it is real change.'));
children.push(spacer(80));
children.push(table(
  ['RSV-B, CpG O/E', 'Slope per year', 'Direction'],
  [
    ['Plain trend', '+0.00065', 'rising'],
    ['Holding lineage constant', '-0.00052', 'falling']
  ],
  [0.46, 0.27, 0.27]
));
children.push(spacer(140));
children.push(p('The two have opposite signs. The apparent rise was entirely lineage replacement. Inside lineages, CpG O/E is falling, which is the same direction the influenza paper reported. This statistical trap is known as Simpson\'s paradox.'));
children.push(p('We then checked whether the within-lineage result survives the known biases in the data. Sequencing is heavily weighted towards Europe, and 2020 to 2021 was disrupted by COVID-19, which changed both RSV circulation and surveillance.'));
children.push(spacer(80));
children.push(table(
  ['Check', 'RSV-A result', 'RSV-B result'],
  [
    ['All data', 'declining', 'declining'],
    ['Excluding 2020-2021', 'disappears (p = 0.97)', 'holds'],
    ['Adjusting for region', 'holds', 'holds'],
    ['Balancing regions', 'disappears (p = 0.16)', 'holds']
  ],
  [0.34, 0.33, 0.33]
));
children.push(spacer(140));
children.push(p('So we trust the RSV-B result and treat the RSV-A one as an artefact of uneven sampling. Reporting it that way, rather than claiming both, is the honest reading.'));

children.push(new Paragraph({ children: [new PageBreak()] }));

// 7. Limitation
children.push(h1('7.  The limitation that shapes everything next'));
children.push(p('This section is the most important one in the document.'));
children.push(p('The paper\'s central comparison is between viruses that achieved sustained transmission in a host and viruses that only caused sporadic infections. As established in section 3, those labels come from each virus\'s known history, not from its sequence.'));
children.push(spacer());
children.push(callout('The problem, stated plainly', [
  'Every RSV sequence in GISAID comes from a human host, and RSV has circulated in humans continuously for decades.',
  'So our dataset contains examples of one class only: sustained, in humans. There are no sporadic RSV examples to compare against.',
  'No model, and no amount of additional RSV sequencing, can learn to tell sustained from sporadic here. The second class simply does not exist in the data.'
], WARNBG));
children.push(spacer());
children.push(p('This was recognised partway through the project. The work continued using a substitute measure: within RSV, some lineages grow and dominate while others fade away, and that competition can be treated as a loose stand-in for transmission success. That substitute produced the finding in section 6.4.'));
children.push(p('But it should be stated clearly that competition between lineages of a virus already established in humans is a different biological question from whether a virus can establish itself in a new host at all. The substitute does not answer the original brief.'));
children.push(p('None of this makes the work wasted. The pipeline is virus-agnostic and will run on any virus family. The finding in section 6.3 stands entirely on its own and does not depend on the transmission question. But answering the question that was actually asked requires a different dataset.'));

// 8. Next
children.push(h1('8.  What is probably next'));
children.push(p('A decision from the supervising team is pending. The options below are what has been proposed.'));

children.push(h2('8.1  Build a paramyxovirus panel with both classes in it'));
children.push(p('The sustained-versus-sporadic contrast does exist within the wider paramyxovirus family, and the sequences are in NCBI, which is public and has no redistribution restrictions. Nipah and Hendra viruses are the textbook sporadic case: they spill over from bats into humans and horses, cause severe disease, and almost always fail to establish onward transmission chains.'));
children.push(spacer(80));
children.push(table(
  ['Virus', 'Host situation', 'Genomes available'],
  [
    ['Nipah', 'bat to human, sporadic', '138'],
    ['Hendra', 'bat to horse/human, sporadic', '92'],
    ['Human RSV', 'sustained, human', '26,399'],
    ['Measles', 'sustained, human', '3,527'],
    ['Human metapneumovirus', 'sustained, human', '2,349'],
    ['Newcastle disease virus', 'sustained, bird', '1,912'],
    ['Human parainfluenza 3', 'sustained, human', '1,547'],
    ['Mumps', 'sustained, human', '1,380'],
    ['Canine distemper', 'sustained, dog', '348'],
    ['Avian metapneumovirus', 'sustained, bird', '215'],
    ['Bovine RSV', 'sustained, cattle', '71']
  ],
  [0.42, 0.4, 0.18]
));

children.push(h2('8.2  The trap to avoid in that design'));
children.push(p('If all of these were simply pooled and fed to a classifier, it would achieve near-perfect accuracy and mean nothing. It would be learning to tell Nipah from RSV, which is trivial, rather than learning what distinguishes sustained from sporadic transmission. Different viruses differ in a thousand ways at once.'));
children.push(p('The influenza paper avoided this because everything it compared was influenza A: one virus, many host lineages. Across different paramyxovirus genera, family relatedness and transmission class are tangled together.'));
children.push(p('The design that controls for it is matched comparisons within the same genus, where relatedness is roughly constant and host is what differs:'));
children.push(bullet('Human RSV versus bovine RSV'));
children.push(bullet('Human metapneumovirus versus avian metapneumovirus'));
children.push(bullet('Measles (human) versus canine distemper (dog)'));
children.push(bullet('Human parainfluenza versus Sendai virus (mouse)'));
children.push(bullet('Nipah and Hendra as the sporadic class'));
children.push(p('Sample sizes are very unequal, tens of sporadic genomes against tens of thousands of sustained ones, so the analysis should compare summaries at the level of virus or lineage rather than classify individual genomes, and the large groups need downsampling.'));

children.push(h2('8.3  Smaller options'));
children.push(bullet('Examine the G gene on its own in more detail, since section 6.2 identified it as the most distinctive region.'));
children.push(bullet('Clarify whether "focus on a specific region" from the supervising team meant a region of the genome or a geographic region. Only the genomic reading has been built so far.'));
children.push(bullet('Follow up the NS1 finding by testing whether the protein sequence itself permits CpG at the positions where it is missing. If it does, the absence must be selection acting on the RNA.'));
children.push(bullet('Build phylogenetic trees to place lineages in evolutionary context. This would need IQ-TREE, which is not currently installed.'));

children.push(new Paragraph({ children: [new PageBreak()] }));

// 9. Getting started
children.push(h1('9.  How to get started with the code'));
children.push(p('The repository is at github.com/Wlimni/RSV. The notebooks are the main pipeline and are written to be read in order; start there rather than with the scripts.'));
children.push(spacer(80));
children.push(table(
  ['Location', 'What it contains'],
  [
    ['notebooks/', 'The documented analysis. Read 01 to 05 in order. Figures appear inline.'],
    ['src/rsv/', 'Core functions the notebooks import: composition measurement, data loading and QC, gene extraction, trend statistics.'],
    ['scripts/', 'Command-line stages: merging batches, building features, per-gene extraction, analyses.'],
    ['R/', 'Figure scripts. The notebooks call these, so you do not normally run them yourself.'],
    ['results/', 'Summary tables and all figures. Safe to share; contains no GISAID records.'],
    ['docs/', 'PROJECT_BRIEF.md, the short write-up prepared for the supervising team.'],
    ['data/', 'GISAID sequences and derived tables. Excluded from the repository; must be downloaded locally.']
  ],
  [0.22, 0.78]
));
children.push(spacer(140));
children.push(p('A practical note on the environment. The machine\'s Anaconda installation has a version conflict that breaks some scientific packages, so the project uses its own isolated environment in a folder called .venv. Use that rather than the system Python. The Jupyter kernel is registered under the name "rsv".'));
children.push(spacer());
children.push(callout('Suggested first hour', [
  '1. Read sections 2 and 3 of this document again, so the vocabulary is familiar.',
  '2. Open notebook 01 and look at how much data there is and where it comes from.',
  '3. Open notebook 03 and find the NS1 comparison table. That is the project\'s main result.',
  '4. Open notebook 04 and follow the argument about why the raw time trend is misleading.',
  '5. Read docs/PROJECT_BRIEF.md for the condensed version and the open questions.'
]));

// 10. Glossary
children.push(h1('10.  Glossary'));
const glossary = [
  ['Accession number', 'The unique identifier GISAID or NCBI gives to each deposited sequence. GISAID ones take the form EPI_ISL_ followed by digits.'],
  ['Alignment', 'Lining two or more sequences up against each other so that corresponding positions match, which lets you find the same gene in different genomes.'],
  ['Ambiguous base', 'A position where sequencing could not determine the letter, usually written N.'],
  ['CpG', 'A C immediately followed by a G along the sequence. The motif the antiviral protein ZAP recognises.'],
  ['Dinucleotide', 'Any pair of adjacent letters in a sequence.'],
  ['GC content', 'The fraction of letters in a sequence that are G or C.'],
  ['Gene', 'A stretch of the genome encoding one protein. RSV has ten.'],
  ['GISAID', 'The international database of virus genome sequences used as our data source. Its terms restrict redistribution.'],
  ['Interferon', 'The alarm signal a human cell releases when it detects a virus. RSV\'s NS1 and NS2 proteins suppress it.'],
  ['Lineage', 'A named branch of a virus family tree. Lineages rise and fall over time.'],
  ['MAFFT', 'The software used to align sequences to the reference genome.'],
  ['NCBI', 'The public US sequence database. Unlike GISAID it has no redistribution restrictions, which is why our reference genomes come from there.'],
  ['Nirsevimab', 'A monoclonal antibody given to infants to prevent severe RSV, in use since 2023. It targets the F protein.'],
  ['O/E ratio', 'Observed divided by expected. How often a dinucleotide actually appears, divided by how often chance would produce it. Below 1 means suppressed.'],
  ['Paramyxovirus', 'The virus family containing RSV, measles, mumps, Nipah and Hendra.'],
  ['QC (quality control)', 'The filtering step that removes incomplete or unreliable sequences before analysis.'],
  ['Reference genome', 'A single well-characterised genome with its genes already annotated, used as a map for locating genes in other genomes.'],
  ['RSV', 'Respiratory Syncytial Virus. Comes in subtypes A and B.'],
  ['Simpson\'s paradox', 'When a trend in combined data disappears or reverses once the data is split into its real groups. Section 6.4 is an example.'],
  ['Sporadic transmission', 'A virus infects a new host species a few times then dies out without establishing a transmission chain.'],
  ['Sustained transmission', 'A virus establishes in a host species and keeps circulating for years.'],
  ['Subtype', 'A major division within a virus species. RSV has subtypes A and B.'],
  ['ZAP', 'Zinc-finger antiviral protein. A human protein that binds viral RNA at CpG sites and targets it for destruction.']
];
children.push(table(['Term', 'Meaning'], glossary, [0.24, 0.76]));

// ---------- assemble ----------
const doc = new Document({
  creator: 'RSV project',
  title: 'RSV Project Onboarding Guide',
  description: 'Onboarding guide for the RSV genomic composition project',
  numbering: {
    config: [
      {
        reference: 'bullets',
        levels: [
          { level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 460, hanging: 240 } } } },
          { level: 1, format: LevelFormat.BULLET, text: '◦', alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 900, hanging: 240 } } } }
        ]
      },
      {
        reference: 'numbers',
        levels: [
          { level: 0, format: LevelFormat.DECIMAL, text: '%1.', alignment: AlignmentType.LEFT,
            style: { paragraph: { indent: { left: 460, hanging: 240 } } } }
        ]
      }
    ]
  },
  sections: [{
    properties: {
      page: {
        size: { width: PAGE_W, height: PAGE_H },
        margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN }
      }
    },
    children
  }]
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(process.argv[2] || 'RSV_Onboarding.docx', buf);
  console.log('wrote', process.argv[2] || 'RSV_Onboarding.docx');
});
