#!/usr/bin/env Rscript
# Script Name: figures.R
# Purpose:     All figures for the RSV composition analysis, in the form used by
#              Ye et al. 2026 (their Figs 2 and 3).
# Inputs:      data/processed/plotdata_genome.tsv   per-sequence points (local, not committed)
#              results/slopes_genome.tsv, slopes_by_gene.tsv, gc_correlation.tsv
# Outputs:     results/figures/fig1..fig4 *.png
# Notes:       Run from the repository root:  Rscript R/figures.R

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pal <- c(A = "#2C7FB8", B = "#C62824")          # RSV-A blue, RSV-B red, used throughout

# Plain research styling: no chartjunk, thin rules, readable at figure size.
theme_set(
  theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = .25, colour = "grey92"),
          panel.border = element_rect(linewidth = .4, colour = "grey35"),
          strip.background = element_rect(fill = "grey95", colour = "grey35", linewidth = .4),
          strip.text = element_text(face = "bold", size = 9.5),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9.5, colour = "grey30"),
          plot.caption = element_text(size = 8, colour = "grey40", hjust = 0),
          legend.key.size = unit(.9, "lines")))

lab <- c(gc_content = "GC content", gc3 = "GC3", gc12 = "GC12",
         per1k_CG = "CpG", per1k_GC = "GpC", per1k_GG = "GpG",
         per1k_CC = "CpC", per1k_TA = "UpA",
         oe_CG = "CpG", oe_GC = "GpC", oe_GG = "GpG", oe_CC = "CpC", oe_TA = "UpA")

## ---- Figure 1: GC content over time, one panel per lineage --------------------
## Mirrors Fig 2a of the paper: points are genomes, the line is the fitted trend.
d <- fread("data/processed/plotdata_genome.tsv")
sl <- fread("results/slopes_genome.tsv")

ann <- sl[feature == "gc_content"]
ann[, panel := paste0(lineage_grp, "  (n = ", formatC(n, big.mark = ","), ")")]
d <- merge(d, ann[, .(Subtype, lineage_grp, panel, slope, r2)],
           by = c("Subtype", "lineage_grp"))
d[, panel := factor(panel, levels = ann[order(Subtype, -n), panel])]
ann[, panel := factor(panel, levels = levels(d$panel))]
ann[, lbl := sprintf("slope %+.2f x 10^-4 / yr\nR2 = %.3f", slope * 1e4, r2)]

f1 <- ggplot(d, aes(dec_date, gc_content)) +
  geom_point(aes(colour = Subtype), size = .25, alpha = .12) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey15",
              linewidth = .6, se = TRUE, fill = "grey70") +
  geom_text(data = ann, aes(label = lbl), x = -Inf, y = Inf,
            hjust = -0.08, vjust = 1.25, size = 2.7, colour = "grey20") +
  facet_wrap(~ panel, scales = "free_y", ncol = 4) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_continuous(breaks = seq(2016, 2026, 4)) +
  labs(title = "GC content over time, fitted separately within each lineage",
       subtitle = "Each point is one genome. Fitting inside lineages keeps lineage replacement out of the trend.",
       x = "Collection date", y = "GC content (coding regions)",
       caption = "Shaded band is the 95% confidence interval of the fit.")
ggsave("results/figures/fig1_lineage_trends.png", f1, width = 10, height = 5.4, dpi = 300)

## ---- Figure 2: GC, GC3 and GC12 slopes side by side ---------------------------
## Mirrors Fig 2c. Third codon sites are mostly synonymous, so a steeper GC3 slope
## points to pressure acting on the RNA rather than on the protein.
g <- sl[feature %chin% c("gc_content", "gc3", "gc12")]
g[, site := factor(lab[feature], levels = c("GC content", "GC12", "GC3"))]
g[, lineage_grp := factor(lineage_grp, levels = ann[order(Subtype, -n), lineage_grp])]

f2 <- ggplot(g, aes(lineage_grp, slope, fill = site)) +
  geom_hline(yintercept = 0, linewidth = .4, colour = "grey20") +
  geom_col(position = position_dodge(.75), width = .68,
           colour = "grey25", linewidth = .22) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), position = position_dodge(.75),
                width = .2, linewidth = .32) +
  facet_grid(~ Subtype, scales = "free_x", space = "free_x",
             labeller = labeller(Subtype = function(x) paste0("RSV-", x))) +
  scale_fill_manual(values = c("GC content" = "grey55", "GC12" = "#80B1D3", "GC3" = "#1F5C8B"),
                    name = NULL) +
  scale_y_continuous(labels = label_number(scale = 1e4, accuracy = .1)) +
  labs(title = "Where along the codon is GC content changing?",
       subtitle = "GC3 is the third codon site, mostly synonymous; GC12 the first and second, mostly not",
       x = NULL, y = expression("Change per year " (10^-4)),
       caption = "Error bars are 95% confidence intervals.") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("results/figures/fig2_gc3_vs_gc12.png", f2, width = 9, height = 4.6, dpi = 300)

## ---- Figure 3: dinucleotide slopes per lineage --------------------------------
## Mirrors Fig 3a and 3d: frequency per 1,000 nt, and the O/E ratio.
mk <- function(feats, ylab, ttl, sub) {
  x <- sl[feature %chin% feats]
  x[, dinuc := factor(lab[feature], levels = c("CpG", "GpC", "GpG", "CpC", "UpA"))]
  x[, lineage_grp := factor(lineage_grp, levels = ann[order(Subtype, -n), lineage_grp])]
  ggplot(x, aes(dinuc, slope, colour = Subtype)) +
    geom_hline(yintercept = 0, linewidth = .4, colour = "grey20") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .18, linewidth = .4) +
    geom_point(size = 1.7) +
    facet_wrap(~ lineage_grp, ncol = 4) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(title = ttl, subtitle = sub, x = NULL, y = ylab)
}
f3a <- mk(c("per1k_CG", "per1k_GC", "per1k_GG", "per1k_CC", "per1k_TA"),
          "Change per 1,000 nt per year",
          "Dinucleotide frequency: change per year within each lineage",
          "Frequencies per 1,000 nucleotides, the unit used in the paper")
f3b <- mk(c("oe_CG", "oe_GC", "oe_GG", "oe_CC", "oe_TA"),
          "Change in O/E per year",
          "Dinucleotide observed/expected ratio: change per year",
          "O/E corrects for base composition, so this is not simply GC content moving")
f3 <- (f3a / f3b) +
  plot_annotation(caption = "Error bars are 95% confidence intervals. UpA is exploratory; it is not used in the source paper.")
ggsave("results/figures/fig3_dinucleotide_slopes.png", f3, width = 10, height = 9, dpi = 300)

## ---- Figure 4: is the dinucleotide signal just GC content moving? -------------
## Mirrors Fig 3b.
cr <- fread("results/gc_correlation.tsv")
cr[, dinuc := factor(lab[feature], levels = c("CpG", "GpC", "GpG", "CpC", "UpA"))]
cr[, lineage_grp := factor(lineage_grp, levels = ann[order(Subtype, -n), lineage_grp])]

f4 <- ggplot(cr, aes(dinuc, r, fill = Subtype)) +
  geom_hline(yintercept = 0, linewidth = .4, colour = "grey20") +
  geom_col(width = .66, colour = "grey25", linewidth = .22) +
  facet_wrap(~ lineage_grp, ncol = 4) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
  labs(title = "Correlation of each dinucleotide frequency with GC content",
       subtitle = "Values near 1 mean the dinucleotide tracks GC content; values near 0 mean it moves independently",
       x = NULL, y = "Pearson r",
       caption = "Computed within each lineage across genomes.")
ggsave("results/figures/fig4_gc_correlation.png", f4, width = 10, height = 5, dpi = 300)

cat("figures written to results/figures/\n")

## ---- Figure 5: composition gene by gene ---------------------------------------
## Context for the time trends: the genome is not uniform, and the interferon
## antagonist genes NS1 and NS2 sit far from the rest in RSV-A.
GENES <- c("NS1","NS2","N","P","M","SH","G","F","M2","L")
gs <- fread("results/gene_summary.tsv")
gs[, gene := factor(gene, levels = GENES)]

q1 <- ggplot(gs, aes(gene, gc, colour = Subtype, group = Subtype)) +
  geom_line(linewidth = .45, alpha = .75) + geom_point(size = 2.1) +
  scale_colour_manual(values = pal) +
  labs(title = "GC content varies strongly between genes", x = NULL, y = "Mean GC content")

q2 <- ggplot(gs, aes(gene, oe_CG, colour = Subtype, group = Subtype)) +
  geom_line(linewidth = .45, alpha = .75) + geom_point(size = 2.1) +
  scale_colour_manual(values = pal) +
  labs(title = "CpG suppression is not uniform along the genome",
       subtitle = "Lower values mean stronger suppression; NS1 and NS2 in RSV-A sit near zero",
       x = "Gene (genome order, 3' to 5')", y = "Mean CpG O/E")

f5 <- (q1 / q2) + plot_layout(guides = "collect")
ggsave("results/figures/fig5_gene_composition.png", f5, width = 8, height = 6.2, dpi = 300)
cat("fig5 written\n")
