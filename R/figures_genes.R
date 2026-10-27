#!/usr/bin/env Rscript
# Script Name: figures_genes.R
# Purpose:     Per-gene composition figures — baseline composition across the 10 RSV genes,
#              and which genes carry the within-lineage CpG O/E trend.
# Inputs:      results/gene_summary.tsv, results/gene_decomposition.tsv
# Outputs:     results/figures/fig9_gene_composition.png, fig10_gene_trends.png
# Notes:       Run after run_gene_analysis.py.  Rscript R/figures_genes.R

suppressPackageStartupMessages({
  library(ggrepel)
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pal  <- c(A = "#2c7fb8", B = "#de2d26")
pal2 <- c(total = "#b2182b", within = "#2166ac")
GENES <- c("NS1", "NS2", "N", "P", "M", "SH", "G", "F", "M2", "L")

th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(linewidth = .4, colour = "grey30"),
        strip.background = element_rect(fill = "grey95", colour = "grey30", linewidth = .4),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))
theme_set(th)

## ---- Fig 9: baseline composition across genes ---------------------------------
s <- fread("results/gene_summary.tsv")
s[, gene := factor(gene, levels = GENES)]

p1 <- ggplot(s, aes(gene, gc, colour = Subtype, group = Subtype)) +
  geom_line(linewidth = .5, alpha = .7) + geom_point(size = 2.4) +
  scale_colour_manual(values = pal) +
  labs(title = "GC content varies strongly between RSV genes",
       x = NULL, y = "Mean GC content")
p2 <- ggplot(s, aes(gene, oe_CG, colour = Subtype, group = Subtype)) +
  geom_line(linewidth = .5, alpha = .7) + geom_point(size = 2.4) +
  scale_colour_manual(values = pal) +
  labs(title = "CpG suppression is not uniform along the genome",
       subtitle = "Lower values = stronger CpG suppression",
       x = "Gene (genome order 3' to 5')", y = "Mean CpG O/E")
f9 <- (p1 / p2) + plot_layout(guides = "collect")
ggsave("results/figures/fig9_gene_composition.png", f9, width = 8, height = 6.5, dpi = 300)

## ---- Fig 10: which genes carry the within-lineage trend -----------------------
d <- fread("results/gene_decomposition.tsv")
d <- d[feature == "oe_CG"]
d[, gene := factor(gene, levels = GENES)]
m <- melt(d, id.vars = c("Subtype", "gene", "p_within_adj"),
          measure.vars = c("slope_total", "slope_within"),
          variable.name = "model", value.name = "slope")
se <- melt(d, id.vars = c("Subtype", "gene"),
           measure.vars = c("se_total", "se_within"), value.name = "se")
m[, se := se$se]
m[, model := factor(fifelse(model == "slope_total", "total", "within"),
                    levels = c("total", "within"))]
m[, star := fifelse(model == "within" & p_within_adj < 0.05, "*", "")]

f10 <- ggplot(m, aes(gene, slope, fill = model)) +
  geom_hline(yintercept = 0, linewidth = .4, colour = "grey20") +
  geom_col(position = position_dodge(.72), width = .64,
           colour = "grey20", linewidth = .22) +
  geom_errorbar(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se),
                position = position_dodge(.72), width = .18, linewidth = .32) +
  geom_text(aes(label = star, y = slope + sign(slope) * 1.5 * se),
            position = position_dodge(.72), size = 4, colour = "grey10") +
  facet_wrap(~ paste0("RSV-", Subtype), ncol = 1, scales = "free_y") +
  scale_fill_manual(values = pal2, name = NULL,
                    labels = c("Total (raw)", "Within-lineage")) +
  scale_y_continuous(labels = label_number(scale = 1e3, accuracy = .5)) +
  labs(title = "Which genes drive the CpG O/E trend?",
       subtitle = "Within-lineage slope isolates evolution from lineage replacement",
       x = "Gene (genome order 3' to 5')",
       y = expression("CpG O/E slope per year " (10^-3)),
       caption = "Error bars = 95% CI. * = within-lineage slope significant after Benjamini-Hochberg correction.") +
  theme(legend.position = "top")
ggsave("results/figures/fig10_gene_trends.png", f10, width = 8.5, height = 7, dpi = 300)

cat("gene figures written to results/figures/\n")

## ---- Fig 11: CpG depletion per gene, RSV-A vs RSV-B ---------------------------
cp <- fread("results/cpg_depletion_by_gene.tsv")
cp[, gene := factor(gene, levels = GENES)]

q1 <- ggplot(cp, aes(gene, cpg_per_gene, fill = Subtype)) +
  geom_col(position = position_dodge(.72), width = .64,
           colour = "grey20", linewidth = .22) +
  scale_fill_manual(values = pal) +
  labs(title = "CpG dinucleotides per gene copy",
       subtitle = "RSV-A has almost eliminated CpG from NS1 and NS2, the two interferon-antagonist genes",
       x = NULL, y = "Mean CpG count per gene")

q2 <- ggplot(cp, aes(gene, pct_zero_cpg, fill = Subtype)) +
  geom_col(position = position_dodge(.72), width = .64,
           colour = "grey20", linewidth = .22) +
  scale_fill_manual(values = pal) +
  labs(title = "Share of sequences carrying no CpG at all",
       x = "Gene (genome order 3' to 5')", y = "% of sequences with zero CpG")

q3 <- ggplot(cp, aes(gc, oe_CG, colour = Subtype)) +
  geom_point(aes(fill = Subtype), size = 2.6, shape = 21, colour = "grey20") +
  ggrepel::geom_text_repel(aes(label = gene), size = 2.8, show.legend = FALSE,
                           max.overlaps = 20, segment.size = .2) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal) +
  labs(title = "CpG suppression is not explained by GC content",
       subtitle = "NS1 has the same GC in both subtypes, but very different CpG O/E",
       x = "Mean GC content", y = "Mean CpG O/E")

f11 <- (q1 / q2 / q3) + plot_layout(guides = "collect", heights = c(1, 1, 1.2))
ggsave("results/figures/fig11_cpg_depletion.png", f11, width = 8.5, height = 11, dpi = 300)
cat("fig11 written\n")
