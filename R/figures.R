#!/usr/bin/env Rscript
# Script Name: figures.R
# Purpose:     Publication-quality figures of RSV genome composition (GC content and
#              GC-dinucleotide O/E ratios) by subtype, lineage and year.
# Inputs:      data/processed/rsv_features.tsv  (per-sequence features; local, git-ignored)
# Outputs:     results/figures/fig1..fig3 *.png, results/yearly_composition_summary.tsv
# Notes:       Run from repo root: Rscript R/figures.R

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
d <- fread("data/processed/rsv_features.tsv")
d <- d[!is.na(year) & year >= 2015 & year <= 2026]

# shared theme + subtype colours
theme_set(theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")))
pal <- c(A = "#2c7fb8", B = "#de2d26")

oe_lab <- c(oe_CG = "CpG", oe_GC = "GpC", oe_GG = "GpG", oe_CC = "CpC")

## Fig 1 — GC content by subtype -------------------------------------------------
f1 <- ggplot(d, aes(Subtype, gc_content, fill = Subtype)) +
  geom_violin(alpha = .5, colour = NA) +
  geom_boxplot(width = .15, outlier.size = .3) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "GC content by RSV subtype", x = NULL, y = "GC content")
ggsave("results/figures/fig1_gc_by_subtype.png", f1, width = 5, height = 4, dpi = 200)

## Fig 2 — four dinucleotide O/E ratios by subtype -------------------------------
m <- melt(d, id.vars = "Subtype", measure.vars = names(oe_lab),
          variable.name = "dinuc", value.name = "oe_val")
m[, dinuc := factor(oe_lab[as.character(dinuc)], levels = unname(oe_lab))]
f2 <- ggplot(m, aes(dinuc, oe_val, fill = Subtype)) +
  geom_boxplot(outlier.size = .2, alpha = .8) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = pal) +
  labs(title = "Dinucleotide observed/expected ratios",
       subtitle = "Dashed line = 1 (as expected); values <1 indicate suppression",
       x = NULL, y = "O/E ratio")
ggsave("results/figures/fig2_oe_dinuc_by_subtype.png", f2, width = 7, height = 4.5, dpi = 200)

## Fig 3 — yearly trends of GC and CpG O/E --------------------------------------
yr <- d[, .(n = .N, gc = mean(gc_content), gc_se = sd(gc_content)/sqrt(.N),
            oeCG = mean(oe_CG), oeCG_se = sd(oe_CG)/sqrt(.N)),
        by = .(Subtype, year)][n >= 20]
fwrite(yr, "results/yearly_composition_summary.tsv", sep = "\t")

p_gc <- ggplot(yr, aes(year, gc, colour = Subtype, fill = Subtype)) +
  geom_ribbon(aes(ymin = gc - gc_se, ymax = gc + gc_se), alpha = .2, colour = NA) +
  geom_line(linewidth = .8) + geom_point(size = 1.3) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(title = "GC content over time", x = NULL, y = "Mean GC content")
p_cg <- ggplot(yr, aes(year, oeCG, colour = Subtype, fill = Subtype)) +
  geom_ribbon(aes(ymin = oeCG - oeCG_se, ymax = oeCG + oeCG_se), alpha = .2, colour = NA) +
  geom_line(linewidth = .8) + geom_point(size = 1.3) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  scale_x_continuous(breaks = seq(2015, 2026, 2)) +
  labs(title = "CpG O/E ratio over time", x = "Collection year", y = "Mean CpG O/E")
f3 <- (p_gc / p_cg) + plot_layout(guides = "collect")
ggsave("results/figures/fig3_yearly_trends.png", f3, width = 7, height = 7, dpi = 200)

cat("Figures written to results/figures/\n")
