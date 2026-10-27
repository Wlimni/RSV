#!/usr/bin/env Rscript
# Script Name: figures_trends.R
# Purpose:     Literature-style figures for the within/between-lineage decomposition,
#              per-lineage trends, lineage growth vs composition, and robustness checks.
# Inputs:      results/decomposition.tsv, per_lineage_slopes.tsv, growth_vs_composition.tsv,
#              robustness.tsv, and data/processed/rsv_features.tsv (local)
# Outputs:     results/figures/fig5..fig8 *.png
# Notes:       Run from repo root: Rscript R/figures_trends.R

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(ggrepel); library(scales)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

pal  <- c(A = "#2c7fb8", B = "#de2d26")
pal2 <- c(total = "#b2182b", within = "#2166ac")

# Journal-style base theme: no chartjunk, clear hierarchy, readable at column width.
th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(linewidth = .3, colour = "grey90"),
        panel.border = element_rect(linewidth = .4, colour = "grey30"),
        strip.background = element_rect(fill = "grey95", colour = "grey30", linewidth = .4),
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0),
        legend.key.size = unit(.9, "lines"))
theme_set(th)

lab_feat <- c(gc_content = "GC content", oe_CG = "CpG O/E", oe_GC = "GpC O/E",
              oe_GG = "GpG O/E", oe_CC = "CpC O/E")

## ---- Fig 5: within vs between decomposition (the key figure) -------------------
dec <- fread("results/decomposition.tsv")
m <- melt(dec, id.vars = c("Subtype", "feature"),
          measure.vars = c("slope_total", "slope_within"),
          variable.name = "model", value.name = "slope")
m[, model := factor(fifelse(model == "slope_total", "total", "within"),
                    levels = c("total", "within"))]
se <- melt(dec, id.vars = c("Subtype", "feature"),
           measure.vars = c("se_total", "se_within"), value.name = "se")
m[, se := se$se]
m[, feature := factor(lab_feat[feature], levels = unname(lab_feat))]

f5 <- ggplot(m, aes(feature, slope, fill = model)) +
  geom_hline(yintercept = 0, linewidth = .4, colour = "grey20") +
  geom_col(position = position_dodge(.7), width = .62, colour = "grey20", linewidth = .25) +
  geom_errorbar(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se),
                position = position_dodge(.7), width = .18, linewidth = .35) +
  facet_wrap(~ paste0("RSV-", Subtype)) +
  scale_fill_manual(values = pal2, name = NULL,
                    labels = c("Total (raw trend)", "Within-lineage (turnover removed)")) +
  scale_y_continuous(labels = label_number(scale = 1e3, accuracy = .1)) +
  labs(title = "Composition trends before and after removing lineage turnover",
       subtitle = "Where the two bars point in opposite directions, the raw trend is driven by lineage replacement, not evolution",
       x = NULL, y = expression("Slope per year " (10^-3)),
       caption = "Error bars = 95% CI. Within-lineage slope is the year coefficient after adding lineage fixed effects.") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
ggsave("results/figures/fig5_within_vs_between.png", f5, width = 9, height = 4.8, dpi = 300)

## ---- Fig 6: per-lineage CpG O/E slopes (forest plot) ---------------------------
sl <- fread("results/per_lineage_slopes.tsv")
sl[, Lineage := factor(Lineage, levels = sl[order(Subtype, slope), Lineage])]
f6 <- ggplot(sl, aes(slope, Lineage, colour = Subtype)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = slope - 1.96 * se, xmax = slope + 1.96 * se),
                 height = 0, linewidth = .5) +
  geom_point(aes(size = n, shape = sig)) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                     labels = c("n.s.", "p < 0.05"), name = NULL) +
  scale_size_continuous(range = c(1.2, 4), name = "sequences") +
  scale_x_continuous(labels = label_number(scale = 1e3, accuracy = .5)) +
  facet_grid(Subtype ~ ., scales = "free_y", space = "free_y") +
  labs(title = "CpG O/E change per year, measured inside each lineage",
       subtitle = "Most lineages sit left of zero: CpG O/E declines as a lineage ages",
       x = expression("Slope per year " (10^-3)), y = NULL,
       caption = "Lineages with >=100 sequences spanning >=4 years. Bars = 95% CI.")
ggsave("results/figures/fig6_per_lineage_slopes.png", f6, width = 7.5, height = 8, dpi = 300)

## ---- Fig 7: lineage growth rate vs composition ---------------------------------
gv <- fread("results/growth_vs_composition.tsv")
cr <- fread("results/growth_composition_corr.tsv")
gl <- melt(gv, id.vars = c("Subtype", "Lineage", "n", "growth_rate"),
           measure.vars = c("oe_CG", "oe_GG"), variable.name = "feature", value.name = "val")
gl[, feature := lab_feat[as.character(feature)]]
lab <- cr[feature %in% c("oe_CG", "oe_GG")]
lab[, feature := lab_feat[feature]]
lab[, txt := sprintf("rho = %+.2f\np = %.3f", rho, p)]

f7 <- ggplot(gl, aes(val, growth_rate)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = .5) +
  geom_point(aes(colour = Subtype, size = n), alpha = .85) +
  geom_text_repel(aes(label = Lineage), size = 2.4, colour = "grey25",
                  max.overlaps = 8, segment.size = .2) +
  geom_text(data = lab, aes(label = txt), x = -Inf, y = Inf,
            hjust = -0.15, vjust = 1.2, size = 3, colour = "grey20") +
  facet_grid(paste0("RSV-", Subtype) ~ feature, scales = "free") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_size_continuous(range = c(1.3, 4.5), name = "sequences") +
  labs(title = "Do faster-growing lineages have different genome composition?",
       subtitle = "Growth rate = logistic change in a lineage's share of its subtype per year",
       x = "Lineage mean O/E ratio", y = "Growth rate (logit / year)",
       caption = "Each point is one lineage. Spearman correlation shown per panel; p-values are uncorrected.")
ggsave("results/figures/fig7_growth_vs_composition.png", f7, width = 9, height = 6.5, dpi = 300)

## ---- Fig 8: robustness of the CpG O/E trend ------------------------------------
rb <- fread("results/robustness.tsv")
rm_ <- melt(rb, id.vars = c("Subtype", "scenario", "n"),
            measure.vars = c("slope_total", "slope_within"),
            variable.name = "model", value.name = "slope")
rse <- melt(rb, id.vars = c("Subtype", "scenario", "n"),
            measure.vars = c("se_total", "se_within"), value.name = "se")
rm_[, se := rse$se]
rm_[, model := factor(fifelse(model == "slope_total", "total", "within"),
                      levels = c("total", "within"))]
rm_[, scenario := factor(scenario, levels = c("full", "excl_2020_2021",
                                              "region_adjusted", "region_balanced"),
                         labels = c("Full data", "Excl. 2020-21",
                                    "Region-adjusted", "Region-balanced"))]
f8 <- ggplot(rm_, aes(slope, scenario, colour = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = slope - 1.96 * se, xmax = slope + 1.96 * se),
                 height = .16, linewidth = .5,
                 position = position_dodge(.55)) +
  geom_point(size = 2.4, position = position_dodge(.55)) +
  facet_wrap(~ paste0("RSV-", Subtype)) +
  scale_colour_manual(values = pal2, name = NULL,
                      labels = c("Total", "Within-lineage")) +
  scale_x_continuous(labels = label_number(scale = 1e3, accuracy = .2)) +
  labs(title = "Robustness of the CpG O/E trend",
       subtitle = "RSV-B's within-lineage decline holds in every scenario; RSV-A's does not",
       x = expression("Slope per year " (10^-3)), y = NULL,
       caption = "Bars = 95% CI. Region-balanced caps each region's contribution within a year.") +
  theme(legend.position = "top")
ggsave("results/figures/fig8_robustness.png", f8, width = 8.5, height = 4.2, dpi = 300)

cat("trend figures written to results/figures/\n")
