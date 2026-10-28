#!/usr/bin/env Rscript
# Script Name: figures_trends.R
# Purpose:     Figures for the within/between-lineage trend split and its robustness checks.
# Inputs:      results/decomposition.tsv, results/robustness.tsv
# Outputs:     results/figures/fig5_within_vs_between.png, fig8_robustness.png
# Notes:       Run from repo root: Rscript R/figures_trends.R

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scales)
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
