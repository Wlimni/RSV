#!/usr/bin/env Rscript
# Script Name: figures_codon.R
# Purpose:     Show that CpG absence in RSV genes is not forced by the protein: compare the
#              real CpG count with the range achievable by synonymous recoding.
# Inputs:      results/codon_cpg_analysis.tsv
# Outputs:     results/figures/fig12_synonymous_cpg.png
# Notes:       Rscript R/figures_codon.R

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
pal   <- c(A = "#2c7fb8", B = "#de2d26")
GENES <- c("NS1","NS2","N","P","M","SH","G","F","M2","L")

th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.border = element_rect(linewidth = .4, colour = "grey30"),
        strip.background = element_rect(fill = "grey95", colour = "grey30", linewidth = .4),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))
theme_set(th)

d <- fread("results/codon_cpg_analysis.tsv")
d[, gene := factor(gene, levels = GENES)]

# Panel A: how much CpG the protein would allow, vs how much is actually there
p1 <- ggplot(d, aes(gene)) +
  geom_linerange(aes(ymin = min_syn, ymax = max_syn), linewidth = 3,
                 colour = "grey85", position = position_dodge(.6)) +
  geom_point(aes(y = obs, fill = Subtype), size = 2.8, shape = 21, colour = "grey20",
             position = position_dodge(.6)) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 5, 20, 100, 400, 800)) +
  labs(title = "Real CpG count versus what the protein would allow",
       subtitle = "Grey bar = range achievable by synonymous recoding; dot = the actual virus",
       x = NULL, y = "CpG per gene (log scale)")

# Panel B: position within that range
p2 <- ggplot(d, aes(gene, frac_of_range, fill = Subtype)) +
  geom_col(position = position_dodge(.72), width = .64,
           colour = "grey20", linewidth = .22) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Position within the achievable range",
       subtitle = "0% = every removable CpG has been removed. All genes sit near the floor.",
       x = NULL, y = "Position in range")

# Panel C: the share of genomes sitting exactly at the floor - this is what singles out NS1
p3 <- ggplot(d, aes(gene, pct_at_min, fill = Subtype)) +
  geom_col(position = position_dodge(.72), width = .64,
           colour = "grey20", linewidth = .22) +
  scale_fill_manual(values = pal) +
  labs(title = "Genomes carrying no removable CpG at all",
       subtitle = "Only RSV-A NS1/NS2 and RSV-B SH are driven completely to zero",
       x = "Gene (genome order 3' to 5')", y = "% of genomes at the floor")

f12 <- (p1 / p2 / p3) + plot_layout(guides = "collect") +
  plot_annotation(caption = "600 year-stratified genomes per subtype. Synonymous recoding keeps the protein identical.")
ggsave("results/figures/fig12_synonymous_cpg.png", f12, width = 8.5, height = 10, dpi = 300)
cat("fig12 written\n")
