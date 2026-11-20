suppressPackageStartupMessages({
  library(ggplot2)
})

# ============================================================
# Mirrored LEfSe LDA bar plot, labels at the BASE of the bars
# (hugging the central axis), NOT at the bar tips. This keeps the
# figure narrow. Bars point to the enriched group.
#   PD-1 enriched      -> bar extends LEFT,  label starts at axis, runs left
#   PD-1 + VP enriched -> bar extends RIGHT, label starts at axis, runs right
# Reads the existing SGB-only LEfSe result table.
# ============================================================
source("group_config.R")
root_dir <- "."
out_dir <- "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grp1 <- "PD-1"
grp2 <- "PD-1 + VP"
gcol <- group_colors

d <- read.table(file.path(out_dir, "lefse_sgb_only_res_diff.tsv"),
                header = TRUE, sep = "\t", check.names = FALSE,
                quote = "", comment.char = "")
d <- d[d$Group %in% c(grp1, grp2), c("Taxa", "Group", "LDA")]
# Shorten the uninformative species name for display (keep the SGB suffix).
d$Taxa <- sub("^Lachnospiraceae_bacterium_MD308", "Lachnospiraceae", d$Taxa)
# Only show taxa with LDA >= 2.
d <- d[d$LDA >= 2, ]
# Highlight the taxa of interest with red labels.
d$lab_color <- ifelse(grepl("Lachnospiraceae|Muribaculum_intestinale|Faecalibaculum_rodentium|Alistipes_shahii", d$Taxa),
                      "red", "grey10")
d$signed_lda <- ifelse(d$Group == grp1, -d$LDA, d$LDA)

# PD-1 block first (bars left), then PD-1 + VP block (bars right);
# within each block sort by LDA descending. First row -> top.
d <- d[order(match(d$Group, c(grp1, grp2)), -d$LDA), ]
d$y <- factor(d$Taxa, levels = rev(d$Taxa))

maxL <- max(d$LDA)
pad  <- 0.03 * maxL
# Names sit on the OPPOSITE side of the axis from their bar:
#   PD-1 (bars left)  -> name just RIGHT of axis, left-aligned, runs right
#   VP   (bars right) -> name just LEFT  of axis, right-aligned, runs left
d$labx <- ifelse(d$Group == grp1, pad, -pad)

# user-specified colours
gcol2 <- c(setNames("#CC79A7", grp1), setNames("#D55E00", grp2))
xlim  <- maxL * 1.06

lab_size <- 2.2

p <- ggplot(d, aes(x = signed_lda, y = y)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_col(aes(fill = Group), width = 0.74) +
  geom_text(data = d[d$Group == grp1, ], aes(x = labx, label = Taxa, color = lab_color),
            hjust = 0, size = lab_size) +
  geom_text(data = d[d$Group == grp2, ], aes(x = labx, label = Taxa, color = lab_color),
            hjust = 1, size = lab_size) +
  scale_color_identity() +
  scale_fill_manual(values = gcol2, name = "") +
  scale_x_continuous(limits = c(-xlim, xlim),
                     breaks = seq(-6, 6, 1), labels = abs) +
  coord_cartesian(clip = "off") +
  labs(x = "LDA score", y = "") +
  theme_bw(11) +
  theme(panel.border = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top")

h <- 6
ggsave(file.path(out_dir, "lefse_SGB_LDA_bar_mirrored_pd1_vs_vp.png"), p,
       width = 3.5, height = h, dpi = 220)
ggsave(file.path(out_dir, "lefse_SGB_LDA_bar_mirrored_pd1_vs_vp.pdf"), p,
       width = 3.5, height = h)
cat("Saved mirrored LEfSe bar plot (", nrow(d), "SGBs, 3.5 x", h, "in )\n")

write.table(d[, c("Taxa", "Group", "LDA", "signed_lda")],
            file.path(out_dir, "lefse_SGB_mirrored_plot_data.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Done!\n")
