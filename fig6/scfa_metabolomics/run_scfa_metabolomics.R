# ============================================================
# Serum SCFA metabolomics — local port of the user's
# scfa_metabolomics_analysis.R (2026-07-29, v1, dataset 1).
# Reads "meta_intensity_fill_plasma.xls" (tab-separated).
# with fallback to data/meta_intensity_fill_plasma.xls.
# Group labels/order/colors unified with group_config.R.
# Group key (user-confirmed 2026-08-03, B/D swap):
#   A=Control, D=PD-1, C=PD-1+VP, B=PD-1+Butyrate,
#   E=PD-1+VP+Butyrate, G=PD-1+E. coli
# Plotting v2 (2026-08-20): every group-comparison bar chart now
# shows ALL C(6,2)=15 pairwise Wilcoxon p-values (red bracket =
# p<0.05, grey = not significant) — not only the significant ones.
# v2 also adds a Total SCFA panel (sum of the 11 SCFA intensities
# per sample) compared across groups with the same all-pairs style.
# ============================================================
suppressMessages({library(ggplot2); library(ggpubr); library(pheatmap); library(vegan); library(rstatix)})

# --- group configuration (group_config.R if present) ----
cfg <- "group_config.R"
if (file.exists(cfg)) {
  source(cfg)
} else {
  group_order <- c("Control", "PD-1", "PD-1 + E. coli", "PD-1 + VP",
                   "PD-1 + Butyrate", "PD-1 + VP + Butyrate")
  group_colors <- c(Control = "#E64B35", "PD-1" = "#4DBBD5",
                    "PD-1 + E. coli" = "#8491B4", "PD-1 + VP" = "#00A087",
                    "PD-1 + Butyrate" = "#3C5488",
                    "PD-1 + VP + Butyrate" = "#F39B7F")
}
gord <- group_order
gcol <- group_colors

# ---------------- input data ----------------
data_file <- "meta_intensity_fill_plasma.xls"

if (!file.exists(data_file)) stop("Data file not found: meta_intensity_fill_plasma.xls")
raw <- read.table(data_file, header = TRUE, sep = "\t", check.names = FALSE, quote = "",
                  comment.char = "", fileEncoding = "UTF-8-BOM")

meta_cols <- colnames(raw)[7:ncol(raw)]
abund <- as.matrix(raw[, meta_cols])
rownames(abund) <- raw$Name
storage.mode(abund) <- "numeric"

# All 39 samples analysed (2026-08-03 final decision), 0 outliers removed.
cat(sprintf("Input: %s\n", data_file))
cat(sprintf("Analysing all %d samples, 0 outliers removed\n", ncol(abund)))

scfa_group_map <- c(A = "Control", B = "PD-1 + Butyrate", C = "PD-1 + VP", D = "PD-1",
                    E = "PD-1 + VP + Butyrate", G = "PD-1 + E. coli")
sample_group <- setNames(scfa_group_map[gsub("_.*", "", colnames(abund))], colnames(abund))

odir <- "."
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Helpers: all-pairs Wilcoxon + bracket annotation
# ============================================================
# Long-format table with every pairwise Wilcoxon test (all group pairs)
pairwise_wilcox <- function(df, adj_method = "BH") {
  glev <- levels(droplevels(df$Group))
  pairs <- combn(glev, 2, simplify = FALSE)
  out <- do.call(rbind, lapply(pairs, function(pr) {
    x <- df$Abundance[df$Group == pr[1]]
    y <- df$Abundance[df$Group == pr[2]]
    wt <- suppressWarnings(wilcox.test(x, y))
    data.frame(Group1 = pr[1], Group2 = pr[2],
               n1 = length(x), n2 = length(y),
               Mean1 = mean(x), Mean2 = mean(y),
               Wilcoxon_p = wt$p.value, stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out$FDR_BH <- p.adjust(out$Wilcoxon_p, method = adj_method)
  out
}

# Draw ALL pairwise comparisons as brackets with p-values above a plot.
# Adjacent pairs are stacked lowest, long-span pairs highest; red = p < 0.05.
add_pairwise_pvals <- function(p, pw, glev, ymax,
                               p_col = "Wilcoxon_p",
                               step_mult = 0.06, txt_size = 2.6,
                               p_dec = 3, plain_line = FALSE,
                               p_prefix = TRUE,
                               line_col = NULL, line_width = 0.35,
                               y_lab_fun = NULL,
                               keep_round = NULL, star_below = NULL) {
  pw <- pw[, c("Group1", "Group2", p_col)]
  ## keep_round: only annotate comparisons whose p-value (rounded to p_dec
  ## decimals) is listed here — used by the focused four-panel figure.
  if (!is.null(keep_round)) {
    keep <- vapply(pw[[p_col]], function(x) {
      is.finite(x) && any(abs(round(x, p_dec) - keep_round) < 1e-8)
    }, logical(1))
    pw <- pw[keep, , drop = FALSE]
  }
  pw$span <- abs(match(pw$Group2, glev) - match(pw$Group1, glev))
  pw <- pw[order(pw$span, pw[[p_col]]), , drop = FALSE]
  n <- nrow(pw)
  if (n == 0) return(p)
  sig <- ifelse(pw[[p_col]] < 0.05, "sig", "ns")
  step <- step_mult * ymax
  tick <- step * 0.35
  hgt  <- ymax + step * seq_len(n)
  seg <- do.call(rbind, lapply(seq_len(n), function(i) {
    x1 <- match(pw$Group1[i], glev); x2 <- match(pw$Group2[i], glev)
    if (plain_line) {
      data.frame(x = x1, xend = x2, y = hgt[i], yend = hgt[i])
    } else {
      data.frame(x = c(x1, x1, x2), xend = c(x2, x1, x2),
                 y = c(hgt[i], hgt[i], hgt[i]),
                 yend = c(hgt[i], hgt[i] - tick, hgt[i] - tick))
    }
  }))
  seg$sig <- rep(sig, each = if (plain_line) 1L else 3L)
  lab <- ifelse(pw[[p_col]] < 10^-p_dec,
                paste0(if (p_prefix) "p < " else "< ", format(10^-p_dec)),
                sprintf(paste0(if (p_prefix) "p = " else "",
                               "%.", p_dec, "f"), pw[[p_col]]))
  ## star_below: below this threshold the numeric label is replaced by "*"
  if (!is.null(star_below)) {
    lab <- ifelse(pw[[p_col]] < star_below, "*", lab)
  }
  txt <- data.frame(
    x = (match(pw$Group1, glev) + match(pw$Group2, glev)) / 2,
    y = hgt + step * 0.5,
    label = lab,
    sig = sig)
  ytop <- ymax + step * (n + 1.5)
  sig_col <- c(sig = "#C0392B", ns = "grey45")
  col_sig <- if (is.null(line_col)) sig_col[["sig"]] else line_col
  col_ns  <- if (is.null(line_col)) sig_col[["ns"]]  else line_col
  p <- p +
    geom_segment(data = seg[seg$sig == "sig", , drop = FALSE],
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = col_sig, linewidth = line_width,
                 inherit.aes = FALSE) +
    geom_segment(data = seg[seg$sig == "ns", , drop = FALSE],
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = col_ns, linewidth = line_width,
                 inherit.aes = FALSE) +
    geom_text(data = txt[txt$sig == "sig", , drop = FALSE],
              aes(x = x, y = y, label = label),
              size = txt_size, color = "black",
              inherit.aes = FALSE) +
    geom_text(data = txt[txt$sig == "ns", , drop = FALSE],
              aes(x = x, y = y, label = label),
              size = txt_size, color = "black",
              inherit.aes = FALSE)
  yscale <- list(limits = c(0, ytop), expand = expansion(mult = c(0, 0.02)))
  if (!is.null(y_lab_fun)) yscale$labels <- y_lab_fun
  p + do.call(scale_y_continuous, yscale)
}

# Pairwise PERMANOVA between every pair of groups
pairwise_permanova <- function(dm, grp, nperm = 999) {
  glev <- unique(as.character(grp))
  pairs <- combn(glev, 2, simplify = FALSE)
  out <- do.call(rbind, lapply(pairs, function(pr) {
    keep <- grp %in% pr
    sub <- as.matrix(dm)[keep, keep]
    g2 <- factor(grp[keep], levels = pr)
    set.seed(1)
    a <- adonis2(as.dist(sub) ~ g2, permutations = nperm)
    data.frame(Group1 = pr[1], Group2 = pr[2],
               R2 = a$R2[1], F = a$F[1], Perm_p = a$`Pr(>F)`[1])
  }))
  out$FDR_BH <- p.adjust(out$Perm_p, method = "BH")
  out
}

# ============================================================
# 1. PCA + PERMANOVA
# ============================================================
pca <- prcomp(t(scale(abund)), center = TRUE, scale. = TRUE)
var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

df_pca <- data.frame(Sample = colnames(abund),
                     Group  = factor(sample_group, levels = gord),
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)

centroids <- aggregate(cbind(PC1, PC2) ~ Group, data = df_pca, FUN = mean)
seg_df <- merge(df_pca[, c("Sample", "Group", "PC1", "PC2")], centroids,
                by = "Group", suffixes = c("", "_centroid"))

dist_mat <- dist(t(scale(abund)))
meta_p <- data.frame(Group = factor(sample_group, levels = gord), row.names = colnames(abund))
set.seed(1)
permanova <- adonis2(dist_mat ~ Group, data = meta_p, permutations = 999)
cat(sprintf("PERMANOVA: R2=%.4f, F=%.3f, p=%.4f\n",
            permanova$R2[1], permanova$F[1], permanova$`Pr(>F)`[1]))

p_pca <- ggplot() +
  geom_segment(data = seg_df,
               aes(PC1, PC2, xend = PC1_centroid, yend = PC2_centroid, color = Group),
               linewidth = 0.25, linetype = "dashed", alpha = 0.5, show.legend = FALSE) +
  geom_point(data = centroids, aes(PC1, PC2, fill = Group),
             size = 5, shape = 23, color = "grey30", stroke = 0.4) +
  geom_point(data = df_pca, aes(PC1, PC2, fill = Group),
             size = 3, shape = 21, color = "grey30", stroke = 0.3) +
  scale_color_manual(values = gcol) +
  scale_fill_manual(values = gcol) +
  labs(title = "PCA – Serum SCFA Metabolomics",
       x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       subtitle = sprintf("PERMANOVA: R\u00b2=%.3f, p=%.3f  |  %d samples, 6 groups",
                          permanova$R2[1], permanova$`Pr(>F)`[1], nrow(df_pca))) +
  theme_bw(13) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))
ggsave(file.path(odir, "scfa_meta_PCA.png"), p_pca, width = 8, height = 6, dpi = 150)
cat("Saved PCA\n")

# 1c. Pairwise PERMANOVA (all 15 group pairs) — saved as a table
pw_permanova <- pairwise_permanova(dist_mat, sample_group)
cat("\nPairwise PERMANOVA (all 15 group pairs):\n")
print(pw_permanova, row.names = FALSE)
write.table(pw_permanova, file.path(odir, "scfa_meta_pairwise_permanova.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved pairwise PERMANOVA table\n")

# ============================================================
# 2. Per-metabolite bar charts (mean ± SD + jitter);
#    ALL 15 pairwise Wilcoxon p-values annotated on every chart
# ============================================================
summary_rows <- list()
plot_list <- list()
pw_combined <- list()
for (i in seq_len(nrow(abund))) {
  met <- rownames(abund)[i]
  v <- abund[i, ]
  df <- data.frame(Sample = names(v),
                   Group = factor(sample_group[names(v)], levels = gord),
                   Abundance = as.numeric(v), stringsAsFactors = FALSE)
  kw <- kruskal.test(Abundance ~ Group, data = df)
  pw <- pairwise_wilcox(df)
  n_sig <- sum(pw$Wilcoxon_p < 0.05)
  enriched <- names(which.max(tapply(as.numeric(v),
                                     factor(sample_group[names(v)], levels = gord),
                                     mean, na.rm = TRUE)))
  summary_rows[[i]] <- data.frame(
    Metabolite = met, KW_chi2 = unname(kw$statistic), KW_df = unname(kw$parameter),
    KW_p = kw$p.value, n_pairs = nrow(pw),
    Sig_pairwise_Wilcoxon = n_sig, Enriched_group = enriched,
    Min_ng_ml = min(as.numeric(v)), Max_ng_ml = max(as.numeric(v)))

  safe_name <- gsub("[ /]", "_", met)
  pw_out <- cbind(Metabolite = met, pw)
  pw_combined[[i]] <- pw_out
  write.table(pw_out, file.path(odir, paste0("scfa_meta_", safe_name, "_pairwise_wilcoxon.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # mean ± SD per group for the bars
  df_sum <- do.call(rbind, lapply(levels(df$Group), function(g) {
    x <- df$Abundance[df$Group == g]
    data.frame(Group = g, Mean = mean(x), SD = sd(x), stringsAsFactors = FALSE)
  }))
  df_sum$Group <- factor(df_sum$Group, levels = gord)
  ymax <- max(c(df$Abundance, df_sum$Mean + df_sum$SD))

  p <- ggplot(df_sum, aes(Group, Mean)) +
    geom_col(width = 0.25, fill = "grey75", color = "grey30", linewidth = 0.4) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.12, linewidth = 0.5) +
    geom_jitter(data = df, aes(Group, Abundance), width = 0.08, size = 1.3,
                alpha = 0.8, shape = 16, color = "black") +
    labs(title = met,
         subtitle = sprintf("Kruskal-Wallis p=%.4f | %d of %d pairwise p<0.05",
                            kw$p.value, n_sig, nrow(pw)),
         x = "", y = "Concentration (ng/mL)") +
    theme_bw(11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  p <- add_pairwise_pvals(p, pw, glev = gord, ymax = ymax)
  plot_list[[i]] <- p
  ggsave(file.path(odir, paste0("scfa_meta_", safe_name, "_bar.png")),
         p, width = 6, height = 7, dpi = 150)
}
summary_tab <- do.call(rbind, summary_rows)
write.table(summary_tab, file.path(odir, "scfa_meta_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Saved %d individual bar charts (all 15 pairwise p-values shown)\n", length(plot_list)))

pw_all <- do.call(rbind, pw_combined)
write.table(pw_all, file.path(odir, "scfa_meta_all_pairwise_wilcoxon.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved combined pairwise-Wilcoxon table\n")


# ============================================================
# 2b. Focused: butyric acid group comparisons (user's key question)
#     Now shows ALL 15 pairwise Wilcoxon p-values
# ============================================================
but <- data.frame(Sample = names(abund["Butyric acid", ]),
                  Group = factor(sample_group, levels = gord),
                  Abundance = as.numeric(abund["Butyric acid", ]))
pw_but <- pairwise_wilcox(but)
cat("\nButyric acid pairwise Wilcoxon, all 15 group pairs (all 39 samples):\n")
print(pw_but, row.names = FALSE)
write.table(pw_but, file.path(odir, "scfa_meta_Butyric_pairwise_wilcoxon.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved butyric pairwise table\n")

but_sum <- do.call(rbind, lapply(levels(but$Group), function(g) {
  x <- but$Abundance[but$Group == g]
  data.frame(Group = g, Mean = mean(x), SD = sd(x), stringsAsFactors = FALSE)
}))
but_sum$Group <- factor(but_sum$Group, levels = gord)
ymax_but <- max(c(but$Abundance, but_sum$Mean + but_sum$SD))

p_but <- ggplot(but_sum, aes(Group, Mean)) +
  geom_col(width = 0.25, fill = "grey75", color = "grey30", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.08, linewidth = 0.5) +
  geom_jitter(data = but, aes(Group, Abundance), width = 0.08, size = 1.3,
              alpha = 0.8, shape = 16, color = "black") +
  labs(title = "Butyric acid by group",
       subtitle = sprintf("Bars = mean, error bars = SD; all 15 pairwise Wilcoxon p-values; red = p<0.05 (%d/%d)",
                          sum(pw_but$Wilcoxon_p < 0.05), nrow(pw_but)),
       x = "", y = "Concentration (ng/mL)") +
  theme_bw(11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 9))
p_but <- add_pairwise_pvals(p_but, pw_but, gord, ymax = ymax_but, txt_size = 3)
ggsave(file.path(odir, "scfa_meta_Butyric_group_comparisons.png"), p_but,
       width = 7.5, height = 6.5, dpi = 150)
cat("Saved butyric group-comparison bar chart (all 15 pairwise p-values)\n")

# ============================================================
# 2c. Total SCFA = sum of all 11 SCFA intensities per sample,
#     compared across groups (all 15 pairwise p-values shown)
# ============================================================
tot <- data.frame(Sample = colnames(abund),
                  Group = factor(sample_group, levels = gord),
                  Abundance = as.numeric(colSums(abund)))
kw_tot <- kruskal.test(Abundance ~ Group, data = tot)
pw_tot <- pairwise_wilcox(tot)
cat("\nTotal SCFA pairwise Wilcoxon, all 15 group pairs (all 39 samples):\n")
print(pw_tot, row.names = FALSE)
write.table(cbind(Metabolite = "Total SCFA", pw_tot),
            file.path(odir, "scfa_meta_Total_SCFA_pairwise_wilcoxon.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved total SCFA pairwise table\n")

tot_sum <- do.call(rbind, lapply(levels(tot$Group), function(g) {
  x <- tot$Abundance[tot$Group == g]
  data.frame(Group = g, Mean = mean(x), SD = sd(x), stringsAsFactors = FALSE)
}))
tot_sum$Group <- factor(tot_sum$Group, levels = gord)
ymax_tot <- max(c(tot$Abundance, tot_sum$Mean + tot_sum$SD))

p_tot <- ggplot(tot_sum, aes(Group, Mean)) +
  geom_col(width = 0.25, fill = "grey75", color = "grey30", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.12, linewidth = 0.5) +
  geom_jitter(data = tot, aes(Group, Abundance), width = 0.08, size = 1.3,
              alpha = 0.8, shape = 16, color = "black") +
  labs(title = "Total SCFA (sum of 11 SCFAs)",
       subtitle = sprintf("Kruskal-Wallis p=%.4f | %d of %d pairwise p<0.05",
                          kw_tot$p.value, sum(pw_tot$Wilcoxon_p < 0.05), nrow(pw_tot)),
       x = "", y = "Concentration (ng/mL)") +
  theme_bw(11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9))
p_tot <- add_pairwise_pvals(p_tot, pw_tot, gord, ymax = ymax_tot)
ggsave(file.path(odir, "scfa_meta_Total_SCFA_group_comparisons.png"), p_tot,
       width = 6, height = 7, dpi = 150)
cat("Saved total SCFA group-comparison bar chart (all 15 pairwise p-values)\n")

# ============================================================
# 2c-2. Focused four-panel SCFA figure (user request 2026-09-09)
#     Acetic acid, Propionic acid, Butyric acid and Total SCFA.
#     Style/group colors match the flow-cytometry figure: hollow bars,
#     SD error bars, and colored jitter points. Only the four user-
#     requested comparisons are annotated in this focused figure.
# ============================================================
flow_group_order <- c("Control", "PD-1", "E. coli + PD-1", "VP + PD-1")
flow_group_colors <- c(
  "Control"        = "#009E73",
  "PD-1"           = "#CC79A7",
  "E. coli + PD-1" = "#56B4E9",
  "VP + PD-1"      = "#D55E00"
)
scfa_to_flow_group <- c(
  "Control"        = "Control",
  "PD-1"           = "PD-1",
  "PD-1 + E. coli" = "E. coli + PD-1",
  "PD-1 + VP"      = "VP + PD-1"
)

## Scientific tick labels "m x 10^e" (plotmath superscript). Keeps the y-axis
## text column narrow and equal width across panels so that panels with large
## concentrations (Acetic acid, Total SCFA: ~40000) are not squeezed.
sci_y_labels <- function(v) {
  parse(text = vapply(v, function(x) {
    if (is.na(x) || x <= 0) return("0")
    e <- floor(log10(x))
    m <- format(x / 10^e)
    paste0(m, " %*% 10^{", e, "}")
  }, character(1)))
}

# Keep only the four biologically targeted comparisons in the focused figure.
scfa_focus_comparisons <- list(
  c("Control", "PD-1"),
  c("PD-1", "E. coli + PD-1"),
  c("PD-1", "VP + PD-1"),
  c("E. coli + PD-1", "VP + PD-1")
)
filter_pw_pairs <- function(pw, pairs) {
  keep <- do.call(rbind, lapply(pairs, function(pr) {
    idx <- (pw$Group1 == pr[1] & pw$Group2 == pr[2]) |
      (pw$Group1 == pr[2] & pw$Group2 == pr[1])
    pw[idx, , drop = FALSE]
  }))
  keep
}

scfa_focus_panel <- function(v, title) {
  df <- data.frame(
    Sample = names(v),
    Group = factor(scfa_to_flow_group[sample_group[names(v)]],
                   levels = flow_group_order),
    Abundance = as.numeric(v),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$Group), , drop = FALSE]
  df_sum <- do.call(rbind, lapply(levels(df$Group), function(g) {
    x <- df$Abundance[df$Group == g]
    data.frame(Group = g, Mean = mean(x),
               SE = sd(x) / sqrt(length(x)), stringsAsFactors = FALSE)
  }))
  df_sum$Group <- factor(df_sum$Group, levels = flow_group_order)
  pw_focus <- filter_pw_pairs(pairwise_wilcox(df), scfa_focus_comparisons)
  ymax_focus <- max(c(df$Abundance, df_sum$Mean + df_sum$SE), na.rm = TRUE)

  cat(sprintf("\n[%s] targeted comparisons:\n", title))
  for (i in seq_len(nrow(pw_focus))) {
    cat(sprintf("  %s vs %s : p = %.4g\n",
                pw_focus$Group1[i], pw_focus$Group2[i], pw_focus$Wilcoxon_p[i]))
  }

  p_focus <- ggplot() +
    geom_col(data = df_sum,
             aes(x = Group, y = Mean, color = Group),
             fill = "white", linewidth = 1.4, width = 0.6) +
    geom_errorbar(data = df_sum,
                  aes(x = Group, ymin = Mean - SE, ymax = Mean + SE,
                      color = Group),
                  width = 0.16, linewidth = 1.0) +
    geom_jitter(data = df,
                aes(x = Group, y = Abundance, color = Group),
                position = position_jitter(width = 0.12, seed = 42),
                size = 1.8, alpha = 0.6) +
    scale_color_manual(values = flow_group_colors, guide = "none") +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, x = NULL, y = "Concentration (ng/mL)") +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", linewidth = 2.0),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text = element_text(size = 14, color = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 14,
                                     color = "black"),
          axis.title = element_text(color = "black"),
          axis.title.y = element_text(size = 16, color = "black"))
  ## Only the P = 0.06 comparisons and the single P = 0.02 result are kept
  ## (user request 2026-09-12); the 0.02 result is shown as "*".
  p_focus <- add_pairwise_pvals(
    p_focus, pw_focus, glev = levels(df$Group), ymax = ymax_focus,
    step_mult = 0.13, txt_size = 5.2,
    p_dec = 2, plain_line = TRUE, p_prefix = FALSE,
    line_col = "black", line_width = 1.0,
    y_lab_fun = sci_y_labels,
    keep_round = c(0.02, 0.06), star_below = 0.05
  )
  p_focus
}

scfa_focus_values <- list(
  "Acetic acid"    = abund["Acetic acid", ],
  "Propionic acid" = abund["Propionic acid", ],
  "Butyric acid"   = abund["Butyric acid", ],
  "Total SCFA"     = colSums(abund)
)
names(scfa_focus_values)[names(scfa_focus_values) == "Total SCFA"] <- "Total SCFA"

p_scfa_focus_list <- lapply(names(scfa_focus_values), function(nm) {
  scfa_focus_panel(scfa_focus_values[[nm]], nm)
})
p_scfa_focus <- ggpubr::ggarrange(plotlist = p_scfa_focus_list,
                                  ncol = 4, nrow = 1)
ggsave(file.path(odir, "scfa_meta_focus_4_metabolites_bar.png"),
       p_scfa_focus, width = 12.6, height = 3.8, dpi = 300)
ggsave(file.path(odir, "scfa_meta_focus_4_metabolites_bar.pdf"),
       p_scfa_focus, width = 12.6, height = 3.8, device = cairo_pdf)
cat("Saved focused four-panel SCFA bar chart (flow-style colors/groups)\n")

# --- extend combined pairwise table with Total SCFA ---
pw_combined[[length(pw_combined) + 1]] <- cbind(Metabolite = "Total SCFA", pw_tot)
pw_all <- do.call(rbind, pw_combined)
write.table(pw_all, file.path(odir, "scfa_meta_all_pairwise_wilcoxon.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- append Total SCFA row to the summary table ---
summary_tab <- rbind(summary_tab, data.frame(
  Metabolite = "Total SCFA",
  KW_chi2 = unname(kw_tot$statistic), KW_df = unname(kw_tot$parameter),
  KW_p = kw_tot$p.value, n_pairs = nrow(pw_tot),
  Sig_pairwise_Wilcoxon = sum(pw_tot$Wilcoxon_p < 0.05),
  Enriched_group = names(which.max(tapply(tot$Abundance, tot$Group, mean))),
  Min_ng_ml = min(tot$Abundance), Max_ng_ml = max(tot$Abundance)))
write.table(summary_tab, file.path(odir, "scfa_meta_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# --- combined figure now includes Total SCFA as the 12th panel ---
plot_list[[length(plot_list) + 1]] <- p_tot
p_all <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 4,
                   common.legend = TRUE, legend = "none")
ggsave(file.path(odir, "scfa_meta_all_metabolites_bar.png"), p_all,
       width = 12, height = 24, dpi = 150)
cat("Saved combined bar chart (12 panels incl. Total SCFA)\n")

# ============================================================
# 2d. Ratio / composition analysis (user request 2026-08-21)
#     Direct butyric acid is NOT significant in PD-1 vs PD-1+VP
#     (Wilcoxon p=0.065, Welch t p=0.087), so explore:
#       (i)  ratios of Acetic / Propionic / Butyric acid
#       (ii) butyric share of total SCFA
#     Key result: Propionic/Butyric ratio is significantly higher
#     in PD-1+VP than PD-1 (3-group KW p=0.007; Dunn BH p.adj=0.005)
#     => butyrate decreases *relative* to propionate.
#     All 15 pairwise Wilcoxon p-values are shown on each chart.
# ============================================================
ac_v <- abund["Acetic acid", ]
pr_v <- abund["Propionic acid", ]
bu_v <- abund["Butyric acid", ]
tot_v <- colSums(abund)

ratio_metrics <- list(
  "Propionic/Butyrate"      = pr_v / bu_v,
  "Butyric/Propionate"      = bu_v / pr_v,
  "Acetic/Butyrate"         = ac_v / bu_v,
  "Acetic/Propionate"       = ac_v / pr_v,
  "Butyric/Total SCFA (%)"  = bu_v / tot_v * 100,
  "Butyric/(Ac+Pr+Bu) (%)"  = bu_v / (ac_v + pr_v + bu_v) * 100
)
ratio_y_lab <- c(
  "Propionic/Butyrate"     = "Ratio (Propionic / Butyric)",
  "Butyric/Propionate"     = "Ratio (Butyric / Propionic)",
  "Acetic/Butyrate"        = "Ratio (Acetic / Butyric)",
  "Acetic/Propionate"      = "Ratio (Acetic / Propionic)",
  "Butyric/Total SCFA (%)" = "Percentage (%)",
  "Butyric/(Ac+Pr+Bu) (%)" = "Percentage (%)"
)

k3 <- c("Control", "PD-1", "PD-1 + VP")
ratio_rows <- list()
ratio_plots <- list()
for (nm in names(ratio_metrics)) {
  v <- ratio_metrics[[nm]]
  df <- data.frame(Sample = names(v),
                   Group = factor(sample_group[names(v)], levels = gord),
                   Abundance = as.numeric(v), stringsAsFactors = FALSE)
  kw6 <- kruskal.test(Abundance ~ Group, data = df)
  df3 <- df[df$Group %in% k3, , drop = FALSE]
  df3$Group <- factor(df3$Group, levels = k3)
  kw3 <- kruskal.test(Abundance ~ Group, data = df3)
  pw <- pairwise_wilcox(df)
  dn6 <- suppressWarnings(dunn_test(df, Abundance ~ Group, p.adjust.method = "BH"))
  dn3 <- suppressWarnings(dunn_test(df3, Abundance ~ Group, p.adjust.method = "BH"))
  row_pd <- which(pw$Group1 == "PD-1" & pw$Group2 == "PD-1 + VP")
  d6_pd <- dn6[dn6$group1 == "PD-1" & dn6$group2 == "PD-1 + VP", , drop = FALSE]
  d3_pd <- dn3[dn3$group1 == "PD-1" & dn3$group2 == "PD-1 + VP", , drop = FALSE]
  ratio_rows[[nm]] <- data.frame(
    Metric = nm,
    Mean_Control = mean(df$Abundance[df$Group == "Control"]),
    Mean_PD1 = mean(df$Abundance[df$Group == "PD-1"]),
    Mean_PD1_VP = mean(df$Abundance[df$Group == "PD-1 + VP"]),
    KW3_p = kw3$p.value, KW6_p = kw6$p.value,
    Wilcoxon_PD1_vs_VP = pw$Wilcoxon_p[row_pd],
    FDR_BH6_PD1_vs_VP = pw$FDR_BH[row_pd],
    Dunn_BH3_PD1_vs_VP = if (nrow(d3_pd) > 0) d3_pd$p.adj else NA_real_,
    Dunn_BH6_PD1_vs_VP = if (nrow(d6_pd) > 0) d6_pd$p.adj else NA_real_)

  safe_name <- gsub("_+$", "", gsub("_+", "_", gsub("[ /()%]", "_", nm)))
  write.table(cbind(Metric = nm, pw),
              file.path(odir, paste0("scfa_meta_ratio_", safe_name, "_pairwise_wilcoxon.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  pw_combined[[length(pw_combined) + 1]] <- cbind(Metabolite = nm, pw)

  df_sum <- do.call(rbind, lapply(levels(df$Group), function(g) {
    x <- df$Abundance[df$Group == g]
    data.frame(Group = g, Mean = mean(x), SD = sd(x), stringsAsFactors = FALSE)
  }))
  df_sum$Group <- factor(df_sum$Group, levels = gord)
  ymax <- max(c(df$Abundance, df_sum$Mean + df_sum$SD))

  p <- ggplot(df_sum, aes(Group, Mean)) +
    geom_col(width = 0.25, fill = "grey75", color = "grey30", linewidth = 0.4) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.12, linewidth = 0.5) +
    geom_jitter(data = df, aes(Group, Abundance), width = 0.08, size = 1.3,
                alpha = 0.8, shape = 16, color = "black") +
    labs(title = nm,
         subtitle = sprintf("Kruskal-Wallis p=%.4f (3 groups) / %.4f (6 groups) | %d of %d pairwise p<0.05",
                            kw3$p.value, kw6$p.value, sum(pw$Wilcoxon_p < 0.05), nrow(pw)),
         x = "", y = ratio_y_lab[[nm]]) +
    theme_bw(11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  p <- add_pairwise_pvals(p, pw, gord, ymax = ymax)
  ratio_plots[[nm]] <- p
  ggsave(file.path(odir, paste0("scfa_meta_ratio_", safe_name, "_bar.png")),
         p, width = 6, height = 7, dpi = 150)
}
ratio_tab <- do.call(rbind, ratio_rows)
write.table(ratio_tab, file.path(odir, "scfa_meta_ratio_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nRatio / composition summary (Control, PD-1, PD-1+VP):\n")
print(format(ratio_tab, digits = 3), row.names = FALSE)
pw_all <- do.call(rbind, pw_combined)
write.table(pw_all, file.path(odir, "scfa_meta_all_pairwise_wilcoxon.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
p_ratio_all <- ggarrange(plotlist = ratio_plots, ncol = 2, nrow = 3)
ggsave(file.path(odir, "scfa_meta_ratios_combined.png"), p_ratio_all,
       width = 12, height = 21, dpi = 150)
cat("Saved 6 ratio/composition bar charts (all 15 pairwise p-values each) + ratio summary table\n")

# ============================================================
# 3. Heatmap (z-score of log10 intensity, group-enriched rows)
# ============================================================
heat_data <- t(scale(t(log10(abund + 1))))
heat_data <- heat_data[!apply(heat_data, 1, anyNA), , drop = FALSE]
col_order <- names(sort(factor(sample_group[colnames(heat_data)], levels = gord)))
heat_data <- heat_data[, col_order, drop = FALSE]

ann_col <- data.frame(Group = factor(sample_group[colnames(heat_data)], levels = gord),
                      row.names = colnames(heat_data))
row_enriched <- sapply(rownames(heat_data), function(m) {
  vals <- abund[m, , drop = TRUE]
  means <- tapply(as.numeric(vals), factor(sample_group[names(vals)], levels = gord),
                  mean, na.rm = TRUE)
  names(which.max(means))
})
row_ann <- data.frame(Enriched = factor(row_enriched, levels = gord),
                      row.names = rownames(heat_data))
ann_colors <- list(Group = gcol, Enriched = gcol)

col_counts <- table(factor(ann_col$Group, levels = gord))
col_gaps <- cumsum(col_counts[col_counts > 0])
row_counts <- table(row_ann$Enriched)
existing_rows <- names(row_counts[row_counts > 0])
row_gaps <- if (length(existing_rows) > 0) cumsum(row_counts[existing_rows]) else numeric(0)

row_order <- order(factor(row_enriched, levels = gord), decreasing = FALSE)
heat_data <- heat_data[row_order, , drop = FALSE]
row_ann   <- row_ann[rownames(heat_data), , drop = FALSE]

n_feat <- nrow(heat_data)
png_h <- max(4, n_feat * 0.6 + 2)
png(file.path(odir, "scfa_meta_heatmap.png"), width = 2400, height = 900 * png_h / 4, res = 300)
pheatmap(heat_data,
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_col = ann_col, annotation_row = row_ann,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
  border_color = NA,
  fontsize = 10, fontsize_row = 10, fontsize_col = 8, angle_col = 45,
  main = "Serum SCFA Profile (z-score of log10 intensity)",
  fontsize_main = 12,
  gaps_row = if (length(row_gaps) > 0) row_gaps[row_gaps < nrow(heat_data)] else NULL,
  gaps_col = col_gaps[col_gaps < ncol(heat_data)]
)
dev.off()
cat("Saved heatmap\n")

cat("\n=== Summary ===\n")
print(summary_tab[, c("Metabolite", "KW_p", "Sig_pairwise_Wilcoxon")], row.names = FALSE)
cat(sprintf("\nPERMANOVA: R\u00b2=%.3f, p=%.3f\n", permanova$R2[1], permanova$`Pr(>F)`[1]))
cat("Done! Output:", odir, "\n")
