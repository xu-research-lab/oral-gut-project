suppressPackageStartupMessages({
  library(microeco)
  library(data.table)
  library(ggplot2)
})

# ============================================================
# Config
# ============================================================
source("group_config.R")
source("annotate_sgb.R")
root_dir <- "."
out_dir <- "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

abund_file <- "merged_abundance_table.tsv"

# Groups: PD-1 = P1..P5, PD-1 + VP = VP1..VP5 (from group_config.R)
grp1 <- "PD-1"
grp2 <- "PD-1 + VP"
samples_grp1 <- names(group_map)[group_map == grp1]
samples_grp2 <- names(group_map)[group_map == grp2]
samples <- c(samples_grp1, samples_grp2)
cat("PD-1 samples:", paste(samples_grp1, collapse = ","), "\n")
cat("PD-1 + VP samples:", paste(samples_grp2, collapse = ","), "\n")

# ============================================================
# 1. Load SGB clades (t__ rows) + build otu / sample / tax tables
# ============================================================
raw <- fread(abund_file, sep = "\t", header = TRUE, check.names = FALSE,
             colClasses = c("character"), quote = "")
setnames(raw, 1, "clade")
raw <- raw[!grepl("^#", clade)]
raw <- raw[grepl("|t__", clade, fixed = TRUE)]   # SGB clades only

clades <- raw$clade
otu <- as.data.frame(raw[, ..samples])
rownames(otu) <- clades
otu[] <- lapply(otu, function(x) suppressWarnings(as.numeric(x)))
# convert percent -> proportion (sum ~100 per sample)
for (s in samples) otu[[s]] <- otu[[s]] / 100
cat("SGB features:", nrow(otu), "samples:", ncol(otu), "\n")
cat("Per-sample total after normalize (min/max):",
    sprintf("%.4f", min(colSums(otu))), "/", sprintf("%.4f", max(colSums(otu))), "\n")

sample_table <- data.frame(
  SampleID = samples,
  Group = factor(ifelse(samples %in% samples_grp1, grp1, grp2), levels = c(grp1, grp2)),
  row.names = samples,
  stringsAsFactors = FALSE
)

# ============================================================
# 2. tax_table from clade hierarchy (full clade as "Taxonomy")
# ============================================================
parse_clade <- function(cl) {
  parts <- strsplit(cl, "|", fixed = TRUE)[[1]]
  rk <- c(Kingdom = "k__", Phylum = "p__", Class = "c__", Order = "o__",
          Family = "f__", Genus = "g__", Species = "s__", Taxonomy = "t__")
  out <- setNames(rep(NA_character_, length(rk)), names(rk))
  for (i in seq_along(parts)) {
    for (nm in names(rk)) {
      if (startsWith(parts[i], rk[nm])) { out[nm] <- parts[i]; break }
    }
  }
  # Taxonomy column = the full SGB clade (t__ row), Species = s__ part if present
  out["Taxonomy"] <- cl
  if (!is.na(out["Species"])) out["Species"] <- sub("_SGB.*$", "", out["Species"])
  as.data.frame(as.list(out), stringsAsFactors = FALSE)
}
tax_list <- do.call(rbind, lapply(clades, parse_clade))
rownames(tax_list) <- clades
cat("tax_table rows:", nrow(tax_list), " cols:", paste(colnames(tax_list), collapse = ","), "\n")

# ============================================================
# 3. microtable + LEfSe
# ============================================================
mt <- microtable$new(otu_table = otu, sample_table = sample_table, tax_table = tax_list)
mt$sample_sums()
mt$cal_abund(rel = FALSE)

cat("\nRunning LEfSe (PD-1 vs PD-1 + VP, SGB level)...\n")
lefse <- trans_diff$new(
  dataset = mt,
  method = "lefse",
  group = "Group",
  taxa_level = "all",
  alpha = 0.05,
  p_adjust_method = "none",
  lefse_subgroup = NULL,
  lefse_min_subsam = 0
)
cat("LEfSe significant features:", nrow(lefse$res_diff), "\n")
if (!is.null(lefse$res_diff) && nrow(lefse$res_diff) > 0) {
  rd <- as.data.frame(lefse$res_diff)
  rd$Level <- ifelse(grepl("|t__", rd$Taxa, fixed = TRUE), "SGB",
              ifelse(grepl("|s__", rd$Taxa, fixed = TRUE) & !grepl("|t__", rd$Taxa, fixed = TRUE), "Species",
              ifelse(grepl("|g__", rd$Taxa, fixed = TRUE), "Genus",
              ifelse(grepl("|f__", rd$Taxa, fixed = TRUE), "Family",
              ifelse(grepl("|o__", rd$Taxa, fixed = TRUE), "Order",
              ifelse(grepl("|c__", rd$Taxa, fixed = TRUE), "Class",
              ifelse(grepl("|p__", rd$Taxa, fixed = TRUE), "Phylum", "Kingdom")))))))
  write.table(rd, file.path(out_dir, "lefse_res_diff_pd1_vs_vp.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Saved res_diff table\n")
  cat("Significant features by level:\n")
  print(table(rd$Level))
  cat("\nSGB-level significant features:\n")
  sgb_hits <- rd[rd$Level == "SGB", ]
  print(sgb_hits[, c("Taxa", "Group", "LDA", "P.unadj")])
  lefse$res_diff <- rd
} else {
  cat("No significant LEfSe features. Nothing to plot.\n")
  quit(save = "no")
}

# ============================================================
# 4. Plots
# ============================================================
cat("\nPlotting LDA bar chart (all levels)...\n")
p1 <- lefse$plot_diff_bar(use_number = 1:20, threshold = 1, width = 0.6, group_order = c(grp1, grp2))
ggsave(file.path(out_dir, "lefse_LDA_bar_pd1_vs_vp.png"), p1, width = 9, height = 8, dpi = 220)
cat("Saved all-level LDA bar chart\n")

# Cladogram: microeco builds the tree from the k/p/c/o/f/g/s hierarchy; the
# MetaPhlAn t__ (SGB) rank and _SGB###-suffixed species names break ggtree's
# fortify. Attempt it; if it fails, note the limitation (cladogram cannot
# render SGB-level clades) and continue with the SGB-only bar chart.
cat("Plotting cladogram...\n")
p2 <- tryCatch(
  lefse$plot_diff_cladogram(use_taxa_num = 100),
  error = function(e) { cat("cladogram failed:", e$message, "\n"); NULL }
)
if (!is.null(p2)) {
  ggsave(file.path(out_dir, "lefse_cladogram_pd1_vs_vp.png"), p2, width = 12, height = 10, dpi = 220)
  cat("Saved cladogram\n")
}

# ============================================================
# 5. SGB-only LEfSe bar chart (flat taxonomy, tips = SGB clades)
# ============================================================
cat("\nRunning SGB-only LEfSe (flat taxonomy) for a clean SGB-level bar chart...\n")
tax_flat <- data.frame(Species = rownames(otu), stringsAsFactors = FALSE)
rownames(tax_flat) <- rownames(otu)
mt_flat <- microtable$new(otu_table = otu, sample_table = sample_table, tax_table = tax_flat)
mt_flat$sample_sums()
mt_flat$cal_abund(rel = FALSE)
lefse_sgb <- trans_diff$new(
  dataset = mt_flat, method = "lefse", group = "Group", taxa_level = "all",
  alpha = 0.05, p_adjust_method = "none", lefse_subgroup = NULL, lefse_min_subsam = 0
)
if (!is.null(lefse_sgb$res_diff) && nrow(lefse_sgb$res_diff) > 0) {
  rd_sgb <- as.data.frame(lefse_sgb$res_diff)
  rd_sgb$Taxa_full <- rd_sgb$Taxa
  rd_sgb$Taxa <- annotate_sgb(rd_sgb$Taxa)
  write.table(rd_sgb, file.path(out_dir, "lefse_sgb_only_res_diff.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("SGB-only significant features:", nrow(rd_sgb), "\n")
  cat("Plotting SGB-only LDA bar chart (annotated labels)...\n")
  lefse_sgb$res_diff <- rd_sgb
  p3 <- lefse_sgb$plot_diff_bar(use_number = 1:40, threshold = 2, width = 0.6,
                                group_order = c(grp1, grp2), keep_full_name = TRUE,
                                keep_prefix = FALSE)
  ggsave(file.path(out_dir, "lefse_SGB_LDA_bar_pd1_vs_vp.png"), p3, width = 9, height = 9, dpi = 220)
  ggsave(file.path(out_dir, "lefse_SGB_LDA_bar_pd1_vs_vp.pdf"), p3, width = 9, height = 9)
  cat("Saved SGB-only LDA bar chart\n")
} else {
  cat("SGB-only LEfSe returned no features.\n")
}

cat("\nDone!\n")
