suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(ggplot2))

outdir  <- "maaslin2_adjusted"
cohorts <- c("All", "CRC", "GC", "EC")

read_term <- function(co) {
  d <- read.table(file.path(outdir, co, "all_results.tsv"), header = TRUE, sep = "\t", check.names = FALSE)
  d <- d[grepl("^Response_6", d$metadata), , drop = FALSE]
  d$group.x <- co
  d
}
comb <- do.call(rbind, lapply(cohorts, read_term))

# tax_id.csv -> SGB species annotation (original label format "Species SGB")
tax <- read.csv("tax_id.csv", row.names = 1, check.names = FALSE)
tax$SGB <- rownames(tax)
ann <- data.frame(
  id      = tax$SGB,
  Kingdom = sub("^k__", "", tax$Kingdom),
  Phylum  = sub("^p__", "", tax$Phylum),
  Class   = sub("^c__", "", tax$Class),
  Order   = sub("^o__", "", tax$Order),
  Family  = sub("^f__", "", tax$Family),
  Genus   = sub("^g__", "", tax$Genus),
  Species = sub("^s__", "", tax$Species),
  group.y = ifelse(!is.na(tax$group) & tax$group == "oral", "shared", "gut native"),
  stringsAsFactors = FALSE
)

comb$id <- sub("^t__", "", comb$feature)
comb <- merge(comb, ann, by = "id", all.x = TRUE)
comb$label <- paste(comb$Species, comb$id)
# Remove GGB placeholder taxa
comb <- comb[!grepl("GGB", comb$Species, ignore.case = TRUE), ]
# Shorten long Candidatus species names (for display labels only)
comb$label <- sub("^Candidatus_Nanosynsacchari_sp_TM7_ANC_38_39_G1_1",
                  "Candidatus_Nanosynsacchari", comb$label)
write.csv(comb, file.path(outdir, "maaslin2_adjusted_response_combined.csv"), row.names = FALSE)

make_fig <- function(stat, thr, tag) {
  # Select taxa: statistic <= threshold in any cohort
  sel <- unique(comb$id[!is.na(comb[[stat]]) & comb[[stat]] <= thr])
  # Exclude SGB6177 (too low prevalence; passed the MaAsLin2 prevalence filter in only some cohorts, so it is incomplete across cohorts)
  sel <- setdiff(sel, "SGB6177")
  pd  <- comb[comb$id %in% sel, ]
  pd$group.x <- factor(pd$group.x, levels = cohorts, ordered = TRUE)

  # Row order: by All-cohort coefficient, ascending (features without an All result go last)
  key <- setNames(pd$coef[pd$group.x == "All"], pd$id[pd$group.x == "All"])
  pd$ordkey <- key[pd$id]
  pd$ordkey[is.na(pd$ordkey)] <- Inf
  pd <- pd[order(pd$ordkey, pd$id, pd$group.x), ]
  pd$label <- factor(pd$label, levels = rev(unique(pd$label)))

  # Add stars: mark * when the threshold for that version is reached
  pd$sig <- ifelse(!is.na(pd[[stat]]) & pd[[stat]] <= thr, "*", "")

  lev  <- levels(pd$label)
  grp  <- pd$group.y[match(lev, pd$label)]
  ycol <- ifelse(grp == "shared", "red", "black")

  p <- ggplot(pd, aes(group.x, label, fill = coef)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    geom_text(aes(label = sig), vjust = 0.75, size = 4) +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8, colour = ycol))

  h <- max(3, length(sel) * 0.13)
  ggsave(file.path(outdir, paste0("Figure1h_maaslin2_adjusted_", tag, ".pdf")), p, height = h, width = 5)
  ggsave(file.path(outdir, paste0("Figure1h_maaslin2_adjusted_", tag, ".png")), p, height = h, width = 5, dpi = 150)
  cat(sprintf("%-6s SGB=%d  rows=%d  starred=%d\n", tag, length(sel), nrow(pd), sum(pd$sig == "*")))
}

make_fig("pval", 0.05, "p005")
make_fig("pval", 0.01, "p001")
make_fig("qval", 0.25, "q025")
make_fig("qval", 0.025, "q0025")
cat("done\n")
