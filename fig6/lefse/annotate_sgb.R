# ============================================================
# Shared SGB label annotation.
#
# MetaPhlAn SGB clades look like
#   k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|
#   f__Lachnospiraceae|g__GGB28873|s__GGB28873_SGB41553|t__SGB41553
#
# For a heatmap/bar-chart y-axis label we want something readable.
#   - If the species segment carries a real name (not a MetaPhlAn
#     genome-group code such as GGB####/CFGB####/OFGB####/FGB####, not a
#     "bacterium_*" placeholder, not "unclassified"), use it.
#   - Otherwise walk up the hierarchy (genus -> family -> order -> class ->
#     phylum -> kingdom) and take the first level with a real annotation.
#   - Always append "|SGB###" (pipe separator, so "species|SGB###" reads as
#     "species of SGB###", avoiding the ambiguity of an underscore).
# ============================================================

# A MetaPhlAn segment is "real annotation" unless it is a genome-group /
# placeholder token.
is_real_annotation <- function(seg) {
  if (is.na(seg) || !nzchar(seg)) return(FALSE)
  # strip rank prefix (k__/p__/...) if present
  seg <- sub("^[kpcofgs]__", "", seg)
  if (grepl("unclassified|uncultured|Incertae_sedis", seg, ignore.case = TRUE)) return(FALSE)
  if (grepl("^bacterium_|^bacteria_|_bacterium$", seg, ignore.case = TRUE)) return(FALSE)
  # genome-group style codes: GGB123, CFGB123, OFGB123, FGB123, SGB123, ...GB123
  if (grepl("GB[0-9]+(_SGB[0-9]+)?$", seg)) return(FALSE)
  TRUE
}

# Parse one clade into named rank segments.
parse_clade_segments <- function(clade) {
  parts <- strsplit(clade, "|", fixed = TRUE)[[1]]
  out <- setNames(rep(NA_character_, 9),
                  c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB", "Other"))
  for (p in parts) {
    rk <- substring(p, 1, 3)           # rank prefix is 3 chars: k__/p__/c__/...
    nm <- if (nchar(p) > 3) substring(p, 4) else ""
    r <- switch(rk,
      "k__" = "Kingdom", "p__" = "Phylum", "c__" = "Class", "o__" = "Order",
      "f__" = "Family", "g__" = "Genus", "s__" = "Species", "t__" = "SGB",
      "Other")
    if (!is.null(r) && r != "Other") out[r] <- nm
  }
  out
}

# Build a readable label for one SGB clade.
annotate_sgb_clade <- function(clade) {
  seg <- parse_clade_segments(clade)
  sgb <- seg[["SGB"]]
  if (is.na(sgb) || !grepl("^SGB", sgb)) sgb <- "SGB??"   # not a t__ clade
  rank_order <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  hit <- NA_character_
  for (rk in rank_order) {
    if (!is.na(seg[[rk]]) && is_real_annotation(seg[[rk]])) { hit <- seg[[rk]]; break }
  }
  if (is.na(hit)) hit <- "unclassified"
  paste0(hit, "|", sgb)
}

# Vectorised wrapper: works on full clade strings (SGB rows) or already
# species-level names (non-SGB rows are returned unchanged).
annotate_sgb <- function(clades) {
  vapply(clades, function(cl) {
    if (!grepl("|t__", cl, fixed = TRUE)) return(cl)   # not an SGB clade
    annotate_sgb_clade(cl)
  }, character(1), USE.NAMES = FALSE)
}
