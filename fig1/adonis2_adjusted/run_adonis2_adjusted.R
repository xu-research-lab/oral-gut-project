suppressMessages(library(vegan))
suppressMessages(library(dplyr))
set.seed(123)

# Read extended metadata (all covariates, incl. smoking + histology)
meta <- read.csv("meta_extended.csv", row.names = 1, header = TRUE, check.names = FALSE)

# Combine MMR/MSI into one variable: dMMR/MSI-H vs pMMR/MSS
meta$MMR_MSI <- NA
meta$MMR_MSI[meta$MSS_MMR %in% c("MSI_dMMR", "dMMR", "MSI-high") | meta$MSI_Status %in% c("MSI-high", "MSI_dMMR", "dMMR")] <- "dMMR/MSI-H"
meta$MMR_MSI[meta$MSS_MMR %in% c("MSS_pMMR", "pMMR", "MSS") | meta$MSI_Status %in% c("MSS", "MSS_pMMR", "pMMR")] <- "pMMR/MSS"
meta$MMR_MSI[meta$MSI_Status %in% c("non MSI-high", "MSI-low")] <- "pMMR/MSS"

# Stage: binary IV vs I-III (detect the stage column by its values, since the
# column name is Chinese and may be read with a non-UTF8 encoding)
stage_col <- names(meta)[sapply(meta, function(x) {
  v <- as.character(x); ("IV" %in% v) && ("III" %in% v) && !("CRC" %in% v)
})][1]
meta$Stage_IV <- NA
meta$Stage_IV[meta[[stage_col]] == "IV"] <- "IV"
meta$Stage_IV[meta[[stage_col]] %in% c("I", "II", "III")] <- "I-III"

# New covariates
meta$smoking   <- factor(ifelse(meta$smoking == 1, "Yes", ifelse(meta$smoking == 2, "No", NA)),
                         levels = c("No", "Yes"))
meta$histology <- factor(meta$histology)

# Factors
fact_vars <- c("Diagnosis", "Treatment_class", "sex", "BMI_category", "MMR_MSI", "Stage_IV", "smoking", "histology")
for (v in fact_vars) if (v %in% colnames(meta)) meta[[v]] <- as.factor(meta[[v]])
meta$Age <- as.numeric(meta$Age)

dist_files <- c(
  oral_bray = "big_cohort_profile_oral_strains_bray-curtis.tsv",
  oral_jaccard = "big_cohort_profile_oral_strains_jaccard.tsv",
  oral_uwuf = "big_cohort_profile_oral_strains_unweighted-unifrac.tsv",
  oral_wuf = "big_cohort_profile_oral_strains_weighted-unifrac.tsv",
  gut_bray = "big_cohort_profile_bray-curtis.tsv",
  gut_jaccard = "big_cohort_profile_jaccard.tsv",
  gut_uwuf = "big_cohort_profile_unweighted-unifrac.tsv",
  gut_wuf = "big_cohort_profile_weighted-unifrac.tsv"
)

run_adonis <- function(dist_m, meta_sub, terms, label, perms = 999, by = "terms") {
  common <- intersect(rownames(meta_sub), rownames(dist_m))
  meta_s <- meta_sub[common, , drop = FALSE]
  keep_rows <- complete.cases(meta_s[, terms, drop = FALSE])
  meta_s <- meta_s[keep_rows, , drop = FALSE]
  if (nrow(meta_s) == 0) return(NULL)

  covs <- setdiff(terms, "Response_6")
  keep_covs <- covs[sapply(covs, function(v) length(unique(meta_s[[v]])) > 1)]
  rhs <- paste(c(keep_covs, "Response_6"), collapse = " + ")

  common <- rownames(meta_s)
  dist_al <- as.dist(as.matrix(dist_m)[common, common])
  res <- tryCatch(
    adonis2(as.formula(paste("dist_al ~", rhs)), data = meta_s, permutations = perms, by = by),
    error = function(e) { message(label, ": ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)
  data.frame(
    analysis = label, n = nrow(meta_s), term = row.names(res),
    Df = res$Df, R2 = res$R2, `F` = res$F, `Pr(>F)` = res$`Pr(>F)`
  )
}

## Model definitions
# Overall (cross-cancer, Diagnosis focus): between-cancer -> host only
#   (smoking is host; histology is tumour-intrinsic -> NOT in the Diagnosis model).
#   Order: covariates first, then Diagnosis, then Response_6, so that
#   Diagnosis is adjusted for covariates and Response_6 for covariates + Diagnosis.
overall_terms <- c("Treatment_class", "Stage_IV", "MMR_MSI", "Age", "sex", "BMI_category",
                   "smoking", "Diagnosis", "Response_6")
# Within-cancer (R vs NR): host + tumour factors (smoking + histology included)
cancer_terms <- list(
  CRC = c("Treatment_class", "Stage_IV", "Age", "sex", "BMI_category", "smoking", "histology", "Response_6"),
  EC  = c("Treatment_class", "Stage_IV", "Age", "sex", "BMI_category", "smoking", "histology", "Response_6"),
  GC  = c("Treatment_class", "Stage_IV", "MMR_MSI", "Age", "sex", "BMI_category", "smoking", "histology", "Response_6")
)

results <- list()
for (dm in names(dist_files)) {
  dist_m <- read.csv(dist_files[dm], sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
  r <- run_adonis(dist_m, meta, overall_terms, paste0(dm, " | overall"), perms = 999)
  if (!is.null(r)) results[[length(results) + 1]] <- r
  for (ct in c("CRC", "GC", "EC")) {
    sub <- meta[meta$Diagnosis == ct, , drop = FALSE]
    if (nrow(sub) == 0) next
    r <- run_adonis(dist_m, sub, cancer_terms[[ct]], paste0(dm, " | ", ct), perms = 9999)
    if (!is.null(r)) results[[length(results) + 1]] <- r
  }
}

tab <- do.call(rbind, results)
rownames(tab) <- NULL
write.csv(tab, "adonis2_adjusted/all_adonis2_adjusted.csv", row.names = FALSE)
resp_tab <- tab[tab$term == "Response_6", ]
write.csv(resp_tab, "adonis2_adjusted/response_marginal_adjusted.csv", row.names = FALSE)
cat("\n== Response_6 rows ==\n"); print(as.data.frame(resp_tab), digits = 5)
