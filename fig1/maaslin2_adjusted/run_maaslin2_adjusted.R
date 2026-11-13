suppressMessages(library(Maaslin2))

load("data_for_GI_gut_cancer.RData")
me <- read.csv("meta_extended.csv", row.names = 1, header = TRUE, check.names = FALSE)

# Merge covariates: Stage_IV and MMR_MSI (from meta_extended); sex from sam1$Gender
stage_col <- names(me)[sapply(me, function(x) { v <- as.character(x); ("IV" %in% v) && ("III" %in% v) && !("CRC" %in% v) })][1]
mmr <- me[rownames(me) %in% sam1$X, c("MSS_MMR", "MSI_Status", stage_col)]
mmr$MMR_MSI <- NA
mmr$MMR_MSI[mmr$MSS_MMR %in% c("MSI_dMMR", "dMMR", "MSI-high") | mmr$MSI_Status %in% c("MSI-high", "MSI_dMMR", "dMMR")] <- "dMMR/MSI-H"
mmr$MMR_MSI[mmr$MSS_MMR %in% c("MSS_pMMR", "pMMR", "MSS") | mmr$MSI_Status %in% c("MSS", "MSS_pMMR", "pMMR")] <- "pMMR/MSS"
mmr$MMR_MSI[mmr$MSI_Status %in% c("non MSI-high", "MSI-low")] <- "pMMR/MSS"
mmr$Stage_IV <- NA
mmr$Stage_IV[mmr[[stage_col]] == "IV"] <- "IV"
mmr$Stage_IV[mmr[[stage_col]] %in% c("I", "II", "III")] <- "I-III"

sam2 <- sam1
sam2$Stage_IV <- mmr$Stage_IV[match(sam2$X, rownames(mmr))]
sam2$MMR_MSI <- mmr$MMR_MSI[match(sam2$X, rownames(mmr))]
sam2$Treatment_class <- me$Treatment_class[match(sam2$X, rownames(me))]
sam2$smoking <- me$smoking[match(sam2$X, rownames(me))]
sam2$histology <- me$histology[match(sam2$X, rownames(me))]
rownames(sam2) <- sam2$X

otu_table2 <- data.frame(t(otu_table1))
rownames(sam2) <- sam2$X
otu_table2 <- otu_table2[rownames(sam2), ]

# Unify factor levels and drop very rare categories to reduce degrees of freedom
sam2$Gender <- as.factor(sam2$Gender)
sam2$Stage_IV <- as.factor(sam2$Stage_IV)
sam2$MMR_MSI <- as.factor(sam2$MMR_MSI)
sam2$smoking <- factor(ifelse(sam2$smoking == 1, "Yes", ifelse(sam2$smoking == 2, "No", NA)))
sam2$histology <- as.factor(sam2$histology)

dir.create("maaslin2_adjusted", showWarnings = FALSE)

# ============================================================
# Unified covariate strategy: include all confounders where possible in every cancer type
#   - Candidate fixed effects: Treatment_class / Stage_IV / MMR_MSI /
#     Gender / smoking / histology
#     (MMR_MSI is derived from MSS_MMR and MSI_Status, using one when the other is missing)
#   - A confounder is automatically dropped if it is constant within a cancer type (only one non-missing level)
#   - Age / BMI_category are random effects (as before)
# ============================================================
all_covs <- c("Treatment_class", "Stage_IV", "MMR_MSI", "Gender", "smoking", "histology")
ref_lookup <- c("Treatment_class,ICI", "Stage_IV,I-III", "MMR_MSI,pMMR/MSS",
                "Gender,M", "smoking,No")

drop_constant <- function(dat, covs) {
  keep <- covs[sapply(covs, function(v) {
    x <- dat[[v]]
    length(unique(na.omit(as.character(x)))) >= 2
  })]
  refs <- ref_lookup[sub(",.*$", "", ref_lookup) %in% keep]
  list(fixed = c("Response_6", keep), reference = refs)
}

run_one <- function(label, dat, random = c("Age", "BMI_category")) {
  spec <- drop_constant(dat, all_covs)
  cat(sprintf("\n[%s] n=%d  fixed=%s\n", label, nrow(dat), paste(spec$fixed, collapse = ", ")))
  Maaslin2(
    input_data = otu_table2,
    input_metadata = dat,
    output = file.path("maaslin2_adjusted", label),
    fixed_effects = spec$fixed,
    random_effects = random,
    reference = spec$reference,
    normalization = "NONE")
}

# Overall: Diagnosis/Age/BMI random; all clinical factors fixed (including MMR_MSI + histology)
run_one("All", sam2[sam2$Response_6 %in% c("R", "NR"), ],
        random = c("Diagnosis", "Age", "BMI_category"))

for (dx in c("CRC", "GC", "EC")) {
  d <- subset(sam2, Diagnosis == dx & Response_6 %in% c("R", "NR"))
  run_one(dx, d)
}

cat("\nDONE\n")
