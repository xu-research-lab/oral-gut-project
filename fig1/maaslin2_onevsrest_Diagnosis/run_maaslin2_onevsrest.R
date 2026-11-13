# ============================================================
# One-vs-rest MaAsLin2 version of the cancer-type differential-abundance analysis.
#   For each cancer type X (CRC / GC / EC):
#       binary variable  X  vs  (the other two cancer types combined)
#   Same cohort (n=270) and covariates as run_maaslin2_diagnosis_withTx.R.
#   Positive coefficient = higher in X than in the rest (X-enriched).
#   Original pairwise (withTx) version and ANCOM-BC2 version are kept unchanged.
# Run from: D:/oral-gut/code_data_all/fig1
# ============================================================
suppressMessages({library(Maaslin2)})

load("data_for_GI_gut_cancer.RData")
me <- read.csv("meta_extended.csv", row.names = 1, header = TRUE, check.names = FALSE)

stage_col <- names(me)[sapply(me, function(x) {
  v <- as.character(x); ("IV" %in% v) && ("III" %in% v) && !("CRC" %in% v)
})][1]
mmr <- me[rownames(me) %in% sam1$X, c("MSS_MMR", "MSI_Status", stage_col)]
mmr$MMR_MSI <- NA
mmr$MMR_MSI[mmr$MSS_MMR %in% c("MSI_dMMR", "dMMR", "MSI-high") |
              mmr$MSI_Status %in% c("MSI-high", "MSI_dMMR", "dMMR")] <- "dMMR/MSI-H"
mmr$MMR_MSI[mmr$MSS_MMR %in% c("MSS_pMMR", "pMMR", "MSS") |
              mmr$MSI_Status %in% c("MSS", "MSS_pMMR", "pMMR")] <- "pMMR/MSS"
mmr$MMR_MSI[mmr$MSI_Status %in% c("non MSI-high", "MSI-low")] <- "pMMR/MSS"
mmr$Stage_IV <- NA
mmr$Stage_IV[mmr[[stage_col]] == "IV"] <- "IV"
mmr$Stage_IV[mmr[[stage_col]] %in% c("I", "II", "III")] <- "I-III"

sam2 <- sam1
sam2$Stage_IV        <- mmr$Stage_IV[match(sam2$X, rownames(mmr))]
sam2$MMR_MSI         <- mmr$MMR_MSI[match(sam2$X, rownames(mmr))]
sam2$Treatment_class <- me$Treatment_class[match(sam2$X, rownames(me))]
sam2$smoking         <- me$smoking[match(sam2$X, rownames(me))]
sam2 <- droplevels(sam2[sam2$Diagnosis %in% c("CRC", "EC", "GC"), ])
covers <- !is.na(sam2$Stage_IV) & !is.na(sam2$MMR_MSI) & !is.na(sam2$Age) &
  !is.na(sam2$Treatment_class) & !is.na(sam2$Gender) & !is.na(sam2$BMI_category)
sam2 <- sam2[covers, ]
rownames(sam2) <- sam2$X
cat("n =", nrow(sam2), " (", paste(names(table(sam2$Diagnosis)), table(sam2$Diagnosis), sep = "=", collapse = ", "), ")\n")

sam2$Stage_IV        <- factor(sam2$Stage_IV, levels = c("I-III", "IV"))
sam2$MMR_MSI         <- factor(sam2$MMR_MSI, levels = c("pMMR/MSS", "dMMR/MSI-H"))
sam2$Gender          <- factor(sam2$Gender, levels = c("M", "F"))
sam2$BMI_category    <- factor(sam2$BMI_category)
sam2$Treatment_class <- factor(sam2$Treatment_class)
sam2$smoking         <- factor(ifelse(sam2$smoking == 1, "Yes", ifelse(sam2$smoking == 2, "No", NA)),
                               levels = c("No", "Yes"))
sam2$Age             <- as.numeric(sam2$Age)

otu_table2 <- data.frame(t(otu_table1))
otu_table2 <- otu_table2[rownames(sam2), ]

covs <- c("Treatment_class", "Stage_IV", "MMR_MSI", "Gender", "Age", "BMI_category", "smoking")
outdir <- "maaslin2_onevsrest"
dir.create(outdir, showWarnings = FALSE)

fit_one <- function(X) {
  sam2$grp <- factor(ifelse(sam2$Diagnosis == X, X, "Other"), levels = c("Other", X))
  refs <- c("grp,Other", "Treatment_class,ICI", "Stage_IV,I-III",
            "MMR_MSI,pMMR/MSS", "Gender,M", "smoking,No")
  od <- file.path(outdir, paste0("MaAsLin2_", X, "_vs_rest"))
  Maaslin2(input_data = otu_table2, input_metadata = sam2, output = od,
           fixed_effects = c("grp", covs), reference = refs,
           normalization = "NONE", plot_heatmap = FALSE, plot_scatter = FALSE)
  d <- read.table(file.path(od, "all_results.tsv"), header = TRUE, sep = "\t", check.names = FALSE)
  d <- d[d$metadata == "grp" & d$value == X, ]
  data.frame(id = sub("^t__", "", d$feature), group = X,
             coef = d$coef, p = d$pval, q = d$qval, stringsAsFactors = FALSE)
}

res <- do.call(rbind, lapply(c("CRC", "GC", "EC"), fit_one))
write.csv(res, file.path(outdir, "maaslin2_onevsrest_all.csv"), row.names = FALSE)

# enriched group per taxon: among positive (X > rest) significant contrasts, smallest q
pos <- res[!is.na(res$q) & res$coef > 0 & res$q < 0.25, ]
best <- do.call(rbind, lapply(split(pos, pos$id), function(z) z[which.min(z$q), ]))
write.csv(best, file.path(outdir, "maaslin2_onevsrest_best.csv"), row.names = FALSE)

cat("\n== per-cancer one-vs-rest significant (q<0.25, coef>0) ==\n")
print(table(res$group[!is.na(res$q) & res$coef > 0 & res$q < 0.25], useNA = "ifany"))
cat("q<0.25 & coef>0:", nrow(pos), " hits | coef>=2:", sum(pos$coef >= 2),
    " | unique taxa:", nrow(best), "\n")

# ---- plot (Fig1e style) ----
suppressMessages({library(ggplot2); library(dplyr)})
tax <- read.csv("tax_id.csv", row.names = 1, check.names = FALSE)
tax$id <- rownames(tax); tax$Species <- sub("^s__", "", tax$Species)
tx <- tax[, c("id", "Species", "group")]; names(tx)[3] <- "oral_group"

pd <- merge(best, tx, by = "id", all.x = TRUE)
pd <- pd %>%
  filter(coef >= 2, !grepl("GGB", Species)) %>%
  mutate(label_color = ifelse(!is.na(oral_group) & oral_group == "oral", "red", "black"),
         group = factor(group, levels = c("CRC", "GC", "EC"))) %>%
  arrange(group, coef) %>%
  mutate(Taxa = paste(Species, id), Taxa = factor(Taxa, levels = unique(Taxa)))

cat("plotted taxa:", nrow(pd), "\n"); print(table(pd$group, useNA = "ifany"))

p <- ggplot(pd, aes(x = coef, y = Taxa, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(name = "Diagnosis",
                    values = c(CRC = "#66c2a5", GC = "#fc8d63", EC = "#e78ac3"),
                    limits = c("CRC", "GC", "EC")) +
  labs(x = "MaAsLin2 coefficient (vs the other two cancer types)", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, colour = pd$label_color),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())

h <- max(3, nrow(pd) * 0.13)
ggsave(file.path(outdir, "Fig1e_maaslin2_onevsrest.pdf"), p, height = h, width = 5)
ggsave(file.path(outdir, "Fig1e_maaslin2_onevsrest.png"), p, height = h, width = 5, dpi = 150)
cat("DONE\n")
