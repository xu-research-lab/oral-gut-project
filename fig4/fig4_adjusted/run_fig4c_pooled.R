# ============================================================
# Fig4c POOLED analysis (small-sample strategy)
#
# Rationale: subgroup-only panels (esp. CRC PMA, n=17) cannot support
# 6 covariates (EPV ~1.3). We therefore fit a single POOLED model across
# all cancer types and derive stratum-specific (Diagnosis x PMA) Response
# contrasts by borrowing strength, adjusting for all measured host AND
# clinical confounders.
#
#   PRIMARY:  log1p(transmitted.abundance.in.Gut) ~
#               Response_6 * Diagnosis * treatment
#             + Age + Gender + BMI_category + smoking      (host)
#             + MMR_MSI + Stage_IV + Treatment_class       (clinical)
#             + (1 | subject_id)
#
# NOTE on PMA: PMA and non-PMA are paired aliquots of the SAME stool
# sample from the SAME subject (116/119 subjects have both; identical
# Response_6). PMA is therefore a technical/measurement factor, NOT a
# confounder. Because the study design deliberately contrasts PMA vs
# non-PMA, stratum-specific effects are of interest, so the interaction
# model is used as the primary analysis; the formal interaction test and
# the more parsimonious collapsed (diagnosis-level) estimates are also
# reported. The subgroup-only models are shown as sensitivity analyses.
#
# Run from: D:/oral-gut/code_data_all/fig4
# ============================================================
suppressMessages({library(lme4); library(lmerTest); library(ggplot2)})
.progfile <- "C:/Users/huang/AppData/Local/Temp/opencode/fig4c_progress.log"
prog <- function(...) { msg <- paste(..., collapse=" "); line <- paste(format(Sys.time(), "%H:%M:%S"), "|", msg)
  cat(line, "\n"); cat(line, "\n", file=.progfile, append=TRUE); flush.console() }
prog("START")

e <- new.env(); load("D:/oral-gut/code_data_all/fig3/Fig3cdef.RData", envir = e)
features <- get("features", e)
outdir <- "fig4_adjusted"; dir.create(outdir, showWarnings = FALSE)
prog("data loaded")

d <- features
d$subject_id <- gsub("_[ABC]$", "", gsub("^PMA_", "", gsub("^Feces_", "", d$id)))
d$MMR_MSI <- factor(ifelse(d$MSI_Status == "MSI-high", "dMMR",
                    ifelse(is.na(d$MSI_Status), NA, "pMMR")), levels = c("pMMR", "dMMR"))
d$Stage   <- factor(ifelse(d$Tumor_Staging %in% c("IV", "IVA"), "IV", "I-III"), levels = c("I-III", "IV"))
d$Tx      <- factor(ifelse(grepl("combined|chemo|XELOX", d$Treatment), "chemo", "ICI"), levels = c("ICI", "chemo"))
d$BMI_cat <- factor(d$BMI_Status, levels = c("normal", "lean", "Obesity"))

# smoking from the matched-cohort clinical metadata (host factor, available for all 96)
pl <- read.csv("D:/oral-gut/addtional analysis/results/clinical_features/patient_level_metadata.csv",
               check.names = FALSE, stringsAsFactors = FALSE)
d$smoking <- factor(pl$smoking_label[match(d$subject_id, as.character(pl$subject_id))], levels = c("No", "Yes"))

d <- d[d$Diagnosis %in% c("CRC", "GC", "EC") & !is.na(d$Response_6) & d$Response_6 %in% c("R", "NR") &
       !grepl("_B$|_C$", d$id) & d$treatment %in% c("non-PMA", "PMA"), ]
d$Response_6 <- factor(d$Response_6, levels = c("NR", "R"))
d$Diagnosis  <- factor(d$Diagnosis, levels = c("CRC", "GC", "EC"))
d$treatment  <- factor(d$treatment, levels = c("non-PMA", "PMA"))
d$logy       <- log1p(d$transmitted.abundance.in.Gut)

host     <- c("Age", "Gender", "BMI_cat", "smoking")
clinical <- c("MMR_MSI", "Stage", "Tx")
panels   <- expand.grid(Diagnosis = c("CRC", "GC", "EC"), treatment = c("non-PMA", "PMA"),
                        stringsAsFactors = FALSE)

# ---- helper: drop factor covariates constant within a stratum ----
drop_const <- function(x, covs) covs[sapply(covs, function(v) length(unique(na.omit(x[[v]]))) > 1)]

# ============================================================
# 1. PRIMARY pooled model: Response x Diagnosis x treatment
#    -> stratum-specific (Diagnosis x PMA) contrasts
# ============================================================
prim_form <- as.formula(paste("logy ~ Response_6 * Diagnosis * treatment +",
                              paste(c(host, clinical), collapse = " + "),
                              "+ (1 | subject_id)"))
m_prim <- lmer(prim_form, data = d)
prog("primary lmer fit; singular =", isSingular(m_prim))

bn <- names(fixef(m_prim))
mk <- function(spec) { L <- setNames(rep(0, length(bn)), bn); for (n in names(spec)) L[n] <- spec[n]; L }
spec <- list(
  CRC_nonPMA = c(Response_6R = 1),
  CRC_PMA    = c(Response_6R = 1, "Response_6R:treatmentPMA" = 1),
  GC_nonPMA  = c(Response_6R = 1, "Response_6R:DiagnosisGC" = 1),
  GC_PMA     = c(Response_6R = 1, "Response_6R:DiagnosisGC" = 1,
                 "Response_6R:treatmentPMA" = 1, "Response_6R:DiagnosisGC:treatmentPMA" = 1),
  EC_nonPMA  = c(Response_6R = 1, "Response_6R:DiagnosisEC" = 1),
  EC_PMA     = c(Response_6R = 1, "Response_6R:DiagnosisEC" = 1,
                 "Response_6R:treatmentPMA" = 1, "Response_6R:DiagnosisEC:treatmentPMA" = 1))
Lmat <- do.call(rbind, lapply(spec, mk))
ct   <- contest(m_prim, L = Lmat, joint = FALSE)
prim <- data.frame(panel = rownames(ct), n_prim = nobs(m_prim),
                   eff = ct$Estimate, se = ct$`Std. Error`,
                   lo = ct$lower, hi = ct$upper, p = ct$`Pr(>|t|)`,
                   row.names = NULL)

# ---- interaction test (does PMA modify the Response effect?) ----
inter <- summary(m_prim)$coefficients
inter_rows <- grep("Response_6R:treatmentPMA|Response_6R:Diagnosis..:treatmentPMA", rownames(inter))
interaction_test <- data.frame(term = rownames(inter)[inter_rows],
                               est = inter[inter_rows, "Estimate"],
                               p = inter[inter_rows, "Pr(>|t|)"], row.names = NULL)

# ---- secondary: collapsed diagnosis-level (PMA as covariate, no interaction) ----
m_par <- lmer(as.formula(paste("logy ~ Response_6 * Diagnosis + treatment +",
                               paste(c(host, clinical), collapse = " + "),
                               "+ (1 | subject_id)")), data = d)
bnp <- names(fixef(m_par))
mkp <- function(spec) { L <- setNames(rep(0, length(bnp)), bnp); for (n in names(spec)) L[n] <- spec[n]; L }
spec_dx <- list(
  CRC = c(Response_6R = 1),
  GC  = c(Response_6R = 1, "Response_6R:DiagnosisGC" = 1),
  EC  = c(Response_6R = 1, "Response_6R:DiagnosisEC" = 1))
ct_dx <- contest(m_par, L = do.call(rbind, lapply(spec_dx, mkp)), joint = FALSE)
collapsed <- data.frame(panel = rownames(ct_dx), eff = ct_dx$Estimate,
                        lo = ct_dx$lower, hi = ct_dx$upper, p = ct_dx$`Pr(>|t|)`, row.names = NULL)

# ============================================================
# 2. SUBGROUP-ONLY sensitivity analyses (per panel)
# ============================================================
sens <- list()
prog("starting subgroup sensitivity")
for (i in seq_len(nrow(panels))) {
  ct_id <- panels$Diagnosis[i]; tt <- panels$treatment[i]
  x <- d[d$Diagnosis == ct_id & d$treatment == tt, ]
  x <- x[complete.cases(x[, c("logy", "Response_6", host, clinical)]), ]
  used <- drop_const(x, c(host, clinical))
  eff <- function(f, lab) {
    m <- tryCatch(lm(f, data = x), error = function(e) NULL)
    if (is.null(m)) return(data.frame(model = lab, eff = NA, p = NA, n = nrow(x)))
    s <- summary(m)$coef
    if (!"Response_6R" %in% rownames(s)) return(data.frame(model = lab, eff = NA, p = NA, n = nrow(x)))
    data.frame(model = lab, eff = s["Response_6R", 1], p = s["Response_6R", 4], n = nrow(x))
  }
  base <- paste("logy ~ Response_6")
  add <- function(v) if (length(v) == 0) as.formula(base) else as.formula(paste(base, "+", paste(v, collapse = "+")))
  r <- rbind(
    eff(as.formula(base), "unadjusted"),
    eff(add(intersect(clinical, used)), "clinical"),
    eff(add(intersect(host, used)), "host"),
    eff(add(used), "clinical_host"))
  if (length(used) > 0) {
    ps <- tryCatch(glm(as.formula(paste("Response_6 ~", paste(used, collapse = "+"))),
                       data = x, family = binomial), error = function(e) NULL)
    if (!is.null(ps)) {
      x$PS <- fitted(ps)
      mps <- tryCatch(lm(logy ~ Response_6 + PS, data = x), error = function(e) NULL)
      if (!is.null(mps)) {
        s <- summary(mps)$coef
        r <- rbind(r, data.frame(model = "ps", eff = s["Response_6R", 1], p = s["Response_6R", 4], n = nrow(x)))
      }
    }
  }
  r$panel <- paste(ct_id, ifelse(tt == "non-PMA", "nonPMA", "PMA"), sep = "_"); r$n_panel <- nrow(x)
  r$covs_used <- paste(used, collapse = "+")
  sens[[length(sens) + 1]] <- r
}
sens <- do.call(rbind, sens)
prog("sensitivity done")
sens_w <- reshape(sens[, c("panel", "model", "eff", "p", "n_panel", "covs_used")],
                  idvar = "panel", timevar = "model", direction = "wide")
names(sens_w) <- gsub("\\.", "_", names(sens_w))

# ============================================================
# 3. Subject-level cluster bootstrap CI for PRIMARY panel contrasts
# ============================================================
set.seed(20240925); B <- 1000
prog("starting bootstrap B =", B)
subj <- unique(d$subject_id)
boot <- matrix(NA_real_, nrow = B, ncol = nrow(Lmat), dimnames = list(NULL, rownames(Lmat)))
for (b in seq_len(B)) {
  s <- sample(subj, length(subj), replace = TRUE)
  dd <- do.call(rbind, lapply(seq_along(s), function(j) {
    z <- d[d$subject_id == s[j], ]; z$subject_id <- paste0(j, "_", z$subject_id); z }))
  mm <- tryCatch(suppressMessages(lmer(prim_form, data = dd,
           control = lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore"))),
           error = function(e) NULL)
  if (is.null(mm)) next
  bb <- fixef(mm)
  for (k in seq_len(nrow(Lmat))) {
    L <- Lmat[k, ]
    if (all(names(L) %in% names(bb))) boot[b, k] <- sum(L * bb[names(L)])
  }
  if (b %% 100 == 0) prog("  bootstrap", b, "/", B)
}
boot_ci <- data.frame(panel = colnames(boot),
                      boot_lo = apply(boot, 2, quantile, .025, na.rm = TRUE),
                      boot_hi = apply(boot, 2, quantile, .975, na.rm = TRUE),
                      boot_med = apply(boot, 2, median, na.rm = TRUE),
                      boot_n = colSums(!is.na(boot)))

# ============================================================
# 4. Combine and save
# ============================================================
order_p <- c("CRC_nonPMA","CRC_PMA","GC_nonPMA","GC_PMA","EC_nonPMA","EC_PMA")
res <- merge(prim, boot_ci, by = "panel", all.x = TRUE)
res <- res[match(order_p, res$panel), ]
write.csv(res, file.path(outdir, "fig4c_pooled_results.csv"), row.names = FALSE)
write.csv(collapsed, file.path(outdir, "fig4c_pooled_diagnosis_level.csv"), row.names = FALSE)
write.csv(interaction_test, file.path(outdir, "fig4c_pooled_interaction_test.csv"), row.names = FALSE)
sens_w <- sens_w[match(order_p, sens_w$panel), ]
write.csv(sens_w, file.path(outdir, "fig4c_pooled_subgroup_sensitivity.csv"), row.names = FALSE)

writeLines(c("PRIMARY POOLED MODEL (Response x Diagnosis x treatment)",
             paste(deparse(formula(m_prim)), collapse = " "),
             paste0("n = ", nobs(m_prim), "; subjects = ", length(unique(d$subject_id)),
                    "; singular = ", isSingular(m_prim)),
             "",
             capture.output(print(summary(m_prim)$coefficients)),
             "",
             "PRIMARY stratum-specific contrasts (model-based):",
             capture.output(print(ct)),
             "",
             "PMA effect-modification test (Response x treatment):",
             capture.output(print(interaction_test, row.names = FALSE)),
             "",
             "SECONDARY collapsed diagnosis-level contrasts (PMA as covariate):",
             capture.output(print(ct_dx)),
             "",
             "SUBGROUP-ONLY SENSITIVITY (eff | p across models):",
             capture.output(print(sens_w, row.names = FALSE))),
           file.path(outdir, "fig4c_pooled_summary.txt"))

cat("\n==== Fig4c POOLED primary (stratum-specific) + bootstrap ====\n")
print(res[, c("panel","n_prim","eff","se","lo","hi","p","boot_lo","boot_hi","boot_n")],
      row.names = FALSE, digits = 3)
cat("\nPMA (Response x treatment) interaction test:\n")
print(interaction_test, row.names = FALSE)
cat("\nSecondary collapsed diagnosis-level estimates:\n")
print(collapsed, row.names = FALSE, digits = 3)
cat("\nSubgroup-only sensitivity:\n")
print(sens_w[, c("panel","n_panel_unadjusted","eff_unadjusted","p_unadjusted",
                 "eff_clinical","p_clinical","eff_host","p_host",
                 "eff_clinical_host","p_clinical_host","eff_ps","p_ps")],
      row.names = FALSE, digits = 3)

# ============================================================
# 5. Figure A: boxplot (raw data, DESCRIPTIVE) with model p-values
#    + Figure B: forest plot of adjusted effects (matches the p-values)
#    + Caption file describing the adjustment
# ============================================================
ann <- data.frame(
  Diagnosis = rep(c("CRC","GC","EC"), each = 2),
  treatment = rep(c("non-PMA","PMA"), 3), stringsAsFactors = FALSE)
ann$panel <- paste(ann$Diagnosis, ifelse(ann$treatment == "non-PMA", "nonPMA", "PMA"), sep = "_")
ann$Diagnosis <- factor(ann$Diagnosis, levels = c("CRC","GC","EC"))
ann$treatment <- factor(ann$treatment, levels = c("non-PMA","PMA"))
ann <- merge(ann, res[, c("panel","p","boot_lo","boot_hi")], by = "panel", all.x = TRUE)
ann$lab <- ifelse(!is.na(ann$p) & ann$p < 0.05, "*", "")

d$treatment <- factor(d$treatment, levels = c("non-PMA","PMA"))
p <- ggplot(d, aes(x = Diagnosis, y = transmitted.abundance.in.Gut, color = Response_6)) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6), size = 1.6) +
  facet_wrap(~treatment) +
  scale_color_manual(name = "Response", values = c("NR" = "#E41A1C", "R" = "#377EB8")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  geom_text(data = ann, aes(x = Diagnosis, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.4, size = 11.0) +
  labs(x = NULL, y = "Transmitted abundance in gut") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14),
        axis.text = element_text(size = 13), axis.title = element_text(size = 14),
        legend.title = element_text(size = 13), legend.text = element_text(size = 13))
ggsave(file.path(outdir, "fig4c_pooled.png"), p, width = 5.5, height = 3.4, dpi = 300)
tryCatch(ggsave(file.path(outdir, "fig4c_pooled.pdf"), p, width = 5.5, height = 3.4),
         error = function(e) { cat("pdf locked -> fig4c_pooled_new.pdf\n")
           ggsave(file.path(outdir, "fig4c_pooled_new.pdf"), p, width = 5.5, height = 3.4) })

# ---- Forest plot (adjusted effects; consistent with the p-values) ----
fp <- res
fp$panel <- factor(fp$panel, levels = rev(order_p))
fp$sig <- ifelse(fp$p < 0.05, "p < 0.05", "ns")
p2 <- ggplot(fp, aes(x = eff, y = panel, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22) +
  geom_point(aes(size = boot_n), shape = 18) +
  geom_text(aes(label = sprintf("p=%.3f", p)), hjust = -0.15, vjust = 0.0, size = 2.6, color = "grey20") +
  scale_color_manual(values = c("p < 0.05" = "#377EB8", "ns" = "grey40"), name = NULL) +
  scale_size_continuous(guide = "none") +
  expand_limits(x = c(min(fp$lo) - 0.6, max(fp$hi) + 0.9)) +
  labs(x = "Adjusted difference in log1p(transmitted abundance), R vs NR",
       y = NULL,
       title = "Adjusted R vs NR effects (linear mixed model)",
       subtitle = "point = adjusted effect; bars = model 95% CI; p = adjusted p") +
  theme_bw() + theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 7.5))
ggsave(file.path(outdir, "fig4c_pooled_forest.png"), p2, width = 5.0, height = 2.8, dpi = 300)
tryCatch(ggsave(file.path(outdir, "fig4c_pooled_forest.pdf"), p2, width = 5.0, height = 2.8),
         error = function(e) { cat("pdf locked -> fig4c_pooled_forest_new.pdf\n")
           ggsave(file.path(outdir, "fig4c_pooled_forest_new.pdf"), p2, width = 5.0, height = 2.8) })

# ---- Caption / methods text ----
writeLines(c(
  "Figure 4c. Oral-gut transmitted abundance by response to immune checkpoint inhibition.",
  "",
  "Panels are shown separately for non-PMA and PMA preparations.",
  "Boxplots and individual points show the RAW (unadjusted) distribution of",
  "transmitted abundance in gut and are DESCRIPTIVE only.",
  "",
  "The reported p-values and confidence intervals come from a single covariate-",
  "adjusted linear mixed model fitted to the pooled cohort (CRC, GC and EC),",
  "specified as:",
  "  log1p(transmitted abundance in gut) ~",
  "    Response_6 * Diagnosis * treatment(PMA)",
  "    + Age + Gender + BMI category + smoking            (host factors)",
  "    + MMR/MSI + Stage (IV vs I-III) + Treatment class   (clinical factors)",
  "    + (1 | subject_id)                                  (random intercept)",
  "",
  "Adjustment set: host factors (age, sex, body-mass index category and smoking)",
  "and clinical/tumour factors (MMR/MSI status, stage and treatment class). PMA",
  "status is a technical/measurement factor (PMA and non-PMA are paired aliquots",
  "of the same stool sample from the same subject) and is included as a factor in",
  "the model; its interaction with response was tested and was not significant.",
  "",
  "Because several cancer-type x preparation strata are small (e.g. CRC PMA, n=17),",
  "the model was fitted to the pooled cohort so that stratum-specific effects borrow",
  "strength for the residual variance. Each stratum retains its own effect estimate.",
  "Confidence intervals for the primary contrasts were additionally verified by a",
  "subject-level cluster bootstrap (1000 resamples).",
  "",
  "Caveat: for the smallest strata the raw boxplot summarises few observations and",
  "is unstable; the adjusted effect estimates and their CIs (forest plot) are the",
  "basis for inference, not the visual separation of the boxes."),
  file.path(outdir, "fig4c_pooled_caption.txt"))

cat("\n== done: fig4c_pooled_results/diagnosis_level/interaction_test/subgroup_sensitivity ==")
cat("\n         fig4c_pooled.png/pdf (boxplot, raw + adjusted p)")
cat("\n         fig4c_pooled_forest.png/pdf (adjusted effect forest)")
cat("\n         fig4c_pooled_caption.txt (figure legend / adjustment description)\n")
