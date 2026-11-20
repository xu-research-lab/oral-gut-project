# ============================================================
# fig2d / fig2e revised: mixed adjusted / unadjusted analysis
#  - cancer-vs-HC/RA  : Wilcoxon (unadjusted; RA/HC have no
#                       clinical metadata, so adjustment is impossible)
#  - cancer-vs-cancer : ANCOVA adjusted for the three clinical
#                       factors  MSI_Status + Tumor_Staging + Treatment_class
#                       (collapsed: dMMR/pMMR, IV vs I-III, chemo vs ICI)
#  - one common Holm correction across all displayed comparisons
# Run from: fig2/
# ============================================================
suppressMessages({library(ggplot2); library(dplyr)})

outdir <- "fig2de_revised"
dir.create(outdir, showWarnings = FALSE)

# ---- 1. data ----
load("shared_abundance_counts.RData")
meta_csv <- "../data/clinical_features/features_filled_with_metadata.csv"
ff <- read.csv(meta_csv, stringsAsFactors = FALSE, check.names = FALSE)

sa <- shared_abundance_counts[, c("id", "treatment", "Diagnosis",
                                  "gut_sum", "oral_sum", "shared_counts")]
names(sa)[names(sa) == "Diagnosis"] <- "saD"
meta <- ff[, c("id", "Diagnosis", "Gender", "Age", "BMI_Status",
               "MSI_Status", "Tumor_Staging", "Treatment", "Treatment_class", "smoking_label")]
d <- merge(sa, meta, by = "id", all.x = TRUE)

# ---- 2. filtering identical to the original fig2def ----
d <- d[!is.na(d$treatment) & d$treatment == "non-PMA", ]
d <- d[!is.na(d$saD) & d$saD != "Treated", ]
d <- d[!grepl("_B$|_C$", d$id), ]
d$Diagnosis <- factor(d$saD, levels = c("CRC", "GC", "EC", "RA", "HC"))
d$log_shared <- log10(d$shared_counts + 1e-6)
d$log_gut    <- log10(d$gut_sum + 1e-6)

# ---- 3. covariates ----
# primary adjustment for the between-cancer (CRC-GC) comparison: HOST factors
d$Age <- as.numeric(d$Age)
d$Gender <- factor(d$Gender)
d$BMI_Status <- factor(d$BMI_Status)
d$smoking_label <- factor(d$smoking_label)
# clinical factors (only for the sensitivity analysis)
d$MMR_MSI  <- factor(ifelse(d$MSI_Status == "MSI-high", "dMMR", "pMMR"), levels = c("pMMR", "dMMR"))
d$Stage_IV <- factor(ifelse(d$Tumor_Staging %in% c("IV", "IVA"), "IV", "I-III"), levels = c("I-III", "IV"))
d$Tx       <- factor(ifelse(grepl("chemo", d$Treatment_class), "chemo", "ICI"), levels = c("ICI", "chemo"))

covs_host  <- c("Age", "Gender", "BMI_Status", "smoking_label")
covs_three <- c("MMR_MSI", "Stage_IV", "Tx")
covs_full  <- c(covs_host, covs_three)   # primary CRC-GC adjustment: host + clinical
cancer <- c("CRC", "GC", "EC")

cat("==== n per group ====\n"); print(table(d$Diagnosis))

# ---- 4. tests ----
wilcox_p <- function(g1, g2, y) {
  a <- d[[y]][d$Diagnosis == g1]; b <- d[[y]][d$Diagnosis == g2]
  wilcox.test(a, b, exact = FALSE)$p.value
}
ancova_p <- function(g1, g2, ylog, covs) {
  sub <- d[d$Diagnosis %in% cancer, ]
  sub <- sub[complete.cases(sub[, c(ylog, "Diagnosis", covs)]), ]
  sub$Diagnosis <- droplevels(sub$Diagnosis)
  form <- as.formula(paste0(ylog, " ~ Diagnosis + ", paste(covs, collapse = " + ")))
  mod <- lm(form, data = sub)
  L <- rep(0, length(coef(mod)))
  i1 <- match(paste0("Diagnosis", g1), names(coef(mod)))
  i2 <- match(paste0("Diagnosis", g2), names(coef(mod)))
  if (!is.na(i1)) L[i1] <- -1
  if (!is.na(i2)) L[i2] <-  1
  s <- summary(multcomp::glht(mod, linfct = t(L)))
  c(p = unname(s$test$pvalues[1]), estimate = unname(s$test$coefficients[1]))
}

analyse <- function(y_raw, y_log, cmps) {
  res <- do.call(rbind, lapply(cmps, function(pp) {
    g1 <- pp[1]; g2 <- pp[2]
    if (g1 %in% cancer && g2 %in% cancer) {
      r <- ancova_p(g1, g2, y_log, covs_full)
      data.frame(group1 = g1, group2 = g2, test = "ANCOVA(adj: host+clinical)",
                 estimate = unname(r["estimate"]), p_raw = unname(r["p"]), stringsAsFactors = FALSE)
    } else {
      data.frame(group1 = g1, group2 = g2, test = "Wilcoxon(unadj)",
                 estimate = NA_real_, p_raw = wilcox_p(g1, g2, y_raw), stringsAsFactors = FALSE)
    }
  }))
  res$p_holm <- p.adjust(res$p_raw, method = "holm")
  # significance based on the (confounder-adjusted) raw p-value, no multiple-testing correction
  res$signif <- ifelse(res$p_raw < 0.001, "***",
                ifelse(res$p_raw < 0.01, "**",
                ifelse(res$p_raw < 0.05, "*", "ns")))
  res
}

cmp_e <- list(c("CRC","HC"), c("EC","HC"), c("GC","HC"),
              c("CRC","GC"), c("EC","RA"), c("RA","GC"))
cmp_d <- list(c("HC","RA"), c("EC","HC"), c("GC","HC"), c("CRC","GC"))

res_e <- analyse("gut_sum", "log_gut", cmp_e)
res_d <- analyse("shared_counts", "log_shared", cmp_d)
cat("\n==== fig2e (gut_sum) ====\n"); print(res_e, row.names = FALSE, digits = 3)
cat("\n==== fig2d (shared_counts) ====\n"); print(res_d, row.names = FALSE, digits = 3)
write.csv(res_e, file.path(outdir, "fig2e_revised_comparisons.csv"), row.names = FALSE)
write.csv(res_d, file.path(outdir, "fig2d_revised_comparisons.csv"), row.names = FALSE)

# ---- 5. sensitivity: CRC-GC adjusted p across covariate sets ----
ancova_contrast <- function(sub, ylog, covs, g1, g2) {
  form <- as.formula(paste0(ylog, " ~ Diagnosis + ", paste(covs, collapse = " + ")))
  mod <- lm(form, data = sub)
  L <- rep(0, length(coef(mod)))
  i1 <- match(paste0("Diagnosis", g1), names(coef(mod)))
  i2 <- match(paste0("Diagnosis", g2), names(coef(mod)))
  if (!is.na(i1)) L[i1] <- -1
  if (!is.na(i2)) L[i2] <-  1
  s <- summary(multcomp::glht(mod, linfct = t(L)))
  p_raw <- unname(s$test$pvalues[1])
  K <- combn(c("CRC","GC","EC"), 2)
  KK <- t(apply(K, 2, function(pp) {
    v <- rep(0, length(coef(mod)))
    a <- match(paste0("Diagnosis", pp[1]), names(coef(mod)))
    b <- match(paste0("Diagnosis", pp[2]), names(coef(mod)))
    if (!is.na(a)) v[a] <- -1; if (!is.na(b)) v[b] <- 1; v
  }))
  st <- summary(multcomp::glht(mod, linfct = KK)); pw <- st$test$pvalues
  names(pw) <- apply(K, 2, function(pp) paste0(pp[1], "-", pp[2]))
  c(p_raw = p_raw, p_tukey = unname(pw[paste0(g1, "-", g2)]))
}
sets <- list(
  none  = character(0),
  three = covs_three,
  host  = c("Age","Gender","BMI_Status","smoking_label"),
  full  = c(covs_three, "Age","Gender","BMI_Status","smoking_label")
)
sens <- data.frame()
for (m in c("log_shared","log_gut")) {
  for (nm in names(sets)) {
    cv <- sets[[nm]]
    sub <- d[d$Diagnosis %in% cancer, ]
    sub <- sub[complete.cases(sub[, c(cv, "Diagnosis")]), ]
    sub$Diagnosis <- droplevels(sub$Diagnosis)
    if (length(cv) == 0) {
      mod0 <- lm(as.formula(paste0(m, " ~ Diagnosis")), sub)
      L <- rep(0, length(coef(mod0)))
      i1 <- match("DiagnosisGC", names(coef(mod0))); i2 <- match("DiagnosisEC", names(coef(mod0)))
      L[i1] <- -1; L[i2] <- 1
      r <- c(p_raw = unname(summary(multcomp::glht(mod0, linfct = t(L)))$test$pvalues[1]), p_tukey = NA)
    } else r <- ancova_contrast(sub, m, cv, "CRC", "GC")
    sens <- rbind(sens, data.frame(metric = m, covset = nm, n = nrow(sub),
                                   CRC_GC_p_raw = unname(r["p_raw"]), CRC_GC_p_tukey = unname(r["p_tukey"])))
  }
}
cat("\n==== sensitivity (CRC vs GC) ====\n"); print(sens, row.names = FALSE, digits = 3)
write.csv(sens, file.path(outdir, "fig2de_revised_sensitivity.csv"), row.names = FALSE)

# ---- 6. plotting (** only significant; no caption; group labels kept) ----
cols <- c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3", "RA" = "#8DA0CB", "HC" = "#A6D854")
add_brackets <- function(p, pd, ycol, res) {
  res <- res[res$p_raw < 0.05, , drop = FALSE]
  if (nrow(res) == 0) return(p)
  ymax <- max(pd[[ycol]], na.rm = TRUE); ymin <- min(pd[[ycol]], na.rm = TRUE); rng <- ymax - ymin
  for (i in seq_len(nrow(res))) {
    y1 <- ymax + rng * (0.05 + 0.10 * (i - 1))
    x1 <- match(res$group1[i], levels(pd$Diagnosis)); x2 <- match(res$group2[i], levels(pd$Diagnosis))
    p <- p +
      annotate("segment", x = x1, xend = x2, y = y1, yend = y1, linewidth = 0.4, linetype = "solid") +
      annotate("segment", x = x1, xend = x1, y = y1, yend = y1 - rng * 0.03, linewidth = 0.4, linetype = "solid") +
      annotate("segment", x = x2, xend = x2, y = y1, yend = y1 - rng * 0.03, linewidth = 0.4, linetype = "solid") +
      annotate("text", x = (x1 + x2) / 2, y = y1 + rng * 0.012, label = res$signif[i], size = 3.4)
  }
  p
}
plot_panel <- function(ycol, ylab, res, fname, w, h) {
  pd <- d[!is.na(d[[ycol]]), ]
  p <- ggplot(pd, aes(Diagnosis, .data[[ycol]])) +
    geom_boxplot(aes(color = Diagnosis), show.legend = FALSE, outlier.shape = NA) +
    geom_jitter(aes(color = Diagnosis), width = 0.05, height = 0, show.legend = FALSE, size = 0.6) +
    theme_bw() + ylab(ylab) + xlab(NULL) + scale_color_manual(values = cols)
  p <- add_brackets(p, pd, ycol, res)
  ggsave(file.path(outdir, fname), plot = p, width = w, height = h); p
}
p_e <- plot_panel("gut_sum", "Cumulative abundance in gut", res_e, "fig2e_revised.pdf", 3, 2.2)
p_d <- plot_panel("shared_counts", "Shared SGB counts", res_d, "fig2d_revised.pdf", 2.8, 2.2)
ggsave(file.path(outdir, "fig2e_revised.png"), p_e, width = 3, height = 2.2, dpi = 300)
ggsave(file.path(outdir, "fig2d_revised.png"), p_d, width = 2.8, height = 2.2, dpi = 300)
cat("\n== done, outputs in", outdir, "==\n")
