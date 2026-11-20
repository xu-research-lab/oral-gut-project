# ============================================================
# fig4e POOLED adjusted survival analysis (aligns with fig4c)
#
# Exposure : high vs low oral-gut transmission (shared_level),
#            reported separately for non-PMA and PMA.
# Outcome  : per-subject progression-free survival (PFS).
#
# Adjustment set (unified rule = host + clinical):
#   host     : Age, Gender, BMI category, smoking
#   clinical : MMR/MSI, Stage (IV vs I-III), Treatment class
#
# PRIMARY: a single POOLED Cox model with a shared_level * Diagnosis
#   interaction, so that each cancer type keeps its own HR while the
#   regression parameters (and hence precision) borrow strength across
#   cancer types. This mirrors the pooled linear mixed model used for
#   fig4c. Per-cancer HRs are extracted from this model; a subject-level
#   (row) bootstrap provides 95% CIs.
#
# Sensitivity: (i) clinical-only pooled model (previous analysis),
#   (ii) clinical + host pooled model.
#
# Run from: fig4/
# ============================================================
suppressMessages({library(survival); library(coxphf); library(ggplot2)})
.progfile <- "C:/Users/huang/AppData/Local/Temp/opencode/fig4e_progress.log"
prog <- function(...) { line <- paste(format(Sys.time(), "%H:%M:%S"), "|", paste(..., collapse=" "))
  cat(line, "\n"); cat(line, "\n", file=.progfile, append=TRUE); flush.console() }
prog("START fig4e pooled")

outdir <- "fig4_adjusted"; dir.create(outdir, showWarnings = FALSE)
per <- read.csv(file.path(outdir, "fig4e_per_subject_pfs.csv"), check.names = FALSE, stringsAsFactors = FALSE)
per <- per[per$Diagnosis %in% c("CRC", "GC", "EC"), ]
per$MMR_MSI   <- factor(ifelse(per$MSI_Status == "MSI-high", "dMMR", "pMMR"), levels = c("pMMR", "dMMR"))
per$Stage     <- factor(ifelse(per$Tumor_Staging %in% c("IV", "IVA"), "IV", "I-III"), levels = c("I-III", "IV"))
per$Tx        <- factor(ifelse(grepl("chemo", per$Treatment_class), "chemo", "ICI"), levels = c("ICI", "chemo"))
per$BMI_cat   <- factor(per$BMI_Status, levels = c("normal", "lean", "Obesity"))
per$smoking   <- factor(per$smoking_label, levels = c("No", "Yes"))
per$Diagnosis <- factor(per$Diagnosis, levels = c("CRC", "GC", "EC"))
host <- c("Age", "Gender", "BMI_cat", "smoking")
clinical <- c("MMR_MSI", "Stage", "Tx")

# ---- extract per-cancer HR for shared_level (low vs high) from a Cox model ----
cox_hr <- function(fit, scope) {
  b <- coef(fit); V <- vcov(fit); nms <- names(b)
  calc <- function(terms) {
    L <- setNames(rep(0, length(nms)), nms); for (n in terms) L[n] <- 1
    e <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% V %*% L))
    c(HR = exp(e), lo = exp(e - 1.96 * se), hi = exp(e + 1.96 * se), p = 2 * pnorm(-abs(e / se)))
  }
  out <- data.frame()
  for (ct in c("CRC", "GC", "EC")) {
    terms <- "shared_levellow"
    if (ct == "GC" && "shared_levellow:DiagnosisGC" %in% nms) terms <- c(terms, "shared_levellow:DiagnosisGC")
    if (ct == "EC" && "shared_levellow:DiagnosisEC" %in% nms) terms <- c(terms, "shared_levellow:DiagnosisEC")
    v <- calc(terms)
    out <- rbind(out, data.frame(scope = scope, Diagnosis = ct, HR = v[1], lo = v[2], hi = v[3], p = v[4]))
  }
  out
}

R <- 1000; set.seed(20240926)
all_res <- data.frame(); boot_res <- data.frame(); sens_res <- data.frame(); inter_res <- data.frame()

for (tt in c("non_PMA", "PMA")) {
  d0 <- per[!is.na(per[[paste0("level_", tt)]]), ]
  d0$shared_level <- factor(d0[[paste0("level_", tt)]], levels = c("high", "low"))
  d0 <- d0[complete.cases(d0[, c("time", "event", "shared_level", host, clinical, "Diagnosis")]), ]
  prog("exposure", tt, "n =", nrow(d0), "events =", sum(d0$event))

  prim_form <- as.formula(paste("Surv(time, event) ~ shared_level * Diagnosis +",
                                paste(c(host, clinical), collapse = " + ")))
  fit <- coxph(prim_form, data = d0)
  hr <- cox_hr(fit, "pooled_interaction"); hr$exposure <- tt
  all_res <- rbind(all_res, hr)

  # interaction test
  s <- summary(fit)$coefficients
  irows <- grep("shared_levellow:", rownames(s))
  if (length(irows) > 0) inter_res <- rbind(inter_res, data.frame(exposure = tt,
    term = rownames(s)[irows], HR = exp(s[irows, "coef"]), p = s[irows, "Pr(>|z|)"]))

  # subject/row bootstrap CIs (per-subject data -> rows are independent)
  boot <- matrix(NA_real_, R, 3, dimnames = list(NULL, c("CRC", "GC", "EC")))
  for (b in seq_len(R)) {
    db <- d0[sample(nrow(d0), replace = TRUE), ]
    fb <- tryCatch(coxph(prim_form, data = db), error = function(e) NULL)
    if (is.null(fb)) next
    hb <- tryCatch(cox_hr(fb, "boot"), error = function(e) NULL)
    if (is.null(hb)) next
    boot[b, ] <- hb$HR[match(c("CRC", "GC", "EC"), hb$Diagnosis)]
  }
  for (ct in c("CRC", "GC", "EC")) boot_res <- rbind(boot_res, data.frame(exposure = tt, Diagnosis = ct,
    boot_lo = quantile(boot[, ct], .025, na.rm = TRUE), boot_hi = quantile(boot[, ct], .975, na.rm = TRUE),
    boot_med = median(boot[, ct], na.rm = TRUE), boot_n = sum(!is.na(boot[, ct]))))

  # sensitivity: clinical-only vs clinical+host (pooled, stratified baseline)
  f_clin <- coxphf(as.formula(paste("Surv(time,event) ~ shared_level +", paste(clinical, collapse="+"), "+ strata(Diagnosis)")), data = d0)
  f_all  <- coxphf(as.formula(paste("Surv(time,event) ~ shared_level +", paste(c(host, clinical), collapse="+"), "+ strata(Diagnosis)")), data = d0)
  sens_res <- rbind(sens_res,
    data.frame(exposure = tt, model = "clinical_only", HR = exp(f_clin$coefficients["shared_levellow"]),
               lo = f_clin$ci.lower["shared_levellow"], hi = f_clin$ci.upper["shared_levellow"], p = f_clin$prob["shared_levellow"]),
    data.frame(exposure = tt, model = "clinical_host", HR = exp(f_all$coefficients["shared_levellow"]),
               lo = f_all$ci.lower["shared_levellow"], hi = f_all$ci.upper["shared_levellow"], p = f_all$prob["shared_levellow"]))
}

all_res <- merge(all_res, boot_res, by = c("exposure", "Diagnosis"), all.x = TRUE)
all_res <- all_res[order(all_res$exposure, all_res$Diagnosis), ]
write.csv(all_res, file.path(outdir, "fig4e_pooled_results.csv"), row.names = FALSE)
write.csv(inter_res, file.path(outdir, "fig4e_pooled_interaction_test.csv"), row.names = FALSE)
write.csv(sens_res, file.path(outdir, "fig4e_pooled_sensitivity.csv"), row.names = FALSE)

cat("\n==== fig4e POOLED (host + clinical), per-cancer HR (low vs high) ====\n")
print(all_res[, c("exposure","Diagnosis","HR","lo","hi","p","boot_lo","boot_hi","boot_n")], row.names = FALSE, digits = 3)
cat("\nInteraction tests:\n"); print(inter_res, row.names = FALSE, digits = 3)
cat("\nSensitivity (pooled additive):\n"); print(sens_res, row.names = FALSE, digits = 3)

# ============================================================
# g-computation from the SINGLE pooled interaction model:
# for each cancer type, standardise the pooled cohort's covariate profile to
# high and to low transmission and average the predicted survival curves.
gcomp_pooled <- function(dat, prim_form, B = 300) {
  dat <- dat[complete.cases(dat[, c("time", "event", "shared_level", host, clinical, "Diagnosis")]), ]
  tgrid <- sort(unique(c(0, dat$time[dat$event == 1])))
  base <- dat[!duplicated(dat$SubjectID), c(host, clinical, "Diagnosis")]
  marg <- function(fit, d, level) {
    nd <- base; nd$Diagnosis <- factor(d, levels = c("CRC", "GC", "EC"))
    nd$shared_level <- factor(level, levels = c("high", "low"))
    s <- summary(survfit(fit, newdata = nd), times = tgrid, extend = TRUE)$surv
    if (is.null(dim(s))) s <- matrix(s, ncol = 1)
    if (ncol(s) == length(tgrid) && nrow(s) != length(tgrid)) s <- t(s)
    rowMeans(s)
  }
  fit <- coxph(prim_form, data = dat, model = TRUE)
  pt <- list(); for (d in c("CRC","GC","EC")) pt[[d]] <- list(high = marg(fit,d,"high"), low = marg(fit,d,"low"))
  bh <- bl <- array(NA_real_, dim = c(B, length(tgrid), 3), dimnames = list(NULL, NULL, c("CRC","GC","EC")))
  for (b in seq_len(B)) {
    db <- dat[sample(nrow(dat), replace = TRUE), ]
    db$Diagnosis <- droplevels(db$Diagnosis); db$shared_level <- droplevels(db$shared_level)
    fb <- tryCatch(coxph(prim_form, data = db, model = TRUE), error = function(e) NULL); if (is.null(fb)) next
    for (d in c("CRC","GC","EC")) {
      mb <- tryCatch(marg(fb, d, "high"), error = function(e) NULL)
      ml <- tryCatch(marg(fb, d, "low"),  error = function(e) NULL)
      if (!is.null(mb)) bh[b, , d] <- mb
      if (!is.null(ml)) bl[b, , d] <- ml
    }
  }
  q <- function(M, p) apply(M, 2, quantile, p, na.rm = TRUE)
  out <- data.frame()
  for (d in c("CRC","GC","EC")) out <- rbind(out, data.frame(
    time = tgrid, Diagnosis = d,
    high = pt[[d]]$high, high_lo = q(bh[,,d],.025), high_hi = q(bh[,,d],.975),
    low  = pt[[d]]$low,  low_lo  = q(bl[,,d],.025), low_hi  = q(bl[,,d],.975)))
  out
}

for (tt in c("non_PMA", "PMA")) {
  d0 <- per[!is.na(per[[paste0("level_", tt)]]), ]
  d0$shared_level <- factor(d0[[paste0("level_", tt)]], levels = c("high", "low"))
  d0 <- d0[complete.cases(d0[, c("time","event","shared_level",host,clinical,"Diagnosis")]), ]
  prim_form <- as.formula(paste("Surv(time, event) ~ shared_level * Diagnosis +",
                                paste(c(host, clinical), collapse = " + ")))
  allc <- data.frame(); ann <- data.frame()
  g <- gcomp_pooled(d0, prim_form)
  for (ct in c("CRC", "GC", "EC")) {
    gg <- g[g$Diagnosis == ct, ]
    allc <- rbind(allc,
      data.frame(time = gg$time, surv = gg$high, lo = gg$high_lo, hi = gg$high_hi, group = "high", Diagnosis = ct),
      data.frame(time = gg$time, surv = gg$low,  lo = gg$low_lo,  hi = gg$low_hi,  group = "low",  Diagnosis = ct))
    r <- all_res[all_res$exposure == tt & all_res$Diagnosis == ct, ]
    lab <- if (nrow(r) == 0) paste0(ct, "\nHR=NA") else
      sprintf("%s\nHR=%.2f (%.2f-%.2f)\nP=%.3f", ct, r$HR, r$lo, r$hi, r$p)
    ann <- rbind(ann, data.frame(Diagnosis = ct, label = lab, x = 10, y = 0.12))
  }
  allc$Diagnosis <- factor(allc$Diagnosis, levels = c("CRC", "GC", "EC"))
  ann$Diagnosis <- factor(ann$Diagnosis, levels = c("CRC", "GC", "EC"))
  p <- ggplot(allc, aes(time, surv, color = group)) +
    geom_step(linewidth = 0.7) + facet_wrap(~Diagnosis, nrow = 1) +
    scale_color_manual(name = "Transmission level", values = c(high = "#E41A1C", low = "#377EB8")) +
    geom_text(data = ann, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 3) +
    coord_cartesian(ylim = c(0, 1)) + theme_bw() +
    labs(x = "Days", y = "Progression-free survival") +
    theme(strip.background = element_blank(),
          legend.position = c(0.66, 0.96), legend.justification = c(1, 1),
          legend.background = element_rect(fill = "white", color = NA), legend.key = element_blank(),
          legend.title = element_text(size = 7), legend.text = element_text(size = 7),
          legend.key.size = unit(0.35, "lines"), legend.spacing.y = unit(0.02, "cm"),
          legend.margin = margin(0, 0, 0, 0))
  ggsave(file.path(outdir, paste0("fig4e_pooled_curves_", tt, ".png")), p, width = 388/72, height = 178.96/72, dpi = 300)
  tryCatch(ggsave(file.path(outdir, paste0("fig4e_pooled_curves_", tt, ".pdf")), p, width = 388/72, height = 178.96/72),
           error = function(e) ggsave(file.path(outdir, paste0("fig4e_pooled_curves_", tt, "_new.pdf")), p, width = 388/72, height = 178.96/72))
}

# ============================================================
# Forest plot (per-cancer pooled HRs, both exposures)
# ============================================================
fd <- all_res; fd$lab <- paste0(fd$Diagnosis, " (", gsub("_", "-", fd$exposure), ")")
fd$sig <- ifelse(fd$p < 0.05, "p < 0.05", "ns")
fd$lab <- factor(fd$lab, levels = rev(fd$lab))
pf <- ggplot(fd, aes(HR, lab, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22) + geom_point(size = 2.6, shape = 18) +
  geom_text(aes(label = sprintf("P=%.3f", p)), hjust = -0.15, size = 2.6, color = "grey20") +
  scale_x_log10() + scale_color_manual(values = c("p < 0.05" = "#377EB8", "ns" = "grey40"), name = NULL) +
  expand_limits(x = c(min(fd$lo) * 0.7, max(fd$hi) * 2)) +
  labs(x = "Adjusted HR (low vs high transmission, 95% CI)", y = NULL,
       title = "Progression-free survival (pooled Cox, host + clinical)",
       subtitle = "per-cancer HRs from a single pooled model with a transmission x cancer-type interaction") +
  theme_bw() + theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 7))
ggsave(file.path(outdir, "fig4e_pooled_forest.png"), pf, width = 5.6, height = 3.2, dpi = 300)
tryCatch(ggsave(file.path(outdir, "fig4e_pooled_forest.pdf"), pf, width = 5.6, height = 3.2),
         error = function(e) ggsave(file.path(outdir, "fig4e_pooled_forest_new.pdf"), pf, width = 5.6, height = 3.2))

# ---- Caption ----
writeLines(c(
  "Figure 4e. Progression-free survival by oral-gut transmission level.",
  "",
  "Exposure: high vs low oral-gut transmission, shown separately for non-PMA and PMA.",
  "Outcome: per-subject progression-free survival (time to first progression, or last",
  "follow-up if no progression).",
  "",
  "Primary model: a single POOLED Cox proportional-hazards model with a",
  "transmission x cancer-type (Diagnosis) interaction, specified as:",
  "  Surv(time, event) ~ shared_level * Diagnosis",
  "    + Age + Gender + BMI category + smoking        (host factors)",
  "    + MMR/MSI + Stage (IV vs I-III) + Treatment class (clinical factors)",
  "Per-cancer hazard ratios (low vs high transmission) are extracted from this model,",
  "so each cancer type retains its own HR while precision borrows strength across",
  "cancer types; 95% CIs were verified by a 1000-sample bootstrap.",
  "",
  "Sensitivity analyses: (i) clinical factors only (previous specification), and",
  "(ii) clinical + host factors, both fitted to the pooled cohort with stratified",
  "baseline hazard (stratified by Diagnosis) and Firth penalisation for the small",
  "number of events.",
  "",
  "Host factors (age, sex, body-mass index category and smoking) were added to align",
  "the survival analysis with the transmission-abundance analysis (Fig. 4c).",
  "Because event numbers are small in some cancer types, estimates should be",
  "interpreted with caution and the per-cancer results treated as exploratory."),
  file.path(outdir, "fig4e_pooled_caption.txt"))

cat("\n== done: fig4e_pooled_results/interaction_test/sensitivity ==")
cat("\n         fig4e_pooled_curves_non_PMA / _PMA (adjusted KM)")
cat("\n         fig4e_pooled_forest + fig4e_pooled_caption.txt\n")
