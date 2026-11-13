# Regenerate PMA / non-PMA adjusted survival curves WITHOUT shaded 95% CI,
# with the legend placed OUTSIDE the panels on the right (for manual adjustment).
# Run from: D:/oral-gut/code_data_all/fig4
suppressMessages({library(survival); library(ggplot2)})

outdir <- "fig4_adjusted"
per <- read.csv(file.path(outdir, "fig4e_per_subject_pfs.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
per <- per[per$Diagnosis %in% c("CRC", "GC", "EC"), ]
per$MMR_MSI <- factor(ifelse(per$MSI_Status == "MSI-high", "dMMR", "pMMR"), levels = c("pMMR", "dMMR"))
per$Stage   <- factor(ifelse(per$Tumor_Staging %in% c("IV", "IVA"), "IV", "I-III"), levels = c("I-III", "IV"))
per$Tx      <- factor(ifelse(grepl("chemo", per$Treatment_class), "chemo", "ICI"), levels = c("ICI", "chemo"))
per$BMI_cat <- factor(per$BMI_Status, levels = c("normal", "lean", "Obesity"))
per$smoking <- factor(per$smoking_label, levels = c("No", "Yes"))
per$Diagnosis <- factor(per$Diagnosis, levels = c("CRC", "GC", "EC"))
host <- c("Age", "Gender", "BMI_cat", "smoking")
clinical <- c("MMR_MSI", "Stage", "Tx")
res <- read.csv(file.path(outdir, "fig4e_pooled_results.csv"), check.names = FALSE)

gcomp_point <- function(dat, fit, d, level, tgrid, base) {
  nd <- base
  nd$Diagnosis <- factor(d, levels = c("CRC", "GC", "EC"))
  nd$shared_level <- factor(level, levels = c("high", "low"))
  s <- summary(survfit(fit, newdata = nd), times = tgrid, extend = TRUE)$surv
  if (is.null(dim(s))) s <- matrix(s, ncol = 1)
  if (ncol(s) == length(tgrid) && nrow(s) != length(tgrid)) s <- t(s)
  rowMeans(s)
}

for (tt in c("non_PMA", "PMA")) {
  d0 <- per[!is.na(per[[paste0("level_", tt)]]), ]
  d0$shared_level <- factor(d0[[paste0("level_", tt)]], levels = c("high", "low"))
  d0 <- d0[complete.cases(d0[, c("time", "event", "shared_level", host, clinical, "Diagnosis")]), ]
  prim_form <- as.formula(paste("Surv(time, event) ~ shared_level * Diagnosis +",
                                paste(c(host, clinical), collapse = " + ")))
  fit <- coxph(prim_form, data = d0, model = TRUE)
  tgrid <- sort(unique(c(0, d0$time[d0$event == 1])))
  base <- d0[!duplicated(d0$SubjectID), c(host, clinical, "Diagnosis")]

  allc <- data.frame(); ann <- data.frame()
  for (ct in c("CRC", "GC", "EC")) {
    allc <- rbind(allc,
      data.frame(time = tgrid, surv = gcomp_point(d0, fit, ct, "high", tgrid, base), group = "high", Diagnosis = ct),
      data.frame(time = tgrid, surv = gcomp_point(d0, fit, ct, "low",  tgrid, base), group = "low",  Diagnosis = ct))
    r <- res[res$exposure == tt & res$Diagnosis == ct, ]
    ann <- rbind(ann, data.frame(
      Diagnosis = ct, x = 10, y = 0.12,
      label = sprintf("%s\nHR=%.2f (%.2f-%.2f)\nP=%.3f", ct, r$HR, r$lo, r$hi, r$p)))
  }
  allc$Diagnosis <- factor(allc$Diagnosis, levels = c("CRC", "GC", "EC"))
  ann$Diagnosis <- factor(ann$Diagnosis, levels = c("CRC", "GC", "EC"))

  p <- ggplot(allc, aes(time, surv, color = group)) +
    geom_step(linewidth = 0.7) +
    facet_wrap(~Diagnosis, nrow = 1) +
    scale_color_manual(name = "Transmission level", values = c(high = "#E41A1C", low = "#377EB8")) +
    geom_text(data = ann, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 3) +
    coord_cartesian(ylim = c(0, 1)) + theme_bw() +
    labs(x = "Days", y = "Progression-free survival") +
    theme(strip.background = element_blank(),
          legend.position = "right",
          legend.background = element_blank(), legend.key = element_blank())

  w <- 470 / 72
  h <- 178.96 / 72
  ggsave(file.path(outdir, paste0("fig4e_pooled_curves_", tt, "_legend_right.pdf")), p, width = w, height = h)
  ggsave(file.path(outdir, paste0("fig4e_pooled_curves_", tt, "_legend_right.png")), p, width = w, height = h, dpi = 300)
  cat("saved", tt, "\n")
}
cat("DONE\n")
