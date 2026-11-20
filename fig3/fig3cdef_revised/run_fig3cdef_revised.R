# ============================================================
# Fig3c-f revised: mixed adjusted / unadjusted analysis
#  - comparisons involving RA/HC : Wilcoxon (unadjusted; RA/HC
#    have no clinical metadata, so adjustment is impossible)
#  - cancer-cancer (only CRC-GC): ANCOVA adjusted for host
#    covariates (Age, Gender, BMI_Status, smoking_label)
#  - one common Holm correction across each panel's comparisons
#  - only significant brackets drawn; no caption; group labels kept
# Run from: fig3/
# ============================================================
suppressMessages({library(ggplot2); library(dplyr)})

outdir <- "fig3cdef_revised"
dir.create(outdir, showWarnings = FALSE)

# ---- 1. data ----
e <- new.env(); load("Fig3cdef.RData", envir = e)
features <- get("features", e)

ff <- read.csv("../data/clinical_features/features_filled_with_metadata.csv",
               stringsAsFactors = FALSE, check.names = FALSE)
smoke <- ff[, c("id", "smoking_label")]
d <- merge(features, smoke, by = "id", all.x = TRUE)

# ---- 2. filtering (non-PMA, valid diagnosis, no _B/_C replicates) ----
d <- d[!is.na(d$treatment) & d$treatment == "non-PMA", ]
d <- d[!is.na(d$Diagnosis) & d$Diagnosis != "Treated", ]
d <- d[!grepl("_B$|_C$", d$id), ]
d$Diagnosis <- factor(d$Diagnosis, levels = c("CRC", "GC", "EC", "RA", "HC"))

# ---- 3. covariates ----
d$Treatment_class <- factor(NA, levels = c("ICI","ICI|chemo","ICI|target","ICI|target|chemo","target|chemo","chemo"))
tt <- d$Treatment
d$Treatment_class[grepl("chemo|XELOX", tt)] <- "chemo"
d$Treatment_class[tt %in% c("PD-1", "PD-L1")] <- "ICI"
d$Treatment_class[tt %in% c("PD-1 combined", "PD-L1 combined", "PD-1→PD-1 combined")] <- "ICI|chemo"
d$Stage_IV <- NA
d$Stage_IV[d$Tumor_Staging %in% c("IV", "IVA")] <- "IV"
d$Stage_IV[d$Tumor_Staging %in% c("II","IIA","IIB","IIC","III","IIIB")] <- "I-III"
d$MMR_MSI <- "pMMR/MSS"; d$MMR_MSI[d$MSI_Status == "MSI-high"] <- "dMMR/MSI-H"

d$Age <- as.numeric(d$Age)
d$Gender <- factor(d$Gender)
d$BMI_Status <- factor(d$BMI_Status)
d$smoking_label <- factor(d$smoking_label)
d$MMR_MSI <- factor(d$MMR_MSI)
d$Stage_IV <- factor(d$Stage_IV)

covs_host <- c("Age", "Gender", "BMI_Status", "smoking_label")
covs_full <- c(covs_host, "MMR_MSI", "Stage_IV", "Treatment_class")  # primary CRC-GC adjustment: host + clinical
cancer <- c("CRC", "GC", "EC")

cat("==== n per group (non-PMA, dedup) ====\n"); print(table(d$Diagnosis))
cat("\n==== covariate availability per group ====\n")
print(sapply(covs_host, function(c) tapply(!is.na(d[[c]]), d$Diagnosis, sum)))
cat("\n==== Treatment_class x Diagnosis (cancer) ====\n")
print(table(d$Diagnosis, d$Treatment_class))

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
                 estimate = unname(r["estimate"]), p_raw = unname(r["p"]),
                 stringsAsFactors = FALSE)
    } else {
      data.frame(group1 = g1, group2 = g2, test = "Wilcoxon(unadj)",
                 estimate = NA_real_, p_raw = wilcox_p(g1, g2, y_raw),
                 stringsAsFactors = FALSE)
    }
  }))
  res$p_holm <- p.adjust(res$p_raw, method = "holm")
  # significance based on the (confounder-adjusted) raw p-value, no multiple-testing correction
  res$signif <- ifelse(res$p_raw < 0.001, "***",
                ifelse(res$p_raw < 0.01, "**",
                ifelse(res$p_raw < 0.05, "*", "ns")))
  res
}

panels <- list(
  list(metric = "transmitted.SGB.counts", ylab = "Transmitted strain counts",
       theme = "bw", w = 3, h = 2.5,
       cmps = list(c("EC","HC"), c("GC","HC"), c("CRC","GC"), c("RA","GC"))),
  list(metric = "Oral.gut.shared.species.counts", ylab = NULL,
       theme = "classic", w = 2, h = 2.5,
       cmps = list(c("HC","RA"), c("EC","HC"), c("GC","HC"), c("CRC","GC"), c("CRC","RA"))),
  list(metric = "transmitted.abundance.in.Gut", ylab = "Transmitted abundance",
       theme = "bw", w = 3, h = 2.5,
       cmps = list(c("CRC","HC"), c("EC","HC"), c("GC","HC"), c("CRC","GC"), c("EC","RA"), c("RA","GC"))),
  list(metric = "shared_abundance.Feces", ylab = NULL,
       theme = "classic", w = 2, h = 2.5,
       cmps = list(c("HC","RA"), c("EC","HC"), c("GC","HC"), c("CRC","GC"), c("EC","RA"), c("RA","GC")))
)

# ---- 5. sensitivity: CRC-GC adjusted p under covariate sets ----
ancova_contrast <- function(sub, ylog, covs, g1, g2) {
  form <- as.formula(paste0(ylog, " ~ Diagnosis + ", paste(covs, collapse = " + ")))
  mod <- lm(form, data = sub)
  L <- rep(0, length(coef(mod)))
  i1 <- match(paste0("Diagnosis", g1), names(coef(mod)))
  i2 <- match(paste0("Diagnosis", g2), names(coef(mod)))
  if (!is.na(i1)) L[i1] <- -1
  if (!is.na(i2)) L[i2] <-  1
  p_raw <- unname(summary(multcomp::glht(mod, linfct = t(L)))$test$pvalues[1])
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

all_res <- data.frame(); sens <- data.frame()
sets <- list(
  none          = character(0),
  host          = covs_host,
  host_MSI      = c(covs_host, "MMR_MSI"),
  host_Stage    = c(covs_host, "Stage_IV"),
  full_clinical = c(covs_host, "MMR_MSI", "Stage_IV", "Treatment_class")
)

for (pl in panels) {
  y <- pl$metric
  d$logy <- log10(d[[y]] + 1)
  res <- analyse(y, "logy", pl$cmps)
  res$metric <- y
  all_res <- rbind(all_res, res)
  cat("\n==== ", y, " ====\n"); print(res[, c("metric","group1","group2","test","estimate","p_raw","p_holm","signif")],
                                        row.names = FALSE, digits = 3)
  write.csv(res, file.path(outdir, paste0(gsub("[^A-Za-z0-9]", "_", y), "_comparisons.csv")), row.names = FALSE)

  for (nm in names(sets)) {
    cv <- sets[[nm]]
    sub <- d[d$Diagnosis %in% cancer, ]
    sub <- sub[complete.cases(sub[, c(cv, "Diagnosis")]), ]
    sub$Diagnosis <- droplevels(sub$Diagnosis)
    if (length(cv) == 0) {
      mod0 <- lm(as.formula(paste0("logy ~ Diagnosis")), sub)
      L <- rep(0, length(coef(mod0)))
      i1 <- match("DiagnosisGC", names(coef(mod0))); i2 <- match("DiagnosisEC", names(coef(mod0)))
      L[i1] <- -1; L[i2] <- 1
      r <- c(p_raw = unname(summary(multcomp::glht(mod0, linfct = t(L)))$test$pvalues[1]), p_tukey = NA)
    } else {
      r <- ancova_contrast(sub, "logy", cv, "CRC", "GC")
    }
    sens <- rbind(sens, data.frame(metric = y, covset = nm, n = nrow(sub),
                                   CRC_GC_p_raw = unname(r["p_raw"]),
                                   CRC_GC_p_tukey = unname(r["p_tukey"])))
  }
}
cat("\n==== sensitivity (CRC vs GC, log10(x+1)) ====\n")
print(sens, row.names = FALSE, digits = 3)
write.csv(all_res, file.path(outdir, "fig3cdef_revised_all_comparisons.csv"), row.names = FALSE)
write.csv(sens, file.path(outdir, "fig3cdef_revised_sensitivity.csv"), row.names = FALSE)

# ---- 6. plotting (original style, only significant brackets) ----
cols <- c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3",
          "RA" = "#8da0cd", "HC" = "#a6d854")

add_brackets <- function(p, pd, ycol, res) {
  res <- res[res$p_raw < 0.05, , drop = FALSE]
  if (nrow(res) == 0) return(p)
  ymax <- max(pd[[ycol]], na.rm = TRUE); ymin <- min(pd[[ycol]], na.rm = TRUE)
  rng <- ymax - ymin
  for (i in seq_len(nrow(res))) {
    y1 <- ymax + rng * (0.05 + 0.10 * (i - 1))
    x1 <- match(res$group1[i], levels(pd$Diagnosis)); x2 <- match(res$group2[i], levels(pd$Diagnosis))
    lt <- "solid"
    p <- p +
      annotate("segment", x = x1, xend = x2, y = y1, yend = y1, linewidth = 0.4, linetype = lt) +
      annotate("segment", x = x1, xend = x1, y = y1, yend = y1 - rng * 0.03, linewidth = 0.4, linetype = lt) +
      annotate("segment", x = x2, xend = x2, y = y1, yend = y1 - rng * 0.03, linewidth = 0.4, linetype = lt) +
      annotate("text", x = (x1 + x2) / 2, y = y1 + rng * 0.012, label = res$signif[i], size = 3.2)
  }
  p
}

plot_panel <- function(pl, res) {
  y <- pl$metric
  pd <- d[!is.na(d[[y]]), ]
  p <- ggplot(pd, aes(Diagnosis, .data[[y]])) +
    geom_boxplot(aes(color = Diagnosis), show.legend = FALSE, outlier.shape = NA) +
    geom_jitter(aes(color = Diagnosis), width = 0.1, show.legend = FALSE, size = 0.7) +
    ylab(pl$ylab) + xlab(NULL) +
    scale_color_manual(values = cols)
  p <- if (pl$theme == "classic") {
    p + theme_classic() + theme(panel.grid = element_blank(),
                                axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else {
    p + theme_bw()
  }
  p <- add_brackets(p, pd, y, res)
  fn <- paste0(gsub("[^A-Za-z0-9]", "_", y), "_revised.pdf")
  ggsave(file.path(outdir, fn), plot = p, width = pl$w, height = pl$h)
  ggsave(file.path(outdir, sub("pdf$", "png", fn)), plot = p, width = pl$w, height = pl$h, dpi = 300)
  p
}

for (pl in panels) {
  res <- all_res[all_res$metric == pl$metric, ]
  plot_panel(pl, res)
}

cat("\n== done, outputs in", outdir, "==\n")
