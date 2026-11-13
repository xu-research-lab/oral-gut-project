suppressMessages(library(tidyverse))
suppressMessages(library(gghalves))
suppressMessages(library(ggpubr))
suppressMessages(library(patchwork))

outdir <- "fig1bfg_adjusted"
dir.create(outdir, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────
e <- new.env()
load("SGB_counts_abundance.RData", envir = e)
df <- as.data.frame(e$SGB_counts_abundance)
df$id <- as.character(df$id)

meta <- read.csv("meta_extended.csv", check.names = FALSE)
names(meta)[names(meta) == "index"] <- "id"
meta$id <- as.character(meta$id)

# ── 2. Build covariates (identical to adonis2_adjusted)─────────
meta$MMR_MSI <- NA
meta$MMR_MSI[meta$MSS_MMR %in% c("MSI_dMMR", "dMMR", "MSI-high") | meta$MSI_Status %in% c("MSI-high", "MSI_dMMR", "dMMR")] <- "dMMR/MSI-H"
meta$MMR_MSI[meta$MSS_MMR %in% c("MSS_pMMR", "pMMR", "MSS") | meta$MSI_Status %in% c("MSS", "MSS_pMMR", "pMMR")] <- "pMMR/MSS"
meta$MMR_MSI[meta$MSI_Status %in% c("non MSI-high", "MSI-low")] <- "pMMR/MSS"

meta$Stage_IV <- NA
meta$Stage_IV[meta$"分期" == "IV"]        <- "IV"
meta$Stage_IV[meta$"分期" %in% c("I", "II", "III")] <- "I-III"
meta$Stage_IV <- factor(meta$Stage_IV, levels = c("I-III", "IV"))

meta$MMR_MSI <- as.factor(meta$MMR_MSI)
meta$Treatment_class <- as.factor(meta$Treatment_class)
meta$sex <- as.factor(meta$sex)
meta$BMI_category <- as.factor(meta$BMI_category)
meta$Age <- as.numeric(meta$Age)
meta$smoking <- factor(ifelse(meta$smoking == 1, "Yes", ifelse(meta$smoking == 2, "No", NA)))
meta$histology <- factor(meta$histology)

# ── 3. Merge ──────────────────────────────────────────────
df <- merge(df, meta[, c("id", "Treatment_class", "Stage_IV", "MMR_MSI", "Age", "sex", "BMI_category", "smoking", "histology")],
            by = "id", all.x = TRUE)

# ── 4. Define model formulas ──────────────────────────────────────
# Overall R vs NR (across cancer types, focus=Response_6): host + tumour factors + Diagnosis
overall_form <- function(metric) {
  as.formula(paste0(metric, " ~ Response_6 + Diagnosis + Treatment_class + Stage_IV + MMR_MSI + Age + sex + BMI_category + smoking + histology"))
}
# Within-cancer R vs NR: host + tumour factors (including smoking + histology)
cancer_form <- function(metric, dx) {
  base <- paste0(metric, " ~ Response_6 + Treatment_class + Stage_IV + Age + sex + BMI_category + smoking + histology")
  if (dx == "GC") base <- paste0(base, " + MMR_MSI")
  as.formula(base)
}

# ── 5. Fit models and extract Response_6 p-values ───────────────────────────
metrics <- c("total observed", "oral-origin", "cumlative abundance")
metrics_log <- c("log_total", "log_oral", "log_cum")
# log transform (cumulative abundance is >0; counts are >=7; add 1e-6 to avoid log(0))
for (i in seq_along(metrics)) {
  df[[metrics_log[i]]] <- log10(df[[metrics[i]]] + 1e-6)
}

drop_const_terms <- function(form, data) {
  tv <- attr(terms(form), "term.labels")
  tv <- tv[sapply(tv, function(v) is.numeric(data[[v]]) || length(unique(na.omit(data[[v]]))) > 1)]
  as.formula(paste(deparse(form[[2]]), "~", paste(tv, collapse = " + ")))
}
get_resp_p <- function(mod) {
  # Use car::Anova Type II if available; otherwise use base summary
  if (requireNamespace("car", quietly = TRUE)) {
    a <- car::Anova(mod, type = "II")
    p <- a["Response_6", "Pr(>F)"]
  } else {
    s <- summary(mod)
    p <- coef(s)["Response_6R", "Pr(>|t|)"]
  }
  if (is.na(p)) p <- 1
  return(p)
}

all_results <- data.frame()

for (m in metrics_log) {
  metric_label <- gsub("log_", "", m)
  # ── Overall model
  form_all <- as.formula(paste0(m, " ~ Response_6 + Diagnosis + Treatment_class + Stage_IV + MMR_MSI + Age + sex + BMI_category + smoking + histology"))
  mod_all <- lm(drop_const_terms(form_all, df), data = df)
  p_all <- get_resp_p(mod_all)
  all_results <- rbind(all_results, data.frame(metric = metric_label, level = "Overall", p_adj = p_all, stringsAsFactors = FALSE))
  
  # ── Within cancer type
  for (dx in c("CRC", "GC", "EC")) {
    sub <- df[df$Diagnosis == dx, ]
    form_ok <- drop_const_terms(cancer_form(m, dx), sub)
    mod <- lm(form_ok, data = sub)
    p <- get_resp_p(mod)
    all_results <- rbind(all_results, data.frame(metric = metric_label, level = dx, p_adj = p, stringsAsFactors = FALSE))
  }
}

# BH correction (overall and per-cancer separately)
all_results$p_adj_fdr <- p.adjust(all_results$p_adj, method = "BH")
write.csv(all_results, file.path(outdir, "adjusted_R_vs_NR_pvals.csv"), row.names = FALSE)
cat("\n== Adjusted Response_6 p-values (FDR within metric) ==\n")
print(all_results, digits = 4)

# ── 5b. Adjusted pairwise between-cancer comparisons (Tukey, from the same overall model)────────
cancer_pair_all <- data.frame()
for (m in metrics_log) {
  metric_label <- gsub("log_", "", m)
  # Between-cancer comparison: focus on Diagnosis -> host only (smoking is host; histology is tumour-intrinsic)
  mod_all <- lm(drop_const_terms(as.formula(paste0(m, " ~ Diagnosis + Age + sex + BMI_category + smoking")), df),
                data = df)
  # Manually build the Tukey contrast matrix to obtain contrast names
  lv <- levels(df$Diagnosis)
  cm <- combn(lv, 2)
  K <- t(apply(cm, 2, function(pp) {
    v <- rep(0, length(coef(mod_all)))
    i1 <- match(paste0("Diagnosis", pp[1]), names(coef(mod_all)))
    i2 <- match(paste0("Diagnosis", pp[2]), names(coef(mod_all)))
    if (!is.na(i1)) v[i1] <- -1
    if (!is.na(i2)) v[i2] <- 1
    v
  }))
  rownames(K) <- paste0(cm[1, ], " - ", cm[2, ])
  gl  <- multcomp::glht(mod_all, linfct = K)
  sm  <- summary(gl)
  pv  <- sm$test$pvalues
  nm  <- rownames(K)
  g1  <- sub(" - .*$", "", nm)
  g2  <- sub("^.* - ", "", nm)
  cancer_pair_all <- rbind(cancer_pair_all, data.frame(
    metric = metric_label, group1 = g1, group2 = g2, p_adj = pv, stringsAsFactors = FALSE
  ))
}
write.csv(cancer_pair_all, file.path(outdir, "adjusted_between_cancer_pvals.csv"), row.names = FALSE)
cat("\n== Adjusted between-cancer (Tukey) p-values ==\n")
print(cancer_pair_all, digits = 4)

# ── 6. Plot (same geoms and colours as the original Fig1bfg; only p-values replaced)──
# Prepare annotation text
make_label <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# Bottom annotation for each panel: per-diagnosis adjusted p
get_annotate <- function(metric_label) {
  rows <- all_results[all_results$metric == metric_label & all_results$level %in% c("CRC","GC","EC"), ]
  overall_row <- all_results[all_results$metric == metric_label & all_results$level == "Overall", ]
  list(
    per_dx = rows,
    overall_p = overall_row$p_adj[1]
  )
}

# Adjusted between-cancer comparisons: draw brackets manually (segment + text) to avoid ggpubr rendering issues on the log scale
make_pair_manual <- function(metric_label, comparisons, metric_col, logscale = FALSE) {
  rows <- cancer_pair_all[cancer_pair_all$metric == metric_label, ]
  out <- list()
  for (i in seq_along(comparisons)) {
    g1 <- comparisons[[i]][1]; g2 <- comparisons[[i]][2]
    hit <- rows[(rows$group1 == g1 & rows$group2 == g2) | (rows$group1 == g2 & rows$group2 == g1), ]
    p <- if (nrow(hit) > 0) hit$p_adj[1] else NA
    if (!is.na(p) && p < 0.05) {
      out[[length(out) + 1]] <- data.frame(
        group1 = g1, group2 = g2, p = p, idx = i,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) return(list())
  man <- do.call(rbind, out)
  xpos <- setNames(seq_along(levels(df$Diagnosis)), levels(df$Diagnosis))
  ymax <- max(df[[metric_col]], na.rm = TRUE)
  rng  <- diff(range(df[[metric_col]], na.rm = TRUE))
  # Bracket geometry: ybar = horizontal bar height, ytip = lower end of bracket leg, ylab = asterisk height
  # Use multiplicative offsets on the log axis (otherwise leg lengths are invisible on an order-of-magnitude scale and collapse to short lines)
  if (logscale) {
    man$ybar <- ymax * 10^(0.30 * man$idx)
    man$ytip <- man$ybar / 10^(0.25)
    man$ylab <- man$ybar * 10^(0.04)
  } else {
    man$ybar <- ymax + rng * 0.12 * man$idx
    man$ytip <- man$ybar - rng * 0.05
    man$ylab <- man$ybar + rng * 0.02
  }
  man$label <- sapply(man$p, make_label)
  man$x1 <- xpos[man$group1]
  man$x2 <- xpos[man$group2]
  man
}

# Final plot: add adjusted p annotations at the bottom
plot_metric <- function(metric_col, metric_label, title, ylab, logscale = FALSE, breakvals = NULL, comparisons = list(c("CRC", "GC"), c("CRC", "EC")), show_legend = FALSE) {
  a <- get_annotate(metric_label)
  man <- make_pair_manual(metric_label, comparisons, metric_col, logscale)
  
  p <- ggplot(df, aes(x = Diagnosis, y = .data[[metric_col]], color = Response_6)) +
    geom_half_violin(
      data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[1]),
      fill = NA, side = "l", show.legend = FALSE, linewidth = 0.2
    ) +
    geom_half_violin(
      data = ~ subset(.x, Response_6 == levels(factor(.x$Response_6))[2]),
      fill = NA, side = "r", show.legend = FALSE, linewidth = 0.2
    ) +
    geom_half_boxplot(
      data = filter(df, Response_6 == "NR"),
      aes(x = Diagnosis, y = .data[[metric_col]]),
      width = 0.25, side = "l", fill = "#E41A1C", coef = 0, alpha = 0.5,
      outlier.shape = NA, linewidth = 0.2
    ) +
    geom_half_boxplot(
      data = filter(df, Response_6 == "R"),
      aes(x = Diagnosis, y = .data[[metric_col]]),
      width = 0.25, side = "r", fill = "#377EB8", coef = 0, alpha = 0.5,
      outlier.shape = NA, linewidth = 0.2
    ) +
    theme_bw() +
    { if (length(man) > 0) {
        layers <- Reduce(c, lapply(seq_len(nrow(man)), function(i) {
          list(
            annotate("segment", x = man$x1[i], xend = man$x2[i],
                     y = man$ybar[i], yend = man$ybar[i],
                     colour = "black", linewidth = 0.4),
            annotate("segment", x = man$x1[i], xend = man$x1[i],
                     y = man$ybar[i], yend = man$ytip[i],
                     colour = "black", linewidth = 0.4),
            annotate("segment", x = man$x2[i], xend = man$x2[i],
                     y = man$ybar[i], yend = man$ytip[i],
                     colour = "black", linewidth = 0.4),
            annotate("text", x = (man$x1[i] + man$x2[i]) / 2,
                     y = man$ylab[i],
                     label = man$label[i], size = 4, colour = "black")
          )
        }))
        layers
      } else NULL } +
    theme(
      strip.background = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5),
      axis.text.x = element_text(size = 12)
    ) +
    ylab(ylab) +
    xlab(NULL) +
    scale_color_manual(values = c("#E41A1C", "#377EB8"), name = "Response") +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = if (show_legend) "right" else "none")
  
  # On the log10 scale, do not use scale_y_log10(limits=...): it filters out the bracket layer; use coord_cartesian instead
  if (logscale) p <- p + scale_y_log10() + coord_cartesian(ylim = c(1e-03, 1e+03))
  if (!is.null(breakvals)) p <- p + scale_y_continuous(breaks = breakvals)
  
  # Bottom: per-cancer adjusted p for R vs NR (replaces the original stat_compare_means Wilcoxon)
  # y position: bottom of the plot
  y_bot <- min(df[[metric_col]], na.rm = TRUE)
  if (logscale) y_bot <- 10^(log10(y_bot) - 0.3)
  else y_bot <- y_bot - diff(range(df[[metric_col]], na.rm = TRUE)) * 0.08
  
  annotate_df <- data.frame(
    Diagnosis = a$per_dx$level,
    y = rep(y_bot, 3),
    label = sapply(a$per_dx$p_adj, make_label),
    stringsAsFactors = FALSE
  )
  p <- p + geom_text(data = annotate_df, aes(x = Diagnosis, y = y, label = label),
                     inherit.aes = FALSE, size = 3.5)
  
  return(p)
}

# Keep y-axis units consistent with the original (raw scale); for cumulative abundance use *100 and log10
df$total_obs  <- df$`total observed`
df$oral_orig  <- df$`oral-origin`
df$cum_pct    <- df$`cumlative abundance` * 100   # *100 to match the original figure

p1 <- plot_metric("total_obs", "total", "Total", "SGB counts", breakvals = c(200, 400, 600),
                  comparisons = list(c("CRC", "GC"), c("CRC", "EC")))
p2 <- plot_metric("oral_orig", "oral", "Putative oral", "SGB counts",
                  comparisons = list(c("CRC", "GC"), c("GC", "EC"), c("CRC", "EC")))
p3 <- plot_metric("cum_pct",   "cum", "Putative oral", "Cumulative abundance (%)", logscale = TRUE,
                  comparisons = list(c("CRC", "GC"), c("GC", "EC")), show_legend = TRUE)

p_all <- p1 + p2 + p3

ggsave(file.path(outdir, "Fig1bfg_adjusted.pdf"), plot = p_all, height = 2.5, width = 8)
ggsave(file.path(outdir, "Fig1bfg_adjusted.png"), plot = p_all, height = 2.5, width = 8, dpi = 300)
cat("\n== Figure saved to", outdir, "==\n")
