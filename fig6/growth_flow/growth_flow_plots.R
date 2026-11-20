#!/usr/bin/env Rscript

## 4 panels from mouse_growth_flow_data.xlsx:
##   1) tumor growth curve (mean +/- SEM per time point, line+points)
##   2-4) flow cytometry: hollow bars (mean) + jitter points + SD error bars
## combined into a single-row figure. Tumor curve significance is based on a
## repeated-measures LMM (log10(volume + 1) ~ group * day + mouse random intercept).
## PD-1 + butyrate group excluded; group names aligned with alpha diversity naming.
## Flow panels: pairwise Wilcoxon (raw P values; no BH correction), only significant pairs labelled.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(lme4)
  library(emmeans)
  library(igraph)
})

root   <- "."
src_xlsx <- "mouse_growth_flow_data.xlsx"
outdir <- "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
xlsx_tmp <- tempfile(fileext = ".xlsx")
file.copy(src_xlsx, xlsx_tmp, overwrite = TRUE)

## ---- group naming aligned with alpha diversity (Shannon/PCoA) ----
group_order <- c("Control", "PD-1", "E. coli + PD-1",
                 "VP + PD-1", "VP + PD-1 + butyrate")
group_colors <- c(
  "Control"              = "#009E73",
  "PD-1"                 = "#CC79A7",
  "E. coli + PD-1"       = "#56B4E9",
  "VP + PD-1"            = "#D55E00",
  "VP + PD-1 + butyrate" = "#E69F00"
)
rename_grp <- c(CRTL = "Control", `PD-1` = "PD-1",
                E.coli = "E. coli + PD-1", `E.coli+PD-1` = "E. coli + PD-1",
                VP = "VP + PD-1", `VP+PD-1` = "VP + PD-1",
                `VP+butyrate` = "VP + PD-1 + butyrate",
                `PD-1+butyrate` = "PD-1 + butyrate")

star_p <- function(p) ifelse(p < 0.001, "***",
                    ifelse(p < 0.01, "**",
                    ifelse(p < 0.05, "*", "ns")))

## ============ read + tidy sheet 1: growth curve ============
g <- read_excel(xlsx_tmp, sheet = 1, col_names = FALSE)
gn <- as.character(unlist(g[1, ]))
gn[gn == "" | is.na(gn)] <- NA
for (i in 2:length(gn)) if (is.na(gn[i])) gn[i] <- gn[i - 1]
growth <- data.frame(
  day   = rep(as.numeric(unlist(g[-1, 1])), times = ncol(g) - 1),
  grp   = rep(gn[-1], each = nrow(g) - 1),
  mouse_col = rep(seq_len(ncol(g) - 1), each = nrow(g) - 1),
  value = as.numeric(as.matrix(g[-1, -1])),
  stringsAsFactors = FALSE
)
growth <- growth[!is.na(growth$value) & !is.na(growth$grp), ]
growth$grp <- unname(rename_grp[growth$grp])
growth <- growth[growth$grp %in% group_order, ]
growth$grp <- factor(growth$grp, levels = group_order)
growth$mouse_id <- paste0(growth$grp, "_M", growth$mouse_col)

## ---- growth-curve significance: repeated-measures LMM -------------------
## Use log10(volume + 1) because tumor-volume distributions are right-skewed;
## mouse_id accounts for repeated measurements within animals.
growth$log_volume <- log10(growth$value + 1)
growth_model <- lme4::lmer(
  log_volume ~ grp * day + (1 | mouse_id),
  data = growth,
  control = lme4::lmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun = 1e5))
)

growth_trends <- as.data.frame(summary(emmeans::emtrends(
  growth_model, ~ grp, var = "day", options = list()
)))
growth_pw_grid <- emmeans::contrast(
  emmeans::emtrends(growth_model, ~ grp, var = "day"),
  method = "pairwise", adjust = "BH"
)
growth_pw <- as.data.frame(summary(growth_pw_grid, infer = TRUE))

## Parse emmeans contrast labels. Parentheses are inserted around names with
## spaces/hyphens, but are not part of the actual group names.
contrast_clean <- gsub("\\(([^)]+)\\)", "\\1", growth_pw$contrast)
contrast_split <- strsplit(contrast_clean, " - ", fixed = TRUE)
growth_pw$group_1 <- unname(vapply(contrast_split, `[`, character(1), 1))
growth_pw$group_2 <- unname(vapply(contrast_split, `[`, character(1), 2))
growth_pw$comparison <- paste(growth_pw$group_1, "vs", growth_pw$group_2)
growth_pw$significance <- unname(star_p(growth_pw$p.value))

## Compact letters for legend labels: groups sharing a letter do not differ
## significantly (BH-adjusted pairwise slope test, P >= 0.05).
growth_levels <- levels(growth$grp)
slope_by_grp <- setNames(growth_trends$day.trend, growth_trends$grp)
nonsig <- matrix(FALSE, nrow = length(growth_levels),
                 ncol = length(growth_levels),
                 dimnames = list(growth_levels, growth_levels))
for (i in seq_len(nrow(growth_pw))) {
  if (!is.na(growth_pw$p.value[i]) && growth_pw$p.value[i] >= 0.05) {
    g1 <- growth_pw$group_1[i]; g2 <- growth_pw$group_2[i]
    if (g1 %in% growth_levels && g2 %in% growth_levels) {
      nonsig[g1, g2] <- TRUE; nonsig[g2, g1] <- TRUE
    }
  }
}
ns_cliques <- igraph::max_cliques(
  igraph::graph_from_adjacency_matrix(nonsig, mode = "undirected", diag = FALSE)
)
## Assign compact letters from high to low estimated growth rate.
clique_order <- order(vapply(ns_cliques, function(ii) max(slope_by_grp[ii]), numeric(1)),
                      decreasing = TRUE)
ns_cliques <- ns_cliques[clique_order]
letter_lists <- setNames(vector("list", length(growth_levels)), growth_levels)
for (k in seq_along(ns_cliques)) {
  members <- growth_levels[ns_cliques[[k]]]
  for (gg in members) letter_lists[[gg]] <- c(letter_lists[[gg]], letters[k])
}
growth_group_letters <- vapply(growth_levels, function(gg) {
  if (length(letter_lists[[gg]]) == 0) return(letters[which(growth_levels == gg)])
  paste(sort(unique(letter_lists[[gg]])), collapse = "")
}, character(1))

## Export statistical details.
growth_slope_out <- growth_trends
names(growth_slope_out)[names(growth_slope_out) == "grp"] <- "group"
growth_slope_out$letters <- unname(growth_group_letters[growth_slope_out$group])
write.csv(growth_slope_out,
          file.path(outdir, "growth_curve_group_slopes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(growth_pw[, c("comparison", "group_1", "group_2", "estimate", "SE",
                       "df", "t.ratio", "p.value", "significance")],
          file.path(outdir, "growth_curve_pairwise_slope_tests.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(summary(growth_model)),
           file.path(outdir, "growth_curve_model_summary.txt"))

cat("\nGrowth-curve slope comparisons (BH-adjusted):\n")
for (i in seq_len(nrow(growth_pw))) {
  cat(sprintf("  %s vs %s: p.adj = %.4g %s\n",
              growth_pw$group_1[i], growth_pw$group_2[i],
              growth_pw$p.value[i], growth_pw$significance[i]))
}
cat("Growth-curve legend letters:\n")
for (gg in growth_levels) {
  cat(sprintf("  %s: %s\n", gg, growth_group_letters[[gg]]))
}

## mean +/- SEM per time point (SEM = SD / sqrt(n))
gsum <- growth %>%
  group_by(grp, day) %>%
  summarise(m = mean(value), s = sd(value) / sqrt(n()), n = n(), .groups = "drop")

## ============ read + tidy sheets 2-4: flow cytometry ============
read_flow <- function(sheet) {
  f <- read_excel(xlsx_tmp, sheet = sheet, col_names = FALSE)
  gn <- as.character(unlist(f[2, ]))
  d <- data.frame(
    grp   = rep(gn, each = nrow(f) - 2),
    value = as.numeric(as.matrix(f[3:nrow(f), ])),
    stringsAsFactors = FALSE
  )
  d <- d[!is.na(d$value) & !is.na(d$grp), ]
  d$grp <- unname(rename_grp[d$grp])
  d <- d[d$grp %in% group_order, ]
  d$grp <- factor(d$grp, levels = group_order)
  d
}
mdsc <- read_flow(2)
gzm  <- read_flow(3)
cd8  <- read_flow(4)

## ============ panel A: growth curve (mean +/- SEM, line + points) ============
## ggplot does not allow significance brackets to be drawn inside a native
## legend. Therefore the native legend is replaced by a custom legend overlay:
## colored keys + labels on the left, and BH-adjusted slope-comparison brackets
## on the right. This keeps brackets "between legend groups" rather than on the
## growth curves themselves.
growth_bracket_pairs <- list(
  c("Control", "PD-1"),
  c("PD-1", "E. coli + PD-1"),
  c("E. coli + PD-1", "VP + PD-1"),
  c("PD-1", "VP + PD-1"),
  c("VP + PD-1", "VP + PD-1 + butyrate")
)
growth_bracket <- do.call(rbind, lapply(growth_bracket_pairs, function(pair) {
  idx <- which(
    (growth_pw$group_1 == pair[1] & growth_pw$group_2 == pair[2]) |
      (growth_pw$group_1 == pair[2] & growth_pw$group_2 == pair[1]))
  if (length(idx) == 1L && is.finite(growth_pw$p.value[idx]) &&
      growth_pw$p.value[idx] < 0.05) growth_pw[idx, , drop = FALSE] else NULL
}))
growth_bracket <- growth_bracket[order(match(growth_bracket$group_1, group_order) +
                                         match(growth_bracket$group_2, group_order)), ]

## custom legend coordinates: one row per group, top to bottom
legend_y <- setNames(rev(seq_along(group_order)), group_order)
legend_df <- data.frame(
  grp = factor(group_order, levels = group_order),
  x0 = 0,
  x1 = 0.45,
  xlab = 0.55,
  y = unname(legend_y[group_order]),
  lab = group_order,
  stringsAsFactors = FALSE
)

## assign bracket lanes greedily so overlapping vertical ranges use new lanes
growth_bracket$y1 <- unname(pmin(legend_y[growth_bracket$group_1],
                                 legend_y[growth_bracket$group_2]))
growth_bracket$y2 <- unname(pmax(legend_y[growth_bracket$group_1],
                                 legend_y[growth_bracket$group_2]))
growth_bracket$lane <- integer(nrow(growth_bracket))
if (nrow(growth_bracket) > 0) {
  used <- list()
  for (k in seq_len(nrow(growth_bracket))) {
    lane <- 1L
    while (TRUE) {
      hit <- vapply(used, function(z) {
        z$lane == lane &&
          growth_bracket$y1[k] <= z$y2 &&
          growth_bracket$y2[k] >= z$y1
      }, logical(1))
      if (!any(hit)) break else lane <- lane + 1L
    }
    growth_bracket$lane[k] <- lane
    used[[length(used) + 1L]] <- list(lane = lane,
                                      y1 = growth_bracket$y1[k],
                                      y2 = growth_bracket$y2[k])
  }
}

legend_bracket_start <- 4.2
legend_lane_width <- 0.9
legend_lane_offsets <- if (nrow(growth_bracket) > 0) {
  lane_ids <- sort(unique(growth_bracket$lane))
  lane_offset <- setNames(rep(0, length(lane_ids)), lane_ids)
  lane_offset[length(lane_offset)] <- -0.30
  lane_offset
} else numeric(0)
legend_xmax <- if (nrow(growth_bracket) > 0) {
  legend_bracket_start +
    (max(growth_bracket$lane) - 1) * legend_lane_width + 0.18
} else 1.35

legend_bracket_df <- do.call(rbind, lapply(seq_len(nrow(growth_bracket)), function(k) {
  xx <- legend_bracket_start +
    (growth_bracket$lane[k] - 1) * legend_lane_width +
    if (length(legend_lane_offsets)) {
      unname(legend_lane_offsets[as.character(growth_bracket$lane[k])])
    } else 0
  data.frame(
    x = xx,
    xend = xx,
    x0 = legend_bracket_start,
    y = growth_bracket$y1[k],
    yend = growth_bracket$y2[k],
    ymid = (growth_bracket$y1[k] + growth_bracket$y2[k]) / 2,
    xstar = xx + 0.20,
    label = paste(strsplit(growth_bracket$significance[k], "")[[1]],
                  collapse = "\n")
  )
}))

## main growth panel: no native legend
pA_base <- ggplot(gsum, aes(x = day, y = m, color = grp)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 0.5, linewidth = 1.0) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Tumor growth curve",
       x = "Days after inoculation", y = "Tumor volume (mm\u00B3)",
       caption = NULL) +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14,
                                  color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        panel.border = element_rect(color = "black", linewidth = 2.0),
        panel.grid = element_blank())

## custom transparent legend with simple vertical comparison lines
pA_legend <- ggplot() +
  geom_segment(data = legend_df,
               aes(x = x0, xend = x1, y = y, yend = y, color = grp),
               linewidth = 1.4, lineend = "round") +
  geom_point(data = legend_df,
             aes(x = (x0 + x1) / 2, y = y, color = grp),
             size = 3.0) +
  geom_text(data = legend_df,
            aes(x = xlab, y = y, label = lab),
            hjust = 0, vjust = 0.5, size = 4.2, color = "black",
            parse = FALSE, inherit.aes = FALSE) +
  scale_color_manual(values = group_colors, guide = "none") +
  coord_cartesian(xlim = c(0, legend_xmax), ylim = c(0.45, 5.55),
                  expand = FALSE) +
  theme_void(base_size = 10) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

if (nrow(legend_bracket_df) > 0) {
  pA_legend <- pA_legend +
    geom_segment(data = legend_bracket_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linewidth = 1.0, inherit.aes = FALSE) +
    geom_text(data = legend_bracket_df,
              aes(x = xstar, y = ymid, label = label),
              size = 4.8, color = "black", lineheight = 0.4,
              inherit.aes = FALSE)
}

## overlay custom legend where the ordinary ggplot legend would be placed
pA <- pA_base +
  annotation_custom(
    grob = ggplotGrob(pA_legend),
    xmin = min(growth$day),
    xmax = min(growth$day) + 0.62 * diff(range(growth$day)),
    ymin = max(gsum$m + gsum$s, na.rm = TRUE) * 0.55,
    ymax = max(gsum$m + gsum$s, na.rm = TRUE) * 0.98
  )

cat("\nTumor-growth legend lines (planned pairs, BH-adjusted):\n")
for (i in seq_len(nrow(growth_bracket))) {
  cat(sprintf("  %s vs %s: p.adj = %.4g %s\n",
              growth_bracket$group_1[i], growth_bracket$group_2[i],
              growth_bracket$p.value[i], growth_bracket$significance[i]))
}

## ============ panels B-D: hollow bars + jitter + SEM + significant comparisons ====
flow_panel <- function(df, title, ylab) {
  s <- df %>%
    group_by(grp) %>%
    summarise(m = mean(value), s = sd(value) / sqrt(n()), .groups = "drop")

  glev <- levels(df$grp)
  ## only the 5 requested comparisons
  comp <- list(
    c("Control", "PD-1"),
    c("PD-1", "E. coli + PD-1"),
    c("E. coli + PD-1", "VP + PD-1"),
    c("PD-1", "VP + PD-1"),
    c("VP + PD-1", "VP + PD-1 + butyrate")
  )
  comp <- comp[sapply(comp, function(cc) all(cc %in% glev))]
  pairs <- sapply(comp, identity)         # 2 x n matrix
  if (is.null(dim(pairs))) pairs <- matrix(pairs, nrow = 2)
  p_raw <- apply(pairs, 2, function(g2) {
    d2 <- df[df$grp %in% g2, ]
    p  <- suppressWarnings(tryCatch(
      wilcox.test(value ~ grp, data = d2, exact = FALSE)$p.value,
      error = function(e) NA_real_))
    p
  })
  sel <- which(p_raw < 0.05 & !is.na(p_raw))    # raw p-values, no correction
  sel <- sel[order((match(pairs[1, sel], glev) + match(pairs[2, sel], glev)) / 2)]

  ymax_data <- max(df$value)
  spacing <- 0.10 * ymax_data
  if (length(sel) > 0) {
    nb <- length(sel)
    ytop <- ymax_data + spacing * (nb + 1.4)     # headroom for all brackets
    y_i  <- ymax_data + spacing * seq_len(nb)
  } else {
    nb <- 0
    ytop <- ymax_data * 1.12
    y_i  <- numeric(0)
  }

  p <- ggplot() +
    geom_col(data = s, aes(x = grp, y = m, color = grp),
             fill = "white", linewidth = 1.4, width = 0.6) +
    geom_errorbar(data = s, aes(x = grp, ymin = m - s, ymax = m + s,
                                color = grp),
                  width = 0.16, linewidth = 1.0) +
    geom_jitter(data = df, aes(x = grp, y = value, color = grp),
                position = position_jitter(width = 0.12, seed = 42),
                size = 1.8, alpha = 0.6) +
    scale_color_manual(values = group_colors, guide = "none") +
    scale_y_continuous(limits = c(0, ytop), expand = c(0, 0)) +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14,
                                    color = "black"),
          axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size = 14, color = "black"),
          panel.border = element_rect(color = "black", linewidth = 2.0),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 14,
                                     color = "black"))

  if (nb > 0) {
    cat(sprintf("\n[%s] significant comparisons (raw Wilcoxon p < 0.05):\n", title))
    tip <- 0.025 * ymax_data
    for (k in seq_len(nb)) {
      i     <- sel[k]
      g1 <- pairs[1, i]; g2 <- pairs[2, i]
      x1 <- match(g1, glev); x2 <- match(g2, glev)
      y  <- y_i[k]
      ## simple straight comparison line, no bracket tips
      p <- p +
        annotate("segment", x = x1 + 0.10, xend = x2 - 0.10, y = y, yend = y,
                 colour = "black", linewidth = 1.0) +
        annotate("text", x = (x1 + x2) / 2, y = y + 0.30 * spacing,
                 label = star_p(p_raw[i]), size = 5.2, colour = "black")
      cat(sprintf("  %s vs %s : p = %.4f  %s\n",
                  g1, g2, p_raw[i], star_p(p_raw[i])))
    }
  }
  p
}

pB <- flow_panel(mdsc, "MDSC", "MDSC (%)")
pC <- flow_panel(gzm,  "IFN-\u03B3+ GranzymeB+ cells", "Cells (% of CD8)")
pD <- flow_panel(cd8,  "CD8+ T cells", "Cells (% of CD45)")

## ============ save individual panels ============
ggsave(file.path(outdir, "growth_curve.png"),    pA, width = 5.6, height = 3.9, dpi = 300)
ggsave(file.path(outdir, "growth_curve.pdf"),    pA, width = 5.6, height = 3.9, device = cairo_pdf)
ggsave(file.path(outdir, "flow_mdsc.png"),       pB, width = 4.8, height = 4.0, dpi = 300)
ggsave(file.path(outdir, "flow_ifn_gzmz.png"),   pC, width = 4.8, height = 4.0, dpi = 300)
ggsave(file.path(outdir, "flow_cd8.png"),        pD, width = 4.8, height = 4.0, dpi = 300)

## ============ flow-cytometry figure: the three flow panels in one row ======
## canvases: flow figure 9 x 3.8 in, focused SCFA figure 10.5 x 3.8 in
flow_row <- pB + pC + pD + plot_layout(nrow = 1)
ggsave(file.path(outdir, "flow_3panel.png"),
       flow_row, width = 9, height = 3.8, dpi = 300)
ggsave(file.path(outdir, "flow_3panel.pdf"),
       flow_row, width = 9, height = 3.8, device = cairo_pdf)

cat("saved:", file.path(outdir, "growth_curve.png"), "\n")
cat("saved:", file.path(outdir, "flow_mdsc.png"), "\n")
cat("saved:", file.path(outdir, "flow_ifn_gzmz.png"), "\n")
cat("saved:", file.path(outdir, "flow_cd8.png"), "\n")
cat("saved:", file.path(outdir, "flow_3panel.png"), "\n")
