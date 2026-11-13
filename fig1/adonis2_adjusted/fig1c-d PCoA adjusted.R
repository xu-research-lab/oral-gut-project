set.seed(123)
suppressMessages(library(vegan))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

# ---------- Read metadata (plotting uses only Diagnosis / Response_6, as in the original) ----------
tryCatch(
  {
      meta_data <- read.csv("meta.csv", row.names = 1, header = TRUE)
  },
  error = function(e) {
      message("Failed to read meta.csv with default encoding, trying GBK encoding...")
      meta_data <- read.csv("meta.csv", row.names = 1, header = TRUE, encoding = "GBK")
  }
)

# ---------- Adjusted adonis2 results (by="margin") ----------
adj <- read.csv("adonis2_adjusted/all_adonis2_adjusted.csv", check.names = FALSE)
adj_oral <- adj[adj$analysis %in% c("oral_bray | overall", "oral_bray | CRC", "oral_bray | EC", "oral_bray | GC"), ]
adj_gut  <- adj[adj$analysis %in% c("gut_bray | overall", "gut_bray | CRC", "gut_bray | EC", "gut_bray | GC"), ]

get_p <- function(df, analysis, term) {
  r <- df[df$analysis == analysis & df$term == term, ]
  if (nrow(r) == 0) return(NA_real_)
  pcol <- grep("^Pr", names(df), value = TRUE)[1]
  as.numeric(r[[pcol]][1])
}
get_r2 <- function(df, analysis, term) {
  r <- df[df$analysis == analysis & df$term == term, ]
  if (nrow(r) == 0) return(NA_real_)
  as.numeric(r[["R2"]][1])
}

# ---------- Oral strains ----------
bray_curtis_oral <- read.csv("big_cohort_profile_oral_strains_bray-curtis.tsv", sep = "\t", row.names = 1, header = TRUE)

common_samples_pcoa <- intersect(rownames(meta_data), rownames(bray_curtis_oral))
bray_curtis_aligned <- as.matrix(bray_curtis_oral)[common_samples_pcoa, common_samples_pcoa]
meta_data_aligned <- meta_data[common_samples_pcoa, , drop = FALSE]

pcoa_result <- cmdscale(as.dist(bray_curtis_aligned), k = 2, eig = TRUE)
eigenvalues <- pcoa_result$eig
variance_explained <- eigenvalues / sum(eigenvalues[eigenvalues > 0])
pco1_var <- round(variance_explained[1] * 100, 1)
pco2_var <- round(variance_explained[2] * 100, 1)
xlab_text <- sprintf("PCo1 (%.1f%%)", pco1_var)
ylab_text <- sprintf("PCo2 (%.1f%%)", pco2_var)

pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("V1", "V2")
pcoa_coords$SampleID <- rownames(pcoa_coords)
meta_data_aligned$SampleID <- rownames(meta_data_aligned)
shared_big <- merge(pcoa_coords, meta_data_aligned, by = "SampleID")

shared_big_subset <- subset(shared_big, Response_6 %in% c("R", "NR"))
pcoa_mean <- aggregate(shared_big_subset[, c("V1", "V2")], by = list(Diagnosis = shared_big_subset$Diagnosis, Response_6 = shared_big_subset$Response_6), FUN = mean)
pcoa_sd <- aggregate(shared_big_subset[, c("V1", "V2")], by = list(Diagnosis = shared_big_subset$Diagnosis, Response_6 = shared_big_subset$Response_6), FUN = sd)
colnames(pcoa_sd)[3:4] <- c("Pco1_sd", "Pco2_sd")
pcoa_plot_data <- merge(pcoa_mean, pcoa_sd, by = c("Diagnosis", "Response_6"))
colnames(pcoa_plot_data)[1:4] <- c("Diagnosis", "Response_6", "Pco1", "Pco2")

# Adjusted annotations (same structure as the original: one Response_6 P per cancer type and one overall Diagnosis P)
adonis_annotations <- data.frame(
  Cancer_Type = c("CRC", "GC", "EC"),
  R2 = sapply(c("oral_bray | CRC", "oral_bray | GC", "oral_bray | EC"), function(a) get_r2(adj_oral, a, "Response_6")),
  P_value = sapply(c("oral_bray | CRC", "oral_bray | GC", "oral_bray | EC"), function(a) get_p(adj_oral, a, "Response_6"))
)
adonis_annotations$text <- sprintf("R2=%.1f%% P=%.3f", adonis_annotations$R2 * 100, adonis_annotations$P_value)
adonis_annotations$x <- c(-0.35, -0.35, -0.35)
adonis_annotations$y <- c(-0.18, -0.22, -0.26)
names(adonis_annotations)[names(adonis_annotations) == "Cancer_Type"] <- "Diagnosis"

overall_diag_result <- data.frame(
  R2 = get_r2(adj_oral, "oral_bray | overall", "Diagnosis"),
  P_value = get_p(adj_oral, "oral_bray | overall", "Diagnosis")
)
overall_diag_label <- sprintf("Diagnosis \nR2=%.1f%% P=%.3f", overall_diag_result$R2 * 100, overall_diag_result$P_value)

pcoa_plot_oral <- ggplot(pcoa_plot_data, aes(x = Pco1, y = Pco2, color = Diagnosis)) +
    geom_errorbarh(aes(xmin = Pco1 - Pco1_sd, xmax = Pco1 + Pco1_sd), height = 0.01, size = 0.3) +
    geom_errorbar(aes(ymin = Pco2 - Pco2_sd, ymax = Pco2 + Pco2_sd), width = 0.01, size = 0.3) +
    geom_point(size = 8) +
    geom_text(aes(label = Response_6), color = "black", size = 3.5) +
    geom_text(data = adonis_annotations, aes(x = x, y = y, label = text, color = Diagnosis),
              size = 3.5, show.legend = FALSE, hjust = 0) +
    # top-left of the panel, identical relative placement in c and d; plain font
    annotate("text", x = -Inf, y = Inf, label = overall_diag_label,
             hjust = -0.08, vjust = 1.15, size = 3.5) +
    labs(x = xlab_text, y = ylab_text) +
    theme_bw() +
    scale_color_manual(name = "Diagnosis", breaks = c("CRC", "GC", "EC"), values = c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3")) +
    ggtitle("Putative oral") +
    theme(plot.title = element_text(hjust = 0.5))

tryCatch(ggsave(pcoa_plot_oral, filename = "adonis2_adjusted/fig1c_putative_oral_adjusted.pdf", width = 3.5, height = 2.5),
         error = function(e) { cat("fig1c pdf locked; writing _v2\n")
           ggsave(pcoa_plot_oral, filename = "adonis2_adjusted/fig1c_putative_oral_adjusted_v2.pdf", width = 3.5, height = 2.5) })

# ---------- Gut microbiome ----------
bray_curtis_gut <- read.csv("big_cohort_profile_bray-curtis.tsv", sep = "\t", row.names = 1, header = TRUE)

common_samples_pcoa_gut <- intersect(rownames(meta_data), rownames(bray_curtis_gut))
bray_curtis_aligned_gut <- as.matrix(bray_curtis_gut)[common_samples_pcoa_gut, common_samples_pcoa_gut]
meta_data_aligned_gut <- meta_data[common_samples_pcoa_gut, , drop = FALSE]

pcoa_result_gut <- cmdscale(as.dist(bray_curtis_aligned_gut), k = 2, eig = TRUE)
eigenvalues_gut <- pcoa_result_gut$eig
variance_explained_gut <- eigenvalues_gut / sum(eigenvalues_gut[eigenvalues_gut > 0])
pco1_var_gut <- round(variance_explained_gut[1] * 100, 1)
pco2_var_gut <- round(variance_explained_gut[2] * 100, 1)
xlab_text_gut <- sprintf("PCo1 (%.1f%%)", pco1_var_gut)
ylab_text_gut <- sprintf("PCo2 (%.1f%%)", pco2_var_gut)

pcoa_coords_gut <- as.data.frame(pcoa_result_gut$points)
colnames(pcoa_coords_gut) <- c("V1", "V2")
pcoa_coords_gut$SampleID <- rownames(pcoa_coords_gut)
meta_data_aligned_gut$SampleID <- rownames(meta_data_aligned_gut)
shared_big_gut <- merge(pcoa_coords_gut, meta_data_aligned_gut, by = "SampleID")

shared_big_subset_gut <- subset(shared_big_gut, Response_6 %in% c("R", "NR"))
pcoa_mean_gut <- aggregate(shared_big_subset_gut[, c("V1", "V2")], by = list(Diagnosis = shared_big_subset_gut$Diagnosis, Response_6 = shared_big_subset_gut$Response_6), FUN = mean)
pcoa_sd_gut <- aggregate(shared_big_subset_gut[, c("V1", "V2")], by = list(Diagnosis = shared_big_subset_gut$Diagnosis, Response_6 = shared_big_subset_gut$Response_6), FUN = sd)
colnames(pcoa_sd_gut)[3:4] <- c("Pco1_sd", "Pco2_sd")
pcoa_plot_data_gut <- merge(pcoa_mean_gut, pcoa_sd_gut, by = c("Diagnosis", "Response_6"))
colnames(pcoa_plot_data_gut)[1:4] <- c("Diagnosis", "Response_6", "Pco1", "Pco2")

adonis_annotations_gut <- data.frame(
  Cancer_Type = c("CRC", "GC", "EC"),
  R2 = sapply(c("gut_bray | CRC", "gut_bray | GC", "gut_bray | EC"), function(a) get_r2(adj_gut, a, "Response_6")),
  P_value = sapply(c("gut_bray | CRC", "gut_bray | GC", "gut_bray | EC"), function(a) get_p(adj_gut, a, "Response_6"))
)
adonis_annotations_gut$text <- sprintf("R2=%.1f%% P=%.3f", adonis_annotations_gut$R2 * 100, adonis_annotations_gut$P_value)
adonis_annotations_gut$x <- c(-0.2, -0.2, -0.2)
adonis_annotations_gut$y <- c(-0.12, -0.16, -0.20)
names(adonis_annotations_gut)[names(adonis_annotations_gut) == "Cancer_Type"] <- "Diagnosis"

overall_diag_result_gut <- data.frame(
  R2 = get_r2(adj_gut, "gut_bray | overall", "Diagnosis"),
  P_value = get_p(adj_gut, "gut_bray | overall", "Diagnosis")
)
overall_diag_label_gut <- sprintf("Diagnosis \nR2=%.1f%% P=%.3f", overall_diag_result_gut$R2 * 100, overall_diag_result_gut$P_value)

pcoa_plot_gut <- ggplot(pcoa_plot_data_gut, aes(x = Pco1, y = Pco2, color = Diagnosis)) +
    geom_errorbarh(aes(xmin = Pco1 - Pco1_sd, xmax = Pco1 + Pco1_sd), height = 0.01, size = 0.3) +
    geom_errorbar(aes(ymin = Pco2 - Pco2_sd, ymax = Pco2 + Pco2_sd), width = 0.01, size = 0.3) +
    geom_point(size = 8) +
    geom_text(aes(label = Response_6), color = "black", size = 3.5) +
    geom_text(data = adonis_annotations_gut, aes(x = x, y = y, label = text, color = Diagnosis),
              size = 3.5, show.legend = FALSE, hjust = 0) +
    annotate("text", x = -Inf, y = Inf, label = overall_diag_label_gut,
             hjust = -0.08, vjust = 1.15, size = 3.5) +
    labs(x = xlab_text_gut, y = ylab_text_gut) +
    theme_bw() +
    scale_color_manual(name = "Diagnosis", breaks = c("CRC", "GC", "EC"), values = c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    ggtitle("Total") +
    theme(plot.title = element_text(hjust = 0.5))

tryCatch(ggsave(pcoa_plot_gut, filename = "adonis2_adjusted/fig1d_gut_adjusted.pdf", width = 3.5, height = 2.5),
         error = function(e) { cat("fig1d pdf locked; writing _v2\n")
           ggsave(pcoa_plot_gut, filename = "adonis2_adjusted/fig1d_gut_adjusted_v2.pdf", width = 3.5, height = 2.5) })