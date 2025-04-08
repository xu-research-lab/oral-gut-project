# --- 1. Install and load necessary packages ---
# If you have not installed these packages, please run the following command first
# install.packages(c("vegan", "dplyr", "knitr", "ggplot2"))

# Set a random seed to ensure reproducibility of results
set.seed(123)

# Load R packages
library(vegan) # For using the adonis2 function for PERMANOVA analysis
library(dplyr) # For data processing and manipulation
library(knitr) # For creating well-formatted tables
library(ggplot2) # For data visualization

# --- 2. Load your data ---
# Note: If your meta.csv file contains Chinese characters and fails to read by default,
# the script will try to use 'GBK' encoding, which is a common way to handle Chinese encoding issues.
tryCatch(
    {
        meta_data <- read.csv("meta.csv", row.names = 1, header = TRUE)
    },
    error = function(e) {
        message("Failed to read meta.csv with default encoding, trying GBK encoding...")
        meta_data <- read.csv("meta.csv", row.names = 1, header = TRUE, encoding = "GBK")
    }
)

################################################################################
# --- First round of analysis: Based on Putative oral strain data ---
################################################################################

cat("============================================================\n")
cat("--- Starting first round of analysis: Based on oral strain data ---\n")
cat("============================================================\n\n")

# --- 2.1 Load oral strain distance matrices ---
bray_curtis_oral <- read.csv("big_cohort_profile_oral_strains_bray-curtis.tsv", sep = "\t", row.names = 1, header = TRUE)
jaccard_oral <- read.csv("big_cohort_profile_oral_strains_jaccard.tsv", sep = "\t", row.names = 1, header = TRUE)
unweighted_unifrac_oral <- read.csv("big_cohort_profile_oral_strains_unweighted-unifrac.tsv", sep = "\t", row.names = 1, header = TRUE)
weighted_unifrac_oral <- read.csv("big_cohort_profile_oral_strains_weighted-unifrac.tsv", sep = "\t", row.names = 1, header = TRUE)

# --- 3.1 Prepare a list of data for loop processing ---
dist_matrices_oral <- list(
    `Bray-Curtis` = bray_curtis_oral,
    `Jaccard` = jaccard_oral,
    `Unweighted UniFrac` = unweighted_unifrac_oral,
    `Weighted UniFrac` = weighted_unifrac_oral
)

# --- Part 1 Analysis: Adonis2 (PERMANOVA) analysis grouped by 'Diagnosis' ---

cat("--- Starting Part 1 analysis: Grouped by 'Diagnosis' ---\n\n")
cat("HOMD oral species \n")
diagnosis_results_oral <- list() # Initialize an empty list to store the results

# Loop through each distance matrix
for (dist_name in names(dist_matrices_oral)) {
    current_dist_matrix <- dist_matrices_oral[[dist_name]]

    # Data alignment: find common samples in metadata and distance matrix
    common_samples <- intersect(rownames(meta_data), rownames(current_dist_matrix))
    meta_subset <- meta_data[common_samples, , drop = FALSE]

    # Preprocessing: remove NA values in the grouping variable 'Diagnosis'
    if (any(is.na(meta_subset[["Diagnosis"]]))) {
        meta_subset <- meta_subset[!is.na(meta_subset[["Diagnosis"]]), , drop = FALSE]
        common_samples <- rownames(meta_subset)
    }

    # Reorganize the distance matrix based on the aligned samples
    dist_aligned <- as.dist(as.matrix(current_dist_matrix)[common_samples, common_samples])

    # Run adonis2 analysis
    adonis_result <- adonis2(dist_aligned ~ Diagnosis, data = meta_subset, permutations = 999)

    # Store key results (R² and P-value) in the list
    diagnosis_results_oral[[dist_name]] <- data.frame(
        Grouping_Variable = "Diagnosis",
        Distance_Metric = dist_name,
        R2 = adonis_result$R2[1],
        P_value = adonis_result$`Pr(>F)`[1]
    )
}

# Combine all results from the list into a data frame and format with kable
diagnosis_table_final_oral <- do.call(rbind, diagnosis_results_oral)
rownames(diagnosis_table_final_oral) <- NULL

print(
    kable(
        diagnosis_table_final_oral,
        format = "pipe",
        col.names = c("Grouping Variable", "Distance Metric", "R²", "P-value"),
        caption = "Table 1: Adonis2 Analysis (Oral Strains, Grouped by Diagnosis)",
        digits = 5
    )
)
cat("\n\n") # Add a newline to separate the two tables


# --- Part 2 Analysis: Adonis2 analysis grouped by 'Response_6' within each cancer type ---

cat("--- Starting Part 2 analysis: Grouped by 'Response_6' within each cancer type ---\n\n")
response_results_oral <- list() # Initialize an empty list to store the results

# Preprocessing: remove rows with NA in the Diagnosis column of the metadata to avoid errors in the following loop
valid_meta_data <- meta_data[!is.na(meta_data$Diagnosis), ]
diagnosis_levels <- unique(valid_meta_data$Diagnosis)

# Outer loop: iterate through each cancer type
for (cancer_type in diagnosis_levels) {
    # Filter the metadata for the current cancer type
    meta_cancer_subset <- valid_meta_data[valid_meta_data$Diagnosis == cancer_type, , drop = FALSE]

    # Inner loop: iterate through each distance matrix
    for (dist_name in names(dist_matrices_oral)) {
        current_dist_matrix <- dist_matrices_oral[[dist_name]]

        # Data alignment and preprocessing
        common_samples <- intersect(rownames(meta_cancer_subset), rownames(current_dist_matrix))
        meta_subset <- meta_cancer_subset[common_samples, , drop = FALSE]

        # Remove NA values from Response_6
        if (any(is.na(meta_subset[["Response_6"]]))) {
            cat(sprintf("NA values found in 'Response_6' for cancer type '%s', removed in '%s' analysis.\n", cancer_type, dist_name))
            meta_subset <- meta_subset[!is.na(meta_subset[["Response_6"]]), , drop = FALSE]
            common_samples <- rownames(meta_subset)
        }

        # Robustness check: ensure there are enough samples and at least 2 groups for analysis
        if (nrow(meta_subset) < 3 || length(unique(meta_subset$Response_6)) < 2) {
            cat(sprintf("After preparing data for cancer type '%s' and distance '%s', there are not enough samples or fewer than 2 groups. Skipping analysis.\n", cancer_type, dist_name))
            next # Skip the current iteration
        }

        dist_aligned <- as.dist(as.matrix(current_dist_matrix)[common_samples, common_samples])

        # Run adonis2
        adonis_result <- adonis2(dist_aligned ~ Response_6, data = meta_subset, permutations = 999)

        # Store results
        response_results_oral[[paste(cancer_type, dist_name)]] <- data.frame(
            Cancer_Type = cancer_type,
            Distance_Metric = dist_name,
            R2 = adonis_result$R2[1],
            P_value = adonis_result$`Pr(>F)`[1]
        )
    }
}

# Generate and print the second table
if (length(response_results_oral) > 0) {
    response_table_final_oral <- do.call(rbind, response_results_oral)
    rownames(response_table_final_oral) <- NULL

    print(
        kable(
            response_table_final_oral,
            format = "pipe",
            col.names = c("Cancer Type (Diagnosis)", "Distance Metric", "R²", "P-value"),
            caption = "Table 2: Adonis2 Analysis (Oral Strains, Grouped by Response_6 within each cancer type)",
            digits = 5
        )
    )
} else {
    print("No results were generated for the second part of the analysis. Please check the data.")
}
cat("\n\n")

# --- Part 3 Analysis: PCoA Visualization (based on Bray-Curtis distance) ---

# 1. Data alignment: ensure samples in the distance matrix and metadata are consistent
common_samples_pcoa <- intersect(rownames(meta_data), rownames(bray_curtis_oral))
bray_curtis_aligned <- as.matrix(bray_curtis_oral)[common_samples_pcoa, common_samples_pcoa]
meta_data_aligned <- meta_data[common_samples_pcoa, , drop = FALSE]

# 2. Perform PCoA analysis
pcoa_result <- cmdscale(as.dist(bray_curtis_aligned), k = 2, eig = TRUE)

# 3. Calculate and prepare variance explained labels for the axes
eigenvalues <- pcoa_result$eig
variance_explained <- eigenvalues / sum(eigenvalues[eigenvalues > 0])
pco1_var <- round(variance_explained[1] * 100, 1)
pco2_var <- round(variance_explained[2] * 100, 1)
xlab_text <- sprintf("PCo1 (%.1f%%)", pco1_var)
ylab_text <- sprintf("PCo2 (%.1f%%)", pco2_var)

# 4. Prepare data frame for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("V1", "V2")
pcoa_coords$SampleID <- rownames(pcoa_coords)
meta_data_aligned$SampleID <- rownames(meta_data_aligned)
shared_big <- merge(pcoa_coords, meta_data_aligned, by = "SampleID")

# 5. Calculate the centroid and standard deviation for each group to draw error bars
shared_big_subset <- subset(shared_big, Response_6 %in% c("R", "NR"))
pcoa_mean <- aggregate(shared_big_subset[, c("V1", "V2")], by = list(Diagnosis = shared_big_subset$Diagnosis, Response_6 = shared_big_subset$Response_6), FUN = mean)
pcoa_sd <- aggregate(shared_big_subset[, c("V1", "V2")], by = list(Diagnosis = shared_big_subset$Diagnosis, Response_6 = shared_big_subset$Response_6), FUN = sd)
colnames(pcoa_sd)[3:4] <- c("Pco1_sd", "Pco2_sd")
pcoa_plot_data <- merge(pcoa_mean, pcoa_sd, by = c("Diagnosis", "Response_6"))
colnames(pcoa_plot_data)[1:4] <- c("Diagnosis", "Response_6", "Pco1", "Pco2")

# 6. [Dynamic Annotation] Dynamically generate Adonis2 annotations from analysis results
adonis_plot_results <- response_table_final_oral[response_table_final_oral$Distance_Metric == "Bray-Curtis", ]
annotation_coords <- data.frame(Cancer_Type = c("CRC", "GC", "EC"), x = c(-0.35, -0.35, -0.35), y = c(-0.18, -0.22, -0.26))
adonis_annotations <- merge(adonis_plot_results, annotation_coords, by = "Cancer_Type")
names(adonis_annotations)[names(adonis_annotations) == "Cancer_Type"] <- "Diagnosis"

overall_diag_result <- diagnosis_table_final_oral[diagnosis_table_final_oral$Distance_Metric == "Bray-Curtis", ]
overall_diag_label <- sprintf("Diagnosis \nR2=%.1f%% P=%.3f", overall_diag_result$R2 * 100, overall_diag_result$P_value)

# 7. Plot using ggplot2
pcoa_plot_oral <- ggplot(pcoa_plot_data, aes(x = Pco1, y = Pco2, color = Diagnosis)) +
    geom_errorbarh(aes(xmin = Pco1 - Pco1_sd, xmax = Pco1 + Pco1_sd), height = 0.01, size = 0.3) +
    geom_errorbar(aes(ymin = Pco2 - Pco2_sd, ymax = Pco2 + Pco2_sd), width = 0.01, size = 0.3) +
    geom_point(size = 8) +
    geom_text(aes(label = Response_6), color = "black", size = 3.5) +
    geom_text(data = adonis_annotations, aes(x = x, y = y, label = paste0("R2=", round(R2 * 100, 1), "%", " P=", format(P_value, digits = 2)), color = Diagnosis), size = 3.5, show.legend = FALSE, hjust = 0) +
    geom_text(aes(x = -0.32, y = 0.3, label = overall_diag_label), inherit.aes = FALSE, size = 3.5, hjust = 0) +
    labs(x = xlab_text, y = ylab_text) +
    theme_bw() +
    scale_color_manual(name = "Diagnosis", breaks = c("CRC", "GC", "EC"), values = c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3")) +
    ggtitle("Putative oral") +
    theme(plot.title = element_text(hjust = 0.5))

print(pcoa_plot_oral)
ggsave(pcoa_plot_oral, filename = "fig1c.png", width = 4, height = 3)
ggsave(pcoa_plot_gut, filename = "fig1c.pdf", width = 3.5, height = 2.5)


################################################################################
# --- Second round of analysis: Based on Gut microbiome data ---
################################################################################

cat("============================================================\n")
cat("--- Starting second round of analysis: Based on gut microbiome data ---\n")
cat("============================================================\n\n")

# --- 2.2 Load gut microbiome distance matrices ---
bray_curtis_gut <- read.csv("big_cohort_profile_bray-curtis.tsv", sep = "\t", row.names = 1, header = TRUE)
jaccard_gut <- read.csv("big_cohort_profile_jaccard.tsv", sep = "\t", row.names = 1, header = TRUE)
unweighted_unifrac_gut <- read.csv("big_cohort_profile_unweighted-unifrac.tsv", sep = "\t", row.names = 1, header = TRUE)
weighted_unifrac_gut <- read.csv("big_cohort_profile_weighted-unifrac.tsv", sep = "\t", row.names = 1, header = TRUE)

# --- 3.2 Prepare a list of data for loop processing ---
dist_matrices_gut <- list(
    `Bray-Curtis` = bray_curtis_gut,
    `Jaccard` = jaccard_gut,
    `Unweighted UniFrac` = unweighted_unifrac_gut,
    `Weighted UniFrac` = weighted_unifrac_gut
)

# --- Part 1 Analysis: Adonis2 (PERMANOVA) analysis grouped by 'Diagnosis' ---

cat("--- Starting Part 1 analysis (Gut Microbiome): Grouped by 'Diagnosis' ---\n\n")
diagnosis_results_gut <- list()

for (dist_name in names(dist_matrices_gut)) {
    current_dist_matrix <- dist_matrices_gut[[dist_name]]
    common_samples <- intersect(rownames(meta_data), rownames(current_dist_matrix))
    meta_subset <- meta_data[common_samples, , drop = FALSE]
    if (any(is.na(meta_subset[["Diagnosis"]]))) {
        meta_subset <- meta_subset[!is.na(meta_subset[["Diagnosis"]]), , drop = FALSE]
        common_samples <- rownames(meta_subset)
    }
    dist_aligned <- as.dist(as.matrix(current_dist_matrix)[common_samples, common_samples])
    adonis_result <- adonis2(dist_aligned ~ Diagnosis, data = meta_subset, permutations = 999)
    diagnosis_results_gut[[dist_name]] <- data.frame(
        Grouping_Variable = "Diagnosis",
        Distance_Metric = dist_name,
        R2 = adonis_result$R2[1],
        P_value = adonis_result$`Pr(>F)`[1]
    )
}

diagnosis_table_final_gut <- do.call(rbind, diagnosis_results_gut)
rownames(diagnosis_table_final_gut) <- NULL
print(
    kable(
        diagnosis_table_final_gut,
        format = "pipe",
        col.names = c("Grouping Variable", "Distance Metric", "R²", "P-value"),
        caption = "Table 3: Adonis2 Analysis (Gut Microbiome, Grouped by Diagnosis)",
        digits = 5
    )
)
cat("\n\n")

# --- Part 2 Analysis: Adonis2 analysis grouped by 'Response_6' within each cancer type ---

cat("--- Starting Part 2 analysis (Gut Microbiome): Grouped by 'Response_6' within each cancer type ---\n\n")
response_results_gut <- list()

for (cancer_type in diagnosis_levels) {
    meta_cancer_subset <- valid_meta_data[valid_meta_data$Diagnosis == cancer_type, , drop = FALSE]
    for (dist_name in names(dist_matrices_gut)) {
        current_dist_matrix <- dist_matrices_gut[[dist_name]]
        common_samples <- intersect(rownames(meta_cancer_subset), rownames(current_dist_matrix))
        meta_subset <- meta_cancer_subset[common_samples, , drop = FALSE]
        if (any(is.na(meta_subset[["Response_6"]]))) {
            cat(sprintf("NA values found in 'Response_6' for cancer type '%s', removed in '%s' analysis.\n", cancer_type, dist_name))
            meta_subset <- meta_subset[!is.na(meta_subset[["Response_6"]]), , drop = FALSE]
            common_samples <- rownames(meta_subset)
        }
        if (nrow(meta_subset) < 3 || length(unique(meta_subset$Response_6)) < 2) {
            cat(sprintf("After preparing data for cancer type '%s' and distance '%s', there are not enough samples or fewer than 2 groups. Skipping analysis.\n", cancer_type, dist_name))
            next
        }
        dist_aligned <- as.dist(as.matrix(current_dist_matrix)[common_samples, common_samples])
        adonis_result <- adonis2(dist_aligned ~ Response_6, data = meta_subset, permutations = 999)
        response_results_gut[[paste(cancer_type, dist_name)]] <- data.frame(
            Cancer_Type = cancer_type,
            Distance_Metric = dist_name,
            R2 = adonis_result$R2[1],
            P_value = adonis_result$`Pr(>F)`[1]
        )
    }
}

if (length(response_results_gut) > 0) {
    response_table_final_gut <- do.call(rbind, response_results_gut)
    rownames(response_table_final_gut) <- NULL
    print(
        kable(
            response_table_final_gut,
            format = "pipe",
            col.names = c("Cancer Type (Diagnosis)", "Distance Metric", "R²", "P-value"),
            caption = "Table 4: Adonis2 Analysis (Gut Microbiome, Grouped by Response_6 within each cancer type)",
            digits = 5
        )
    )
} else {
    print("No results were generated for the second part of the second round of analysis. Please check the data.")
}
cat("\n\n")

# --- Part 3 Analysis: PCoA Visualization (based on Bray-Curtis distance) ---

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

adonis_plot_results_gut <- response_table_final_gut[response_table_final_gut$Distance_Metric == "Bray-Curtis", ]
annotation_coords_gut <- data.frame(Cancer_Type = c("CRC", "GC", "EC"), x = c(-0.2, -0.2, -0.2), y = c(-0.12, -0.16, -0.20))
adonis_annotations_gut <- merge(adonis_plot_results_gut, annotation_coords_gut, by = "Cancer_Type")
names(adonis_annotations_gut)[names(adonis_annotations_gut) == "Cancer_Type"] <- "Diagnosis"

overall_diag_result_gut <- diagnosis_table_final_gut[diagnosis_table_final_gut$Distance_Metric == "Bray-Curtis", ]
overall_diag_label_gut <- sprintf("Diagnosis \nR2=%.1f%% P=%.3f", overall_diag_result_gut$R2 * 100, overall_diag_result_gut$P_value)

pcoa_plot_gut <- ggplot(pcoa_plot_data_gut, aes(x = Pco1, y = Pco2, color = Diagnosis)) +
    geom_errorbarh(aes(xmin = Pco1 - Pco1_sd, xmax = Pco1 + Pco1_sd), height = 0.01, size = 0.3) +
    geom_errorbar(aes(ymin = Pco2 - Pco2_sd, ymax = Pco2 + Pco2_sd), width = 0.01, size = 0.3) +
    geom_point(size = 8) +
    geom_text(aes(label = Response_6), color = "black", size = 3.5) +
    geom_text(data = adonis_annotations_gut, aes(x = x, y = y, label = paste0("R2=", round(R2 * 100, 1), "%", " P=", format(P_value, digits = 2)), color = Diagnosis), size = 3.5, show.legend = FALSE, hjust = 0) +
    geom_text(aes(x = -0.2, y = 0.17, label = overall_diag_label_gut), inherit.aes = FALSE, size = 3.5, hjust = 0) +
    labs(x = xlab_text_gut, y = ylab_text_gut) +
    theme_bw() +
    scale_color_manual(name = "Diagnosis", breaks = c("CRC", "GC", "EC"), values = c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3")) +
    ggtitle("Total") +
    theme(plot.title = element_text(hjust = 0.5))

print(pcoa_plot_gut)
ggsave(pcoa_plot_gut, filename = "fig1d.png", width = 3.4, height = 2.7)
ggsave(pcoa_plot_gut, filename = "fig1d.pdf", width = 3.4, height = 2.7)
