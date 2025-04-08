library(vegan) 
library(dplyr) 
library(ggplot2) 
library(knitr) 


set.seed(123)
data_types <- c("non_PMA", "PMA")
sample_groups <- c("oral", "gut")
dist_metrics_names <- c(
    "Bray-Curtis" = "bray-curtis",
    "Jaccard" = "jaccard",
    "Unweighted UniFrac" = "unweighted-unifrac",
    "Weighted UniFrac" = "weighted-unifrac"
)


meta_data_full <- read.csv("fig4a.csv", sep = ",", row.names = 1, header = TRUE)


meta_data_full <- meta_data_full %>%
    filter(Diagnosis %in% c("CRC", "GC", "EC"))


meta_oral <- meta_data_full
meta_oral$Group <- "oral"

meta_gut <- meta_data_full
rownames(meta_gut) <- gsub("Saliva", "Feces", rownames(meta_gut))
meta_gut$Group <- "gut"

meta_data <- rbind(meta_oral, meta_gut)

all_results <- data.frame()


for (dtype in data_types) {
    dist_matrices <- list()
    for (i in 1:length(dist_metrics_names)) {
        dist_name <- names(dist_metrics_names)[i]
        file_suffix <- dist_metrics_names[i]
        file_path <- sprintf("total_abundance_profile_%s_%s.tsv", dtype, file_suffix)

        if (file.exists(file_path)) {
            dist_matrices[[dist_name]] <- read.csv(file_path, sep = "\t", row.names = 1, header = TRUE)
        } else {
            warning(sprintf("Warning: File %s does not exist, skipping related analysis.", file_path))
        }
    }

    if (length(dist_matrices) == 0) {
        cat(sprintf("No distance matrix files found for type %s, skipping this data type.\n\n", dtype))
        next
    }
    for (sgroup in sample_groups) {
        meta_group_subset <- meta_data %>% filter(Group == sgroup)
        for (dist_name in names(dist_matrices)) {
            current_dist_matrix <- dist_matrices[[dist_name]]
            common_samples <- intersect(rownames(meta_group_subset), rownames(current_dist_matrix))
            meta_subset <- meta_group_subset[common_samples, , drop = FALSE]
            dist_aligned <- as.dist(as.matrix(current_dist_matrix)[common_samples, common_samples])
            adonis_result <- adonis2(dist_aligned ~ Diagnosis, data = meta_subset, permutations = 999)
            result_row <- data.frame(
                Data_Type = dtype,
                Sample_Group = sgroup,
                Analysis = "Diagnosis",
                Cancer_Type = "All",
                Distance_Metric = dist_name,
                R2 = adonis_result$R2[1],
                P_value = adonis_result$`Pr(>F)`[1]
            )
            all_results <- rbind(all_results, result_row)
        }
        diagnosis_levels <- unique(meta_group_subset$Diagnosis)

        for (cancer_type in diagnosis_levels) {
            meta_cancer_subset <- meta_group_subset %>% filter(Diagnosis == cancer_type)

            for (dist_name in names(dist_matrices)) {
                current_dist_matrix <- dist_matrices[[dist_name]]
                common_samples <- intersect(rownames(meta_cancer_subset), rownames(current_dist_matrix))
                meta_subset <- meta_cancer_subset[common_samples, , drop = FALSE]
                if (any(is.na(meta_subset[["Response_6"]]))) {
                    meta_subset <- meta_subset[!is.na(meta_subset[["Response_6"]]), , drop = FALSE]
                    common_samples <- rownames(meta_subset)
                }
                if (nrow(meta_subset) < 3 || length(unique(meta_subset$Response_6)) < 2) {
                    next
                }

                dist_aligned <- as.dist(as.matrix(current_dist_matrix)[common_samples, common_samples])
                adonis_result <- adonis2(dist_aligned ~ Response_6, data = meta_subset, permutations = 999)
                result_row <- data.frame(
                    Data_Type = dtype,
                    Sample_Group = sgroup,
                    Analysis = "Response_6",
                    Cancer_Type = cancer_type,
                    Distance_Metric = dist_name,
                    R2 = adonis_result$R2[1],
                    P_value = adonis_result$`Pr(>F)`[1]
                )
                all_results <- rbind(all_results, result_row)
            }
        }
    }
}

if (nrow(all_results) > 0) {
    write.csv(all_results, "adonis2_comprehensive_results.csv", row.names = FALSE, quote = TRUE, fileEncoding = "UTF-8")

    display_results <- all_results %>%
        mutate(P_value = case_when(
            P_value < 0.001 ~ "< 0.001",
            TRUE ~ as.character(round(P_value, 3))
        ))

    print(
        kable(
            display_results,
            format = "pipe",
            col.names = c("Data Type", "Sample Group", "Analysis Variable", "Cancer Type", "Distance Metric", "R²", "P-value"),
            caption = "Adonis2 Comprehensive Analysis Results",
            digits = 4
        )
    )
} else {
    print("Failed to generate any analysis results. Please check your file paths and input data.")
}

cat("\n--- Starting to plot combined PCoA graph ---\n")

pcoa_plot_data_list <- list()
pcoa_annotation_list <- list()

for (dtype in data_types[1]) {
    for (sgroup in sample_groups) {
        cat(sprintf("Preparing data for PCoA plot: %s - %s\n", dtype, sgroup))

        bray_file <- sprintf("total_abundance_profile_%s_bray-curtis.tsv", dtype)
        if (!file.exists(bray_file)) {
            warning(sprintf("File %s not found, skipping PCoA plotting.", bray_file))
            next
        }
        bray_curtis_matrix <- read.csv(bray_file, sep = "\t", row.names = 1, header = TRUE)

        meta_group_subset <- meta_data %>% filter(Group == sgroup)
        common_samples <- intersect(rownames(meta_group_subset), rownames(bray_curtis_matrix))
        bray_curtis_aligned <- as.matrix(bray_curtis_matrix)[common_samples, common_samples]
        meta_aligned <- meta_group_subset[common_samples, , drop = FALSE]

        pcoa_result <- cmdscale(as.dist(bray_curtis_aligned), k = 2, eig = TRUE)
        pcoa_coords <- as.data.frame(pcoa_result$points)
        colnames(pcoa_coords) <- c("Pco1", "Pco2")

        pcoa_eig <- round(pcoa_result$eig / sum(pcoa_result$eig) * 100, 1)

        plot_df_full <- cbind(pcoa_coords, meta_aligned[rownames(pcoa_coords), ])
        plot_df_subset <- subset(plot_df_full, Response_6 %in% c("R", "NR"))

        pcoa_mean <- aggregate(plot_df_subset[, c("Pco1", "Pco2")], by = list(Diagnosis = plot_df_subset$Diagnosis, Response_6 = plot_df_subset$Response_6), FUN = mean)
        pcoa_sd <- aggregate(plot_df_subset[, c("Pco1", "Pco2")], by = list(Diagnosis = plot_df_subset$Diagnosis, Response_6 = plot_df_subset$Response_6), FUN = sd)
        colnames(pcoa_sd)[3:4] <- c("Pco1_sd", "Pco2_sd")

        pcoa_plot_data <- merge(pcoa_mean, pcoa_sd, by = c("Diagnosis", "Response_6"))
        pcoa_plot_data$Data_Type <- dtype
        pcoa_plot_data$Sample_Group <- sgroup

        pcoa_plot_data_list[[paste(dtype, sgroup)]] <- pcoa_plot_data

        adonis_response_results <- all_results %>%
            filter(Data_Type == dtype, Sample_Group == sgroup, Distance_Metric == "Bray-Curtis", Analysis == "Response_6")

        adonis_diag_result <- all_results %>%
            filter(Data_Type == dtype, Sample_Group == sgroup, Distance_Metric == "Bray-Curtis", Analysis == "Diagnosis") %>%
            slice(1) 

        adonis_response_results$label <- paste0("R2=", round(adonis_response_results$R2 * 100, 1), "%, P=", format(adonis_response_results$P_value, digits = 2))

        adonis_response_results <- adonis_response_results %>%
            mutate(
                x = -0.4, 
                y = case_when(
                    Cancer_Type == "CRC" ~ -0.24,
                    Cancer_Type == "GC" ~ -0.28,
                    Cancer_Type == "EC" ~ -0.32,
                    TRUE ~ -0.34 
                ),
                Data_Type = dtype,
                Sample_Group = sgroup
            )

        diag_label_df <- data.frame(
            label = sprintf("Diagnosis \nR2=%.1f%%, P=%.3f", adonis_diag_result$R2 * 100, adonis_diag_result$P_value),
            x = -0.4, 
            y = 0.22, 
            Data_Type = dtype,
            Sample_Group = sgroup
        )

        pcoa_annotation_list[[paste(dtype, sgroup)]] <- list(
            response_ann = adonis_response_results,
            diag_ann = diag_label_df,
            axis_labels = c(paste0("PCo1 [", pcoa_eig[1], "%]"), paste0("PCo2 [", pcoa_eig[2], "%]")),
            Data_Type = dtype,
            Sample_Group = sgroup
        )
    }
}


if (length(pcoa_plot_data_list) > 0) {
    combined_plot_data <- do.call(rbind, pcoa_plot_data_list)

    combined_response_annotations <- do.call(rbind, lapply(pcoa_annotation_list, `[[`, "response_ann"))

    combined_diag_annotations <- do.call(rbind, lapply(pcoa_annotation_list, `[[`, "diag_ann"))

    pcoa_combined_plot <- ggplot(combined_plot_data, aes(x = Pco1, y = Pco2, color = Diagnosis)) +
        facet_grid(~Sample_Group) +
        geom_errorbarh(aes(xmin = Pco1 - Pco1_sd, xmax = Pco1 + Pco1_sd), height = 0.01, linewidth = 0.3) +
        geom_errorbar(aes(ymin = Pco2 - Pco2_sd, ymax = Pco2 + Pco2_sd), width = 0.01, linewidth = 0.3) +
        geom_point(size = 8) +
        geom_text(aes(label = Response_6), color = "black", size = 3.5) +

        geom_text(
            data = combined_response_annotations,
            aes(x = x, y = y, label = paste0(label), color = Cancer_Type), 
            size = 3.5, show.legend = FALSE, hjust = 0
        ) +
        geom_text(
            data = combined_diag_annotations,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 3.5, hjust = 0
        ) +
        theme_bw() +
        scale_color_manual(
            name = "Diagnosis",
            breaks = c("CRC", "GC", "EC"),
            values = c("CRC" = "#66c2a5", "GC" = "#fc8d62", "EC" = "#e78ac3")
        ) +
        labs(x = "PCo1", y = "PCo2") +
        theme(
            strip.background = element_blank()
        )

    print(pcoa_combined_plot)

    ggsave("fig4_a.png", plot = pcoa_combined_plot, width = 6.4, height = 3.5, dpi = 300)
    ggsave("fig4_a.pdf", plot = pcoa_combined_plot, width = 6, height = 2.5)
}
