library(dplyr)
library(ggplot2)
library(forcats)
library(tidyr)
library(ggpubr) 
library(scales) 
library(ggbeeswarm) 


load("Extended_Data_Fig1b.RData")
compare1 <- compare1[!is.na(compare1$id), ]
taxonomy_annotation$SGB <- gsub("t__", "", taxonomy_annotation$SGB)


oral_gut_abundance <- oral_gut_abundance[oral_gut_abundance$oral > 0, ]
oral_gut_abundance$shared_group <- ifelse(oral_gut_abundance$gut > 0 & oral_gut_abundance$oral > 0, "Present in both oral and gut samples", "Present only in oral samples")


oral_gut_abundance$species <- taxonomy_annotation$species[match(oral_gut_abundance$SGB, taxonomy_annotation$SGB)]
oral_gut_abundance$species <- gsub("s__", "", oral_gut_abundance$species)
oral_gut_abundance$species[oral_gut_abundance$species == "Candidatus_Nanosynsacchari_sp_TM7_ANC_38_39_G1_1"] <- "Candidatus_Nanosynsacchari"
oral_gut_abundance$species[oral_gut_abundance$species == "Candidatus_Saccharibacteria_unclassified_SGB19850"] <- "Candidatus_Saccharibacteria"
oral_gut_abundance$new_species <- paste(oral_gut_abundance$species, oral_gut_abundance$SGB, sep = " ")
oral_gut_abundance <- oral_gut_abundance %>%
    filter(new_species %in% compare1$id)
oral_gut_abundance$new_species <- gsub("_unclassified", "", oral_gut_abundance$new_species)
oral_gut_abundance$new_species <- gsub("_group", "", oral_gut_abundance$new_species)

oral_gut_abundance$Diagnosis[oral_gut_abundance$Diagnosis == "Untreated"] <- "RA"
oral_gut_abundance <- oral_gut_abundance %>%
    filter(Diagnosis %in% c("CRC", "EC", "GC", "HC", "RA"))


create_oral_abundance_dotplot <- function(df, group_filter, facet_by = NULL) {

    plot_data <- df %>%
        filter(group == group_filter)
    temp_species <- plot_data$new_species[plot_data$shared_group == "Present in both oral and gut samples"]
    plot_data <- plot_data %>%
        filter(new_species %in% temp_species) 

    species_order <- plot_data %>%
        filter(shared_group == "Present in both oral and gut samples") %>%
        group_by(new_species) %>%
        summarise(median_oral = median(oral, na.rm = TRUE)) %>%
        arrange(desc(median_oral)) %>%
        pull(new_species)
    plot_data$new_species <- factor(plot_data$new_species, levels = rev(species_order))


    plot_data$shared_group <- factor(plot_data$shared_group, levels = c("Present in both oral and gut samples", "Present only in oral samples"))


    median_data <- plot_data %>%
        group_by(new_species, shared_group) %>%
        summarise(median_val = median(oral, na.rm = TRUE), .groups = "drop")


    stat_test_groups <- if (!is.null(facet_by)) c("new_species", facet_by) else "new_species"
    stat_data <- plot_data %>%
        group_by(across(all_of(stat_test_groups))) %>%
        summarise(p_value = tryCatch(wilcox.test(oral ~ shared_group)$p.value, error = function(e) NA), .groups = "drop") %>%
        filter(!is.na(p_value))

    if (!is.null(facet_by)) {
        stat_data <- stat_data %>%
            group_by(across(all_of(facet_by))) %>%
            mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
            ungroup()
    } else {
        stat_data <- stat_data %>%
            mutate(p_adj = p.adjust(p_value, method = "BH"))
    }

    signif_data <- stat_data %>%
        filter(p_adj < 0.05) %>%
        mutate(p.signif = case_when(
            p_adj < 0.001 ~ "***",
            p_adj < 0.01 ~ "**",
            p_adj < 0.05 ~ "*",
            TRUE ~ ""
        ))


    total_species_count <- n_distinct(plot_data$new_species)
    significant_species_count <- nrow(signif_data)
    subtitle_text <- paste(
        "Total species shown:", total_species_count,
        "| Significant (p_adj < 0.05):", significant_species_count
    )

    min_x_val <- min(plot_data$oral[plot_data$oral > 0], na.rm = TRUE)
    max_x_val <- max(plot_data$oral, na.rm = TRUE)

    x_breaks <- if (is.finite(min_x_val) && is.finite(max_x_val) && min_x_val > 0) {
        10^seq(floor(log10(min_x_val)), ceiling(log10(max_x_val)), by = 1)
    } else {
        10^seq(-4, 2, by = 1)
    }

    x_star_pos <- min(x_breaks) * 0.2

    if (nrow(signif_data) > 0) {
        signif_data$x.position <- x_star_pos
    }



    p <- ggplot(plot_data, aes(x = new_species, y = oral, color = shared_group)) +
        geom_jitter(
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8, jitter.height = 0),
            size = 0.5,
            shape = 16,
            alpha = 0.8
        ) +
        stat_summary(
            aes(group = shared_group),
            fun = median,
            geom = "line",
            position = position_dodge(width = 0.8),
            linewidth = 0.4
        ) +
        stat_summary(
            aes(fill = shared_group),
            fun = median,
            geom = "point",
            shape = 23,
            size = 0.8,
            stroke = 0,
            position = position_dodge(width = 0.8),
            show.legend = FALSE
        ) +
        scale_y_log10(
            breaks = x_breaks,
            labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        {
            if (nrow(signif_data) > 0) {
                geom_text(
                    data = signif_data,
                    aes(x = new_species, y = x.position, label = p.signif),
                    inherit.aes = FALSE, size = 1.5, vjust = 0.5
                )
            }
        } +
        scale_color_manual(
            name = "Type",
            values = c("Present in both oral and gut samples" = "#6baed6", "Present only in oral samples" = "#fdbe85") 
        ) +
        scale_fill_manual(
            name = "Type",
            values = c("Present in both oral and gut samples" = "#08519c", "Present only in oral samples" = "#a63603") 
        ) +
        labs(
            title = paste("Oral Abundance Comparison in Group:", group_filter),
            subtitle = subtitle_text,
            x = NULL,
            y = "Abundance (%)"
        ) +
        coord_flip(clip = "off") +
        guides(color = guide_legend(nrow = 2, override.aes = list(size = 2))) +

        theme_bw() +
        theme(
            plot.title = element_text(hjust = 1),
            plot.subtitle = element_text(hjust = 1),
            axis.ticks.x = element_line(),
            legend.position = "top",
            legend.spacing.y = unit(0, "cm"), 
            legend.margin = margin(t = 0, b = 0), 
        )


    if (!is.null(facet_by)) {
        p <- p + facet_grid(reformulate(termlabels = ".", response = facet_by), scales = "free_x", space = "free_x")
    }

    return(p)
}



plot_non_pma_abundance <- create_oral_abundance_dotplot(oral_gut_abundance, "non-PMA")
if (!is.null(plot_non_pma_abundance)) {
    ggsave("Extended_Data_Fig1_b.png", plot = plot_non_pma_abundance, width = 7, height = 11, dpi = 300)
    ggsave("Extended_Data_Fig1_b.pdf", plot = plot_non_pma_abundance, width = 7, height = 13)
}
