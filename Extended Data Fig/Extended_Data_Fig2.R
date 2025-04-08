library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)


load("Extended_Data_Fig2")

original_diagnosis_levels <- c("CRC", "GC", "EC", "RA untreated", "HC")

final_diagnosis_levels <- c("CRC", "GC", "EC", "RA", "HC")


diagnosis_colors <- c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb", "#a6d854")
names(diagnosis_colors) <- final_diagnosis_levels
plot_data_full <- trans2 %>%
    filter(treatment.y == "non-PMA") %>%
    filter(Diagnosis.x %in% original_diagnosis_levels) %>%
    mutate(Diagnosis.x = recode(Diagnosis.x, "RA untreated" = "RA")) %>%
    mutate(Diagnosis.x = factor(Diagnosis.x, levels = final_diagnosis_levels)) %>%
    filter(!is.na(centered.min) & !is.na(gut) & !is.na(species) & !is.na(trans))


plot_data_full$species <- gsub("s__", "", plot_data_full$species)
plot_data_full$species <- gsub("_SGB.*", "", plot_data_full$species)
plot_data_full$species <- gsub("_unclassified", "", plot_data_full$species)


plot_data_full$facet_label <- paste(plot_data_full$species, plot_data_full$SGB, sep = "\n")
plot_data_full$facet_label


plot_data_filtered <- plot_data_full %>%
    group_by(facet_label) %>%
    add_tally() %>%
    filter(n >= 6) %>%
    ungroup() %>%
    mutate(trend_group = ifelse(Diagnosis.x == "HC", "HC", "Cancer"))


count_order <- plot_data_filtered %>%
    distinct(facet_label, n) %>%
    arrange(desc(n)) %>%
    pull(facet_label)

plot_data_filtered$facet_label <- factor(plot_data_filtered$facet_label, levels = count_order)


p_Extended_Data_Fig2 <- ggplot(plot_data_filtered, aes(x = centered.min, y = gut)) +
    geom_point(aes(color = Diagnosis.x, shape = trans), alpha = 0.8, size = 2.5, stroke = 1) +

    stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 0.4, color = "firebrick") +

    scale_color_manual(name = "Diagnosis", values = diagnosis_colors) +
    scale_shape_manual(name = "Trans", values = c(16, 17, 15)) +

    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
        breaks = c(0.01, 0.1, 1, 10),
        labels = c("0.01", "0.1", "1", "10")
    ) +
    stat_cor(method = "pearson", label.y.npc = "top", show.legend = FALSE) +
    facet_wrap(~facet_label, scales = "free") +
    xlab("nGD") +
    ylab("Gut abundance (%)") +
    theme_bw() +
    theme(
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 8),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold", size = 8)
    )

print(p_Extended_Data_Fig2)
ggsave("Extended_Data_Fig2.png", p_Extended_Data_Fig2, width = 12, height = 10, dpi = 300)
ggsave("Extended_Data_Fig2.pdf", p_Extended_Data_Fig2, width = 11.5, height = 10)
