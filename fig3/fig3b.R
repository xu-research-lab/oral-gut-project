# Load necessary R packages
# If not already installed, use install.packages("package_name") to install
library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggnewscale)

# Load your data file
# Please ensure the "fig3b.RData" file is in your R working directory
load("fig3b.RData")

# --- 1. Data Preparation and Filtering ---

# Define the diagnosis categories to be filtered from the original data
original_diagnosis_levels <- c("CRC", "GC", "EC", "RA untreated", "HC")

# Define the final label order for the legend
final_diagnosis_levels <- c("CRC", "GC", "EC", "RA", "HC")

# Define colors consistent with the category order
diagnosis_colors <- c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb", "#a6d854")
# Match colors to the new label names to ensure correct color assignment
names(diagnosis_colors) <- final_diagnosis_levels
# Step 1a: Create the base data frame for plotting
plot_data_full <- trans2 %>%
    filter(treatment.y == "non-PMA") %>%
    filter(Diagnosis.x %in% original_diagnosis_levels) %>%
    mutate(Diagnosis.x = recode(Diagnosis.x, "RA untreated" = "RA")) %>%
    mutate(Diagnosis.x = factor(Diagnosis.x, levels = final_diagnosis_levels)) %>%
    filter(!is.na(centered.min) & !is.na(gut) & !is.na(species) & !is.na(trans))


# Next, clean the 'species' column, removing taxonomic prefixes and SGB suffixes
plot_data_full$species <- gsub("s__", "", plot_data_full$species)
plot_data_full$species <- gsub("_SGB.*", "", plot_data_full$species)
plot_data_full$species <- gsub("_unclassified", "", plot_data_full$species)

# Then, create the 'facet_label' for faceting
# Combine the cleaned species name with the extracted SGB information
plot_data_full$facet_label <- paste(plot_data_full$species, plot_data_full$SGB, sep = "\n")
plot_data_full$facet_label

# --- 2. Plot and save the first figure (using full data) ---

p_combined <- ggplot(plot_data_full, aes(x = centered.min, y = gut, color = Diagnosis.x)) +
    geom_point(aes(shape = trans), alpha = 0.95, size = 2.5, stroke = 1) +
    stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, linewidth = 0.6) +
    scale_color_manual(name = "Diagnosis", values = diagnosis_colors) +
    scale_shape_manual(name = "same strain", values = c(16, 17, 15)) +
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
        breaks = c(0.01, 0.1, 1, 10),
        labels = c("0.01", "0.1", "1", "10")
    ) +
    # The stat_cor here inherits aes(color=...), so it calculates by disease group
    stat_cor(method = "pearson", label.y.npc = 0.32, label.x.npc = 0, show.legend = FALSE, size = 3.2) +
    xlab("nGD") +
    ylab("Gut abundance (%)") +
    theme_bw() +
    theme(legend.spacing.y = unit(0.1, "cm"))


# Print the first plot
print(p_combined)
# Save the first plot
ggsave("fig3b.png", p_combined, width = 4.2, height = 2.8)
ggsave("fig3b.pdf", p_combined, width = 4.2, height = 2.7)
