library(reshape2)
library(ggplot2)
library(ggpubr)

jaccard<-read.table("../data/Figure2/all_profiled_metagenome_jaccard.tsv",head=T,row.names=1)
meta<-read.table("../data/Figure2/metadata.tsv",header = T,sep="\t")
# Prepare metadata by selecting relevant columns and filtering samples
meta <- meta[!grepl("PMA_|Treated",meta$ID), ]
meta <- subset(meta, !grepl("_B|_C", meta$ID))

# meta <- subset(meta, treatment != "PMA")
# meta$Row.names <- gsub("Feces_|_A|_B|_C", "", meta$Row.names)

# Standardize diagnosis labels
meta$Diagnosis <- gsub("Untreated", "RA", meta$Diagnosis)

# meta <- subset(meta, Diagnosis %in% c("EC", "GC", "CRC", "RA", "HC"))

plot_weighted_unifrac <- function(weighted_matrix, meta_data, 
                                  oral_prefix = "Saliva_", gut_prefix = "Feces_", 
                                  remove_suffixes = c("_A", "_B", "_C", "d_metagenome"),
                                  output_file = "weighted_unifrac.pdf", 
                                  plot_width = 3.5, plot_height = 2.5,
                                  exclude_groups = "Treated",
                                  diagnosis_levels = c("CRC", "GC", "EC", "RA", "HC")) {
  
  # Parameter validation
  if (!is.matrix(weighted_matrix) && !is.data.frame(weighted_matrix)) {
    stop("weighted_matrix must be a matrix or data frame")
  }
  if (!all(c("ID", "Diagnosis") %in% colnames(meta_data))) {
    stop("meta_data must contain 'ID' and 'Diagnosis' columns")
  }
  
  # Convert matrix to long format
  weighted_melt <- melt(as.matrix(weighted_matrix))
  colnames(weighted_melt) <- c("oral", "gut", "dis")
  
  # Clean sample names by removing prefixes and suffixes
  clean_names <- function(x, prefix) {
    x <- gsub(prefix, "", x)
    for (suffix in remove_suffixes) {
      x <- gsub(suffix, "", x, fixed = TRUE)
    }
    return(x)
  }
  
  weighted_melt$oral <- clean_names(weighted_melt$oral, oral_prefix)
  weighted_melt$gut <- clean_names(weighted_melt$gut, gut_prefix)
  
  # Add group classification (intra vs inter individual)
  weighted_melt$group <- ifelse(weighted_melt$oral == weighted_melt$gut, 
                                "intra-individual", 
                                "inter-individuals")
  
  # Prepare metadata with cleaned names
  meta <- meta_data[, c("ID", "Diagnosis")]
  meta$ID <- clean_names(meta$ID, gut_prefix)
  
  # Merge distance data with metadata
  weighted_meta <- merge(weighted_melt, meta, by.x = "oral", by.y = "ID", all.x = TRUE)
  weighted_meta <- merge(weighted_meta, meta, by.x = "gut", by.y = "ID", all.x = TRUE)
  
  # Filter to keep only pairs with matching diagnoses
  weighted_meta <- subset(weighted_meta, Diagnosis.x == Diagnosis.y)
  
  # Standardize diagnosis labels and set factor levels
  weighted_meta$Diagnosis.x <- gsub("Untreated", "RA", weighted_meta$Diagnosis.x)
  weighted_meta$Diagnosis.x <- factor(weighted_meta$Diagnosis.x, levels = diagnosis_levels)
  
  # Filter out excluded groups if specified
  plot_data <- weighted_meta
  if (!is.null(exclude_groups)) {
    plot_data <- subset(weighted_meta, !Diagnosis.x %in% exclude_groups)
  }
  
  # Create the boxplot
  p <- ggplot(plot_data, aes(Diagnosis.x, dis, color = group)) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(label = "p.signif", show.legend = FALSE, vjust = 1) +
    theme_bw() +
    ylim(c(0.85, 1.07)) +
    theme(legend.position = "top") +
    labs(x = NULL, y = "Jaccard (oral-gut)")
  
  # Save the plot
  ggsave(output_file, plot = p, height = plot_height, width = plot_width)
  
  # Return both the processed data and the plot object
  return(list(data = weighted_meta, plot = p))
}

plot_weighted_unifrac(jaccard,meta,output_file = "../results/oral_gut_jaccard.pdf")