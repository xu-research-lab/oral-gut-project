#!/usr/bin/env Rscript
#
# merge_ngd_analysis.R
# Description: Merge and analyze nGD results from StrainPhlAn and Panphlan methods
# Usage: Rscript merge_ngd_analysis.R

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

###############################################################################
# Step 1: Load data
###############################################################################

# Load StrainPhlAn distance data and panphlan data
load("../data/data processing/distance_strainphlan_panphlan.RData")

###############################################################################
# Step 2: Merge results from both methods and select minimal distance as final genetic distance
###############################################################################
# Merge distance from both methods
dis_merged <- full_join(
  dis_panphlan %>% 
    select(ID1, ID2, SGB, centered,bin_cluster,Species) %>%
    rename(centered_panphlan = centered),
  
  dis_strainphlan %>% 
    select(ID1, ID2, SGB, centered,Species) %>%
    rename(centered_strainphlan = centered),
  by = c("ID1", "ID2", "SGB")
)

dis_merged <- dis_merged %>%
  mutate(centered.min = case_when(
    !is.na(centered_strainphlan) & !is.na(centered_panphlan) ~ pmin(centered_strainphlan, centered_panphlan),
    !is.na(centered_strainphlan) ~ centered_strainphlan,
    !is.na(centered_panphlan) ~ centered_panphlan,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(centered.min))

###############################################################################
# Step 3: Filter distance pairs for oral-gut and technical replicate comparisons
###############################################################################

# Filter distance pairs to include:
# 1. Oral-gut pairs (Saliva-Feces comparisons)
# 2. Paired samples with/without PMA treatment or technical replicates
dis_merged<-subset(dis_merged,(grepl("Saliva_",ID1)&grepl("Feces_",ID2))|
                     (grepl("Saliva_",ID2)&grepl("Feces_",ID1))|
                     gsub("PMA_|_A|_B|_C","",ID1) == ID2|
                     gsub("PMA_|_A|_B|_C","",ID2) == ID1)


#Clean up sample names and define groups for merged distance
dis_merged <- dis_merged %>%
  mutate(
    # Clean sample IDs by removing prefixes and suffixes
    clean_id1 = gsub("Feces_|Saliva_", "", ID1),
    clean_id2 = gsub("Feces_|Saliva_", "", ID2),
    
    # Define comparison groups based on sample relationships
    group = ifelse(
      clean_id1 == clean_id2 & 
        ((grepl("^Feces_", ID1) & grepl("^Saliva_", ID2)) | 
           (grepl("^Saliva_", ID1) & grepl("^Feces_", ID2))),
      "within individual oral_gut",
      ifelse(
        gsub("PMA_|_A|_B|_C","",clean_id1) == clean_id2|
          gsub("PMA_|_A|_B|_C","",clean_id2) == clean_id1,
        "within individual",
        "between individuals"
      )
    )
  ) %>%
  select(-clean_id1, -clean_id2)


###############################################################################
# Step 4: Calculate transmission thresholds for each SGB using Youden's index
###############################################################################

cat("Calculating transmission thresholds for each SGB...\n")

calculate_youden_threshold <- function(sgb_data, bin_name) {
  # Extract distances for within-individual and between-individual comparisons
  within_dist <- sgb_data$centered.min[sgb_data$group == "within individual"]
  between_dist <- sgb_data$centered.min[sgb_data$group == "between individuals"]
  
  # Remove NA and infinite values
  within_dist <- within_dist[is.finite(within_dist)]
  between_dist <- between_dist[is.finite(between_dist)]
  
  if (length(within_dist) < 5 || length(between_dist) < 5) {
    warning(paste("Insufficient data for SGB:", bin_name))
    return(NA)
  }
  
  # Define density estimation range
  density_range <- c(mean(within_dist, na.rm = TRUE), 
                     mean(between_dist, na.rm = TRUE))
  
  # Estimate densities
  da <- density(within_dist, from = density_range[1], to = density_range[2], na.rm = TRUE)
  db <- density(between_dist, from = density_range[1], to = density_range[2], na.rm = TRUE)
  
  # Find Youden's index (point where densities are closest)
  youden_index <- da$x[which.min(abs(da$y - db$y))]
  
  return(youden_index)
}

# Calculate thresholds for each SGB
thresholds <- dis_merged %>%
  group_by(SGB) %>%
  group_modify(~ {
    data.frame(threshold = calculate_youden_threshold(.x, .y$SGB))
  }) %>%
  ungroup()

# Merge thresholds back to the main dataframe
dis_merged<- dis_merged %>%
  left_join(thresholds, by = "SGB")

###############################################################################
# Step 5: Determine transmission events
###############################################################################
trans_merged<-subset(dis_merged,group=="within individual oral_gut")

trans_merged <- trans_merged %>%
  mutate(
    trans = case_when(
      centered.min<=threshold ~ "YES",
      centered.min>threshold ~ "NO"
    )
  )
cat("Determining transmission events...\n")

save(dis_merged,trans_merged,file="../results/processed_merge_dis.RData")


