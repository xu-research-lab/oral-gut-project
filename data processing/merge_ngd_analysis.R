#!/usr/bin/env Rscript
#
# merge_ngd_analysis.R
# Description: Merge and analyze nGD results from StrainPhlAn and assembly-based methods
# Usage: Rscript merge_ngd_analysis.R

# Load required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# Set seed for reproducibility
set.seed(123)

###############################################################################
# Step 1: Load and preprocess data
###############################################################################

cat("Loading and preprocessing data...\n")

# Load StrainPhlAn distance data
load("distance.RData")
dis_strainphlan <- dis
rm(dis)

# Load assembly-based distance data
load("kept_dis.RData")
dis_assembly <- dis
rm(dis)

# Clean up sample names and define groups for StrainPhlAn data
dis_strainphlan <- dis_strainphlan %>%
  mutate(
    # Clean sample IDs by removing prefixes and suffixes
    clean_id1 = gsub("Feces_|PMA_|Saliva_|Dental|_A|_B|_C", "", ID1),
    clean_id2 = gsub("Feces_|PMA_|Saliva_|_A|_B|_C", "", ID2),
    
    # Define groups based on sample relationships
    group = case_when(
      clean_id1 == clean_id2 ~ "within individual",
      gsub("Feces_", "Saliva_", ID1) == ID2 ~ "oral_gut",
      TRUE ~ "between individuals"
    )
  ) %>%
  select(-clean_id1, -clean_id2)

# Extract oral-gut transmission pairs
trans_strainphlan <- dis_strainphlan %>%
  filter(group == "oral_gut")
  
# Clean up sample names and define groups for StrainPhlAn data
trans_assembly <- trans_assembly %>%
  mutate(
    # Clean sample IDs by removing prefixes and suffixes
    clean_id1 = gsub("Feces_|PMA_|Saliva_|Dental|_A|_B|_C", "", ID1),
    clean_id2 = gsub("Feces_|PMA_|Saliva_|_A|_B|_C", "", ID2),
    
    # Define groups based on sample relationships
    group = case_when(
      clean_id1 == clean_id2 ~ "within individual",
      gsub("Feces_", "Saliva_", ID1) == ID2 ~ "oral_gut",
      TRUE ~ "between individuals"
    )
  ) %>%
  select(-clean_id1, -clean_id2)

# Extract oral-gut transmission pairs
trans_assembly <- dis_strainphlan %>%
  filter(group == "oral_gut")



###############################################################################
# Step 2: Merge results from both methods and select minimal distance as final genetic distance
###############################################################################

cat("Merging results from both methods...\n")

# Merge transmission results
trans_merged <- full_join(
  trans_panphlan %>% 
    select(ID1, ID2, SGB, centered) %>%
    rename(centered_panphlan = centered),
  
  trans_strainphlan %>% 
    select(ID1, ID2, SGB, centered) %>%
    rename(centered_strainphlan = centered),
  by = c("ID1", "ID2", "SGB")
)

trans_merged<-trans_merged%>%mutate(centered.min=min(centered_panphlan,centered_strainphlan))


###############################################################################
# Step 3: Calculate transmission thresholds for each species using Youden's index
###############################################################################

cat("Calculating transmission thresholds for each species...\n")

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
thresholds <- trans_merged %>%
  group_by(SGB) %>%
  group_modify(~ {
    data.frame(threshold = calculate_youden_threshold(.x, .y$SGB))
  }) %>%
  ungroup()

# Merge thresholds back to the main dataframe
trans_merged<- trans_merged %>%
  left_join(thresholds, by = "SGB") %>%
  rename(threshold_strainphlan = threshold)

###############################################################################
# Step 4: Determine transmission events
###############################################################################

cat("Determining transmission events...\n")

trans_merged <- trans_merged %>%
  mutate(
    # Check if distances are below thresholds
    trans_panphlan = ifelse(centered_panphlan <= thresholds, TRUE, FALSE),
    trans_strainphlan = ifelse(centered_strainphlan <= thresholds, TRUE, FALSE),
    
    # Handle NA values
    trans_panphlan = ifelse(is.na(trans_panphlan), FALSE, trans_panphlan),
    trans_strainphlan = ifelse(is.na(trans_strainphlan), FALSE, trans_strainphlan),
    
    # Combine results from both methods
    trans_count = as.numeric(trans_panphlan) + as.numeric(trans_strainphlan),
    
    # Final transmission call
    transmission = case_when(
      trans_count >= 1 ~ "YES",
      trans_count == 0 ~ "NO",
      TRUE ~ "UNCERTAIN"
    )
  )

