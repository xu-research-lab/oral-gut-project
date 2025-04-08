library(dplyr)
library(purrr)
library(ggplot2)
library(aplot)
library(ggpubr)
library(reshape2)
library(patchwork)
library(scatterpie)
library(tidyr)
library(survival)
library(survminer)


load("Figure4e.RData")


shared_abd <- subset(paired_sum, Group.2 == "oral-gut shared")
features <- merge(features, shared_abd[, c("Group.1", "x")], by.x = "sampleid", by.y = "Group.1", all.x = TRUE)


abundance_threshold <- 2.09818
features <- mutate(features, shared_level = case_when(
  is.na(transmitted.abundance.in.gut) ~ NA_character_, # 明确处理NA值
  transmitted.abundance.in.gut <= abundance_threshold ~ "low",
  transmitted.abundance.in.gut > abundance_threshold ~ "high"
))


shared_group <- features[, c("sampleid", "treatment", "shared_level")]

shared_group <- subset(shared_group, !grepl("_B|_C", sampleid))
shared_group$sampleid <- gsub("Feces_", "", shared_group$sampleid)
shared_group$sampleid <- gsub("PMA_", "", shared_group$sampleid)
shared_group$sampleid <- gsub("_A|_B|_C", "", shared_group$sampleid)
shared_group <- unique(shared_group)


shared_group <- shared_group %>%
  mutate(treatment = gsub("-", "_", treatment))


wide_data <- shared_group %>%
  filter(!is.na(treatment) & !is.na(shared_level)) %>%
  pivot_wider(
    names_from = treatment,
    names_glue = "shared_level_{treatment}",
    values_from = shared_level
  )



wide_data_to_merge <- wide_data %>%
  select(sampleid, starts_with("shared_level_"))


cols_to_remove <- names(wide_data_to_merge)[!names(wide_data_to_merge) %in% "sampleid"]
if (any(cols_to_remove %in% names(melt_L))) {
  melt_L <- melt_L[, !names(melt_L) %in% cols_to_remove]
}

melt_L <- merge(melt_L, wide_data_to_merge, by.x = "SubjectID", by.y = "sampleid", all.x = TRUE)

melt_L1 <- unique(melt_L)
melt_L1 <- subset(melt_L1, grepl("PD", Treatment))

colnames(melt_L1)[colnames(melt_L1) == "Diagnosis.x"] <- "Diagnosis"

fit <- survfit(Surv(Time, Status) ~ Diagnosis + shared_level_non_PMA, data = melt_L1, id = SubjectID)

fit1 <- survfit(Surv(Time, Status) ~ Diagnosis + shared_level_PMA, data = melt_L1, id = SubjectID)

survival_plot <- ggsurvplot(
  fit = fit, data = melt_L1, pval = TRUE, facet.by = "Diagnosis",
  legend.title = "Transmission level", legend.position = "right"
) +
  ylab(label = "Progression-free survival") +
  xlab(label = "Days") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"))

ggsave("/home/hby/huiyi/new/Figure4e.pdf", plot = survival_plot, height = 2.5, width = 6)
