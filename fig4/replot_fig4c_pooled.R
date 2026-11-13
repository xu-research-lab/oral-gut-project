# Re-plot fig4c_pooled.pdf only (uses saved pooled results), marking only
# significant comparisons with "*".
suppressMessages({library(ggplot2)})

outdir <- "fig4_adjusted"
e <- new.env(); load("D:/oral-gut/code_data_all/fig3/Fig3cdef.RData", envir = e)
features <- get("features", e)

d <- features
d$subject_id <- gsub("_[ABC]$", "", gsub("^PMA_", "", gsub("^Feces_", "", d$id)))
d$MMR_MSI <- factor(ifelse(d$MSI_Status == "MSI-high", "dMMR",
                    ifelse(is.na(d$MSI_Status), NA, "pMMR")), levels = c("pMMR", "dMMR"))
d$Stage   <- factor(ifelse(d$Tumor_Staging %in% c("IV", "IVA"), "IV", "I-III"), levels = c("I-III", "IV"))
d$Tx      <- factor(ifelse(grepl("combined|chemo|XELOX", d$Treatment), "chemo", "ICI"), levels = c("ICI", "chemo"))
d$BMI_cat <- factor(d$BMI_Status, levels = c("normal", "lean", "Obesity"))
pl <- read.csv("D:/oral-gut/addtional analysis/results/clinical_features/patient_level_metadata.csv",
               check.names = FALSE, stringsAsFactors = FALSE)
d$smoking <- factor(pl$smoking_label[match(d$subject_id, as.character(pl$subject_id))], levels = c("No", "Yes"))

d <- d[d$Diagnosis %in% c("CRC", "GC", "EC") & !is.na(d$Response_6) & d$Response_6 %in% c("R", "NR") &
       !grepl("_B$|_C$", d$id) & d$treatment %in% c("non-PMA", "PMA"), ]
d$Response_6 <- factor(d$Response_6, levels = c("NR", "R"))
d$Diagnosis  <- factor(d$Diagnosis, levels = c("CRC", "GC", "EC"))
d$treatment  <- factor(d$treatment, levels = c("non-PMA", "PMA"))

res <- read.csv(file.path(outdir, "fig4c_pooled_results.csv"), check.names = FALSE, stringsAsFactors = FALSE)
order_p <- c("CRC_nonPMA","CRC_PMA","GC_nonPMA","GC_PMA","EC_nonPMA","EC_PMA")

ann <- data.frame(
  Diagnosis = rep(c("CRC","GC","EC"), each = 2),
  treatment = rep(c("non-PMA","PMA"), 3), stringsAsFactors = FALSE)
ann$panel <- paste(ann$Diagnosis, ifelse(ann$treatment == "non-PMA", "nonPMA", "PMA"), sep = "_")
ann$Diagnosis <- factor(ann$Diagnosis, levels = c("CRC","GC","EC"))
ann$treatment <- factor(ann$treatment, levels = c("non-PMA","PMA"))
ann <- merge(ann, res[, c("panel","p")], by = "panel", all.x = TRUE)
ann$lab <- ifelse(!is.na(ann$p) & ann$p < 0.05, "*", "")

p <- ggplot(d, aes(x = Diagnosis, y = transmitted.abundance.in.Gut, color = Response_6)) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.6, position = position_dodge(width = 0.6)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6), size = 1.6) +
  facet_wrap(~treatment) +
  scale_color_manual(name = "Response", values = c("NR" = "#E41A1C", "R" = "#377EB8")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  geom_text(data = ann, aes(x = Diagnosis, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.4, size = 11.0) +
  labs(x = NULL, y = "Transmitted abundance in gut") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 14),
        axis.text = element_text(size = 13), axis.title = element_text(size = 14),
        legend.title = element_text(size = 13), legend.text = element_text(size = 13))
ggsave(file.path(outdir, "fig4c_pooled.pdf"), p, width = 5.5, height = 3.4)
ggsave(file.path(outdir, "fig4c_pooled.png"), p, width = 5.5, height = 3.4, dpi = 300)
cat("done fig4c replot\n")
