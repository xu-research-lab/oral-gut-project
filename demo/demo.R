# Minimal demo of the analysis pipeline on a small real-data subset.
# Run from the demo/ directory: Rscript demo.R
suppressMessages(library(ggplot2))

d <- read.csv("demo_data.csv", check.names = FALSE, stringsAsFactors = FALSE)
d$Diagnosis <- factor(d$Diagnosis, levels = c("CRC", "GC", "EC"))
d$Gender <- factor(d$Gender)
d$log_shared <- log10(d$shared_counts + 1e-6)

# Representative analysis: shared SGB counts across cancer types,
# adjusted for age and sex (simplified version of the full host-adjusted model).
fit <- lm(log_shared ~ Diagnosis + Age + Gender, data = d)
print(summary(fit))

a <- anova(fit)
print(a)
cat("\nDiagnosis ANOVA p =", signif(a["Diagnosis", "Pr(>F)"], 3), "\n")
cf <- summary(fit)$coefficients
cat("CRC vs GC adjusted p =", signif(cf["DiagnosisGC", "Pr(>|t|)"], 3), "\n")

p <- ggplot(d, aes(Diagnosis, shared_counts, color = Diagnosis)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.8) +
  theme_bw() +
  labs(x = NULL, y = "Shared SGB counts")

ggsave("demo_output.pdf", p, width = 4, height = 3)
ggsave("demo_output.png", p, width = 4, height = 3, dpi = 300)
cat("\nWrote demo_output.pdf and demo_output.png\n")
