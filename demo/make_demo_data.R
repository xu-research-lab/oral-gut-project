# Create a small real-data subset for the demo (18 samples: 6 CRC, 6 GC, 6 EC).
# Run from the demo/ directory: Rscript make_demo_data.R
load("../fig2/shared_abundance_counts.RData")
ff <- read.csv("../data/clinical_features/features_filled_with_metadata.csv",
               check.names = FALSE, stringsAsFactors = FALSE)

sa <- as.data.frame(shared_abundance_counts)[, c("id", "treatment", "Diagnosis", "gut_sum", "shared_counts")]
meta <- ff[, c("id", "Gender", "Age", "BMI_Status",
               "MSI_Status", "Tumor_Staging", "Treatment_class", "smoking_label")]
d <- merge(sa, meta, by = "id")
d <- d[d$treatment == "non-PMA" & d$Diagnosis %in% c("CRC", "GC", "EC"), ]
d <- d[!grepl("_B$|_C$", d$id), ]
d <- d[complete.cases(d[, c("shared_counts", "Age", "Gender", "BMI_Status",
                            "MSI_Status", "Tumor_Staging", "Treatment_class",
                            "smoking_label")]), ]

pick <- function(x, n = 6) {
  x <- x[order(x$Gender), ]
  m <- x[x$Gender == "M", ]
  f <- x[x$Gender == "F", ]
  out <- rbind(head(m, 3), head(f, 3))
  if (nrow(out) < n) {
    rest <- x[!x$id %in% out$id, ]
    out <- rbind(out, head(rest, n - nrow(out)))
  }
  out[seq_len(min(n, nrow(out))), ]
}

set.seed(1)
sel <- do.call(rbind, lapply(c("CRC", "GC", "EC"), function(dx) pick(d[d$Diagnosis == dx, ])))
sel <- sel[, c("id", "Diagnosis", "shared_counts", "gut_sum",
               "Age", "Gender", "BMI_Status", "MSI_Status", "Tumor_Staging",
               "Treatment_class", "smoking_label")]
write.csv(sel, "demo_data.csv", row.names = FALSE)
cat("wrote demo_data.csv:", nrow(sel), "samples\n")
print(table(sel$Diagnosis, sel$Gender))
