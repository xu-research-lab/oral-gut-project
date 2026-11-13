# One-vs-rest plot sized to the Fig.1 right-column aspect (w/h = 230.77/450),
# so it fills the column and its x-axis aligns with panel h.
suppressMessages({library(ggplot2); library(dplyr)})
outdir <- "maaslin2_onevsrest"
best <- read.csv(file.path(outdir, "maaslin2_onevsrest_best.csv"), check.names = FALSE)
tax  <- read.csv("tax_id.csv", row.names = 1, check.names = FALSE)
tax$id <- rownames(tax); tax$Species <- sub("^s__", "", tax$Species)
tx <- tax[, c("id", "Species", "group")]; names(tx)[3] <- "oral_group"

pd <- merge(best, tx, by = "id", all.x = TRUE)
pd <- pd %>%
  filter(coef >= 1, !grepl("GGB", Species)) %>%
  mutate(label_color = ifelse(!is.na(oral_group) & oral_group == "oral", "red", "black"),
         group = factor(group, levels = c("CRC", "GC", "EC"))) %>%
  arrange(group, coef) %>%
  mutate(Taxa = paste(Species, id), Taxa = factor(Taxa, levels = unique(Taxa)))

cat("plotted taxa:", nrow(pd), "\n"); print(table(pd$group, useNA = "ifany"))

p <- ggplot(pd, aes(x = coef, y = Taxa, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(name = "Diagnosis",
                    values = c(CRC = "#66c2a5", GC = "#fc8d63", EC = "#e78ac3"),
                    limits = c("CRC", "GC", "EC")) +
  labs(x = "MaAsLin2 coefficient (X vs rest)", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, colour = pd$label_color),
        axis.title.x = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())

target_asp <- 230.77 / 450
h <- 5 / target_asp
cat("height (in):", round(h, 3), " aspect:", round(5/h, 4), "\n")
ggsave(file.path(outdir, "Fig1e_maaslin2_onevsrest_coef1_bl.pdf"), p, height = h, width = 5)
ggsave(file.path(outdir, "Fig1e_maaslin2_onevsrest_coef1_bl.png"), p, height = h, width = 5, dpi = 150)
cat("DONE\n")
