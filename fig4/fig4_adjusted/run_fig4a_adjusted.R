# ============================================================
# fig4a revised: covariate-adjusted PERMANOVA (adonis2), cancer-only
#   unified with fig4e:
#     - Diagnosis analysis (all cancer types COMBINED):
#         adjust for ALL available factors
#         Age + Gender + BMI_Status + smoking_label + MSI + Stage + Treatment
#     - Response_6 analysis (WITHIN each cancer, separated):
#         adjust for the 3 factors  MSI_Status + Tumor_Staging + Treatment_class
#   coding collapsed as in fig4e: MMR d/p, Stage IV vs I-III, Tx chemo/ICI
# Run from: D:/oral-gut/code_data_all/fig4
# ============================================================
suppressMessages({library(vegan); library(dplyr); library(ggplot2)})
set.seed(123)
outdir <- "fig4_adjusted"; dir.create(outdir, showWarnings = FALSE)

data_types <- c("non_PMA", "PMA")
sample_groups <- c("oral", "gut")
dist_metrics_names <- c("Bray-Curtis"="bray-curtis", "Jaccard"="jaccard",
                        "Unweighted UniFrac"="unweighted-unifrac", "Weighted UniFrac"="weighted-unifrac")
covs_all   <- c("Age","Gender","BMI_Status","smoking_label","MMR_MSI","Stage","Tx")
covs_three <- c("MMR_MSI","Stage","Tx")

# ---- metadata + covariates ----
meta_data_full <- read.csv("fig4a.csv", sep=",", row.names=1, header=TRUE, check.names=FALSE)
meta_data_full <- meta_data_full[meta_data_full$Diagnosis %in% c("CRC","GC","EC"), ]
rn <- rownames(meta_data_full)
meta_data_full$treatment <- ifelse(grepl("PMA", rn), "PMA", "non-PMA")
sid <- sub("^Saliva_","",rn); sid <- sub("^PMA_","",sid); sid <- sub("_.*$","",sid)
meta_data_full$subject_id <- sid

ff <- read.csv("D:/oral-gut/addtional analysis/results/clinical_features/features_filled_with_metadata.csv",
               check.names=FALSE, stringsAsFactors=FALSE)
ff$key <- paste(ff$subject_id, ff$treatment)
ff <- ff[!duplicated(ff$key), ]
ff$MMR_MSI <- ifelse(ff$MSI_Status=="MSI-high","dMMR","pMMR")
ff$Stage   <- ifelse(ff$Tumor_Staging %in% c("IV","IVA"),"IV","I-III")
ff$Tx      <- ifelse(grepl("chemo",ff$Treatment_class),"chemo","ICI")
meta_data_full$key <- paste(meta_data_full$subject_id, meta_data_full$treatment)
meta_data_full$row_id <- rn
meta <- merge(meta_data_full, ff[, c("key","Age","BMI_Status","MMR_MSI","Stage","Tx","smoking_label")],
              by="key", all.x=TRUE, sort=FALSE)
meta <- meta[match(rn, meta$row_id), ]; rownames(meta) <- meta$row_id
meta$Age <- as.numeric(meta$Age)
for (v in c("Gender","BMI_Status","smoking_label","MMR_MSI","Stage","Tx")) meta[[v]] <- factor(meta[[v]])
cat("meta rows:", nrow(meta), " complete(all):", sum(complete.cases(meta[,covs_all])), "\n")

meta_oral <- meta
meta_gut  <- meta; rownames(meta_gut) <- gsub("Saliva","Feces", rownames(meta_gut))

adonis_term <- function(D, dat, term, covs) {
  keep <- complete.cases(dat[, c(term, covs), drop=FALSE])
  dat2 <- dat[keep, , drop=FALSE]
  if (nrow(dat2) < 5 || length(unique(dat2[[term]])) < 2) return(NULL)
  smp <- intersect(rownames(dat2), rownames(D)); dat2 <- dat2[smp, , drop=FALSE]
  if (nrow(dat2) < 5 || length(unique(dat2[[term]])) < 2) return(NULL)
  D2 <- as.dist(as.matrix(D)[smp, smp])
  rhs <- if (length(covs) > 0) paste(c(covs, term), collapse=" + ") else term
  ad <- tryCatch(adonis2(as.formula(paste("D2 ~", rhs)), data=dat2, permutations=999), error=function(e) NULL)
  if (is.null(ad) || !(term %in% rownames(ad))) return(NULL)
  c(n=nrow(dat2), R2=ad[term,"R2"], P=ad[term,"Pr(>F)"])
}

all_res <- data.frame()
for (dtype in data_types) {
  dm <- list()
  for (i in seq_along(dist_metrics_names)) {
    fp <- sprintf("total_abundance_profile_%s_%s.tsv", dtype, dist_metrics_names[i])
    if (file.exists(fp)) dm[[names(dist_metrics_names)[i]]] <- read.csv(fp, sep="\t", row.names=1, check.names=FALSE)
  }
  for (sgroup in sample_groups) {
    md <- if (sgroup=="oral") meta_oral else meta_gut
    for (dist_name in names(dm)) {
      D <- dm[[dist_name]]
      u <- adonis_term(D, md, "Diagnosis", character(0))
      a <- adonis_term(D, md, "Diagnosis", covs_all)          # combined -> ALL factors
      all_res <- rbind(all_res, data.frame(Data_Type=dtype, Sample_Group=sgroup, Analysis="Diagnosis",
        Cancer_Type="All", Distance_Metric=dist_name, adj_set="all",
        n=ifelse(is.null(u),NA,u["n"]), R2_unadj=ifelse(is.null(u),NA,u["R2"]), P_unadj=ifelse(is.null(u),NA,u["P"]),
        R2_adj=ifelse(is.null(a),NA,a["R2"]), P_adj=ifelse(is.null(a),NA,a["P"])))
      for (ct in c("CRC","GC","EC")) {
        mdc <- md[md$Diagnosis==ct, , drop=FALSE]
        u <- adonis_term(D, mdc, "Response_6", character(0))
        a <- adonis_term(D, mdc, "Response_6", covs_three)    # separated -> 3 factors
        all_res <- rbind(all_res, data.frame(Data_Type=dtype, Sample_Group=sgroup, Analysis="Response_6",
          Cancer_Type=ct, Distance_Metric=dist_name, adj_set="three",
          n=ifelse(is.null(u),NA,u["n"]), R2_unadj=ifelse(is.null(u),NA,u["R2"]), P_unadj=ifelse(is.null(u),NA,u["P"]),
          R2_adj=ifelse(is.null(a),NA,a["R2"]), P_adj=ifelse(is.null(a),NA,a["P"])))
      }
    }
  }
}
write.csv(all_res, file.path(outdir,"fig4a_adjusted_adonis2_results.csv"), row.names=FALSE)
cat("\n==== Bray-Curtis key results ====\n")
print(all_res[all_res$Distance_Metric=="Bray-Curtis",
      c("Data_Type","Sample_Group","Analysis","Cancer_Type","n","R2_unadj","P_unadj","R2_adj","P_adj")],
      row.names=FALSE, digits=3)

# ---- PCoA (non-PMA, Bray) with adjusted annotations ----
pcoa_data <- list(); pcoa_ann <- list()
for (sgroup in sample_groups) {
  md <- if (sgroup=="oral") meta_oral else meta_gut
  D <- read.csv("total_abundance_profile_non_PMA_bray-curtis.tsv", sep="\t", row.names=1, check.names=FALSE)
  smp <- intersect(rownames(md), rownames(D)); Dm <- as.matrix(D)[smp, smp]
  pc <- cmdscale(as.dist(Dm), k=2, eig=TRUE)
  co <- as.data.frame(pc$points); colnames(co) <- c("Pco1","Pco2")
  eig <- round(pc$eig/sum(pc$eig)*100,1)
  df <- cbind(co, md[rownames(co), ]); dfsub <- subset(df, Response_6 %in% c("R","NR"))
  agg <- aggregate(dfsub[,c("Pco1","Pco2")], by=list(Diagnosis=dfsub$Diagnosis, Response_6=dfsub$Response_6), FUN=mean)
  sda <- aggregate(dfsub[,c("Pco1","Pco2")], by=list(Diagnosis=dfsub$Diagnosis, Response_6=dfsub$Response_6), FUN=sd)
  colnames(sda)[3:4] <- c("Pco1_sd","Pco2_sd")
  pl <- merge(agg, sda, by=c("Diagnosis","Response_6")); pl$Sample_Group <- sgroup
  pcoa_data[[sgroup]] <- pl
  ga <- all_res[all_res$Data_Type=="non_PMA" & all_res$Sample_Group==sgroup & all_res$Analysis=="Diagnosis" & all_res$Distance_Metric=="Bray-Curtis",]
  ra <- all_res[all_res$Data_Type=="non_PMA" & all_res$Sample_Group==sgroup & all_res$Analysis=="Response_6" & all_res$Distance_Metric=="Bray-Curtis",]
  ra$label <- sprintf("R2=%.1f%%, P=%.3f", ra$R2_adj*100, ra$P_adj); ra$y <- c(CRC=-0.24, GC=-0.28, EC=-0.32)[ra$Cancer_Type]
  pcoa_ann[[sgroup]] <- list(resp=ra, diag=ga, axes=c(paste0("PCo1 [",eig[1],"%]"), paste0("PCo2 [",eig[2],"%]")))
}
combined <- do.call(rbind, pcoa_data)
resp_ann <- do.call(rbind, lapply(pcoa_ann, `[[`, "resp"))
diag_ann <- do.call(rbind, lapply(names(pcoa_ann), function(g){ x <- pcoa_ann[[g]]
  data.frame(Sample_Group=g, label=sprintf("Diagnosis (adj)\nR2=%.1f%%, P=%.3f", x$diag$R2_adj*100, x$diag$P_adj), x=-0.4, y=0.22) }))
pcp <- ggplot(combined, aes(Pco1, Pco2, color=Diagnosis)) +
  facet_grid(~Sample_Group) +
  geom_errorbarh(aes(xmin=Pco1-Pco1_sd, xmax=Pco1+Pco1_sd), height=0.01, linewidth=0.3) +
  geom_errorbar(aes(ymin=Pco2-Pco2_sd, ymax=Pco2+Pco2_sd), width=0.01, linewidth=0.3) +
  geom_point(size=8) + geom_text(aes(label=Response_6), color="black", size=3.5) +
  geom_text(data=resp_ann, aes(x=-0.4, y=y, label=label, color=Cancer_Type), size=3.5, show.legend=FALSE, hjust=0) +
  geom_text(data=diag_ann, aes(x=x, y=y, label=label), inherit.aes=FALSE, size=3.5, hjust=0) +
  theme_bw() + scale_color_manual(name="Diagnosis", breaks=c("CRC","GC","EC"),
    values=c("CRC"="#66c2a5","GC"="#fc8d62","EC"="#e78ac3")) +
  labs(x="PCo1", y="PCo2") + theme(strip.background=element_blank())
ggsave(file.path(outdir,"fig4a_revised_pcoa.pdf"), pcp, width=6, height=2.5)
ggsave(file.path(outdir,"fig4a_revised_pcoa.png"), pcp, width=6, height=2.5, dpi=300)
cat("\n== fig4a done: combined->all factors, per-cancer->MSI/Stage/Tx ==\n")
