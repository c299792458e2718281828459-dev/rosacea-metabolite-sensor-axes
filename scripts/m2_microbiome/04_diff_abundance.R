# 04_diff_abundance.R — genus-level differential abundance (rosacea vs control)
# Requires a sample-group mapping file: data/prjna1189573/sample_groups.tsv
#   columns: sample (library name) <tab> group (Rosacea|Control)
# Falls back to published Joura Table 2 values if mapping absent.
library(dada2)

setwd("data/prjna1189573")

seqtab <- readRDS("../../output/m2_microbiome/seqtab_nochim.rds")
taxa <- read.delim("../../output/m2_microbiome/taxonomy.tsv", row.names = 1)

# remove plant-organelle ASVs (chloroplast/mitochondria)
is_org <- grepl("Chloroplast|Mitochondria", taxa[colnames(seqtab), "Order"]) | 
          grepl("Chloroplast|Mitochondria", taxa[colnames(seqtab), "Family"])
is_org[is.na(is_org)] <- FALSE
seqtab2 <- seqtab[, !is_org, drop = FALSE]
cat("ASVs removed (plant organelles):", sum(is_org), "\n")
genus <- taxa[colnames(seqtab2), "Genus"]
genus[is.na(genus) | genus == "NA"] <- "Unclassified"
gcounts <- t(rowsum(t(seqtab2), group = genus, reorder = TRUE))
grel <- sweep(gcounts, 1, rowSums(gcounts), "/")

# median genus table (all samples)
gm <- sort(apply(grel, 2, median), decreasing = TRUE)
cat("Top 15 genera (median %):\n"); print(round(head(gm, 15), 2))
write.table(grel, "../../output/m2_microbiome/genus_relabund.tsv", sep = "\t", quote = FALSE)

# compare with Joura published values
pub <- c(Cutibacterium = 35.14, Neisseria = 18.65, Staphylococcus = 13.97,
         Corynebacterium = 10.56, Bacteroides = 2.94, Faecalibacterium = 1.28,
         Prevotella = 0.47, Blautia = 0.83, Ruminococcus = 0.69, Subdoligranulum = 0.28)
cat("\nOur median vs Joura published (rosacea):\n")
for (g in names(pub)) {
  ours <- ifelse(g %in% names(gm), round(gm[g], 2), NA)
  cat(sprintf("  %-18s ours=%5.2f  published=%5.2f\n", g, ours, pub[g]))
}

# --- differential abundance if group labels available -------------------------
mapfile <- "sample_groups.tsv"
if (file.exists(mapfile)) {
  grp <- read.delim(mapfile, stringsAsFactors = FALSE)
  rownames(grp) <- grp$sample
  grp <- grp[rownames(grel), ]
  stopifnot(all(grp$group %in% c("Rosacea", "Control")))
  res_list <- list()
  for (g in colnames(grel)) {
    x <- grel[grp$group == "Rosacea", g]
    y <- grel[grp$group == "Control", g]
    if (sum(x > 0) >= 3 && sum(y > 0) >= 3) {
      wt <- wilcox.test(x, y)
      res_list[[length(res_list) + 1]] <- data.frame(
        genus = g, median_Rosacea = median(x) * 100, median_Control = median(y) * 100,
        p_value = wt$p.value)
    }
  }
  res <- do.call(rbind, res_list)
  res$fdr <- p.adjust(res$p_value, method = "BH")
  res <- res[order(res$p_value), ]
  write.table(res, "../../output/m2_microbiome/genus_diff_abundance.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nDifferential genera (Wilcoxon, BH-FDR):\n")
  print(res[res$fdr < 0.05, ], row.names = FALSE)
} else {
  cat("\nNOTE: no sample_groups.tsv — differential testing skipped; use published Table 2.\n")
}
cat("done\n")
