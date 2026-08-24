# 03_phylogenetic_summary.R — genus-level relative abundance from DADA2 ASV table
# Compares our re-analysis with Joura 2024 published Table 2 (skin).
library(dada2)

setwd("data/prjna1189573")
seqtab <- readRDS("../../output/m2_microbiome/seqtab_nochim.rds")
taxa <- read.delim("../../output/m2_microbiome/taxonomy.tsv", row.names = 1)

# aggregate to genus
genus <- taxa[colnames(seqtab), "Genus"]
genus[is.na(genus) | genus == "NA"] <- "Unclassified"
gtable <- t(rowsum(t(seqtab), group = genus, reorder = TRUE))
grel <- sweep(gtable, 1, rowSums(gtable), "/") * 100

# sample library names (from rownames)
write.table(grel, "../../output/m2_microbiome/genus_relabund.tsv", sep = "\t", quote = FALSE)
write.table(gtable, "../../output/m2_microbiome/genus_counts.tsv", sep = "\t", quote = FALSE)

# dominant genera (all 31 skin samples pooled median)
gm <- apply(grel, 2, median)
gm <- sort(gm, decreasing = TRUE)
cat("Top 12 genera (median % across 31 skin samples):\n")
print(round(head(gm, 12), 2))

# compare with Joura Table 2 published medians
pub <- c(Cutibacterium = 35.14, Neisseria = 18.65, Staphylococcus = 13.97,
         Corynebacterium = 10.56, Bacteroides = 2.94, Faecalibacterium = 1.28,
         Prevotella = 0.47, Blautia = 0.83, Ruminococcus = 0.69, Subdoligranulum = 0.28)
cat("\nOur median vs Joura published (rosacea patients):\n")
for (g in names(pub)) {
  ours <- ifelse(g %in% names(gm), round(gm[g], 2), NA)
  cat(sprintf("  %-18s ours=%5.2f  published=%5.2f\n", g, ours, pub[g]))
}
cat("done\n")
