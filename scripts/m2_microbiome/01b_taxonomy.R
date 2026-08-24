# 01b_taxonomy.R — assignTaxonomy on raw seqtab (DNA colnames retained)
library(dada2)
setwd("data/prjna1189573")
seqtab_nochim <- readRDS("../../output/m2_microbiome/seqtab_nochim.rds")
cat("seqtab dims:", dim(seqtab_nochim), "\n")
# keep only valid DNA sequences >= 50 nt
ok <- nchar(colnames(seqtab_nochim)) >= 50 & !grepl("[^ACGT]", colnames(seqtab_nochim))
cat("valid sequences:", sum(ok), "/", length(ok), "\n")
seqtab_filt <- seqtab_nochim[, ok, drop = FALSE]
taxa <- assignTaxonomy(seqtab_filt, "../silva/silva_train_acgt.fa.gz",
                       multithread = 2, tryRC = FALSE, verbose = TRUE)
cat("taxonomy done\n")
rownames(taxa) <- colnames(seqtab_filt)
write.table(t(seqtab_filt), "../../output/m2_microbiome/asv_table.tsv", sep = "\t", quote = FALSE)
write.table(taxa, "../../output/m2_microbiome/taxonomy.tsv", sep = "\t", quote = FALSE)
cat("saved\n")
