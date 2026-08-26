# 01d_idtaxa_v34.R — IdTaxa with V3-V4-region reference (primer-extracted)
library(DECIPHER)
library(Biostrings)
setwd("data/prjna1189573")
seqtab <- readRDS("../../output/m2_microbiome/seqtab_nochim.rds")
ref <- readDNAStringSet("../silva/silva_train_v34.fa.gz")
tax_str <- paste0("Root;", sub("^[^ ]+ ", "", names(ref)))
trfile <- "../../output/m2_microbiome/idtaxa_v34_training.rds"
if (file.exists(trfile)) {
  cat("loading existing v34 training set\n")
  trainingSet <- readRDS(trfile)
} else {
  cat("training v34 classifier (", length(ref), "seqs)...\n")
  trainingSet <- LearnTaxa(ref, tax_str)
  saveRDS(trainingSet, trfile)
}
asv_dna <- DNAStringSet(colnames(seqtab)); names(asv_dna) <- colnames(seqtab)
ids <- IdTaxa(asv_dna, trainingSet, type = "collapsed", strand = "top", threshold = 0)
saveRDS(ids, "../../output/m2_microbiome/idtaxa_v34_ids.rds")
cat("classified:", length(ids), "ASVs\n")
ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxa <- matrix(NA_character_, nrow=length(ids), ncol=7,
               dimnames=list(colnames(seqtab), ranks))
for (i in seq_along(ids)) {
  s <- ids[[i]]
  toks <- strsplit(s, "; ", fixed=TRUE)[[1]]
  nms <- sub(" [%[].*$", "", toks)
  nms <- nms[nms != "Root"]
  n <- min(length(nms), 7)
  taxa[i, seq_len(n)] <- nms[seq_len(n)]
}
taxa[grepl("^unclassified", taxa)] <- NA
write.table(taxa, "output/m2_microbiome/taxonomy_v34.tsv", sep="\t", quote=FALSE)
genus <- taxa[,"Genus"]
print(head(sort(table(genus), decreasing=TRUE), 12))
cat("saved\n")
