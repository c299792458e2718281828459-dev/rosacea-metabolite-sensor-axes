# 01c_idtaxa.R — taxonomy assignment with DECIPHER IdTaxa (SILVA v138.1 reference)
library(DECIPHER)
library(Biostrings)

setwd("data/prjna1189573")

seqtab <- readRDS("../../output/m2_microbiome/seqtab_nochim.rds")
cat("seqtab dims:", dim(seqtab), "\n")

# --- Load reference (7-rank SILVA, pure ACGT, 400bp) -------------------------
ref <- readDNAStringSet("../silva/silva_train_acgt.fa.gz")
tax_names <- names(ref)
tax_str <- paste0("Root;", sub("^[^ ]+ ", "", tax_names))  # DECIPHER requires Root; prefix
cat("reference sequences:", length(ref), "\n")

# --- Train classifier (skip if already trained) --------------------------------
trfile <- "../../output/m2_microbiome/idtaxa_training.rds"
if (file.exists(trfile)) {
  cat("loading existing training set\n")
  trainingSet <- readRDS(trfile)
} else {
  cat("training classifier...\n")
  trainingSet <- LearnTaxa(ref, tax_str)
  saveRDS(trainingSet, trfile)
}

# --- Classify ASVs --------------------------------------------------------------
asv_dna <- DNAStringSet(colnames(seqtab))
names(asv_dna) <- colnames(seqtab)
ids <- IdTaxa(asv_dna, trainingSet, type = "collapsed", strand = "both", threshold = 0)
saveRDS(ids, "../../output/m2_microbiome/idtaxa_ids.rds")
cat("classified:", length(ids), "ASVs\n")

# --- Format as DADA2-style taxonomy table (positional: Root,K,P,C,O,F,G,S) -----
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa <- matrix(NA_character_, nrow = length(ids), ncol = 7,
               dimnames = list(colnames(seqtab), ranks))
for (i in seq_along(ids)) {
  x <- ids[[i]]
  if (is.null(x)) next
  tv <- x$taxon
  if (is.null(tv) || length(tv) == 0) next
  # drop Root; take up to 7 following levels
  tv <- tv[tv != "Root"]
  n <- min(length(tv), 7)
  taxa[i, seq_len(n)] <- tv[seq_len(n)]
  # "unclassified_*" entries become NA
}
taxa[grepl("^unclassified_", taxa)] <- NA
write.table(taxa, "../../output/m2_microbiome/taxonomy.tsv", sep = "\t", quote = FALSE)
write.table(t(seqtab), "../../output/m2_microbiome/asv_table.tsv", sep = "\t", quote = FALSE)
cat("saved taxonomy.tsv and asv_table.tsv\n")
