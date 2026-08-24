# 01_dada2_pipeline.R — PRJNA1189573 skin samples (Joura et al. 2024)
# DADA2: filter -> error models -> ASV inference -> merge -> chimera removal -> taxonomy
# Case/control labels pending (author request); this script builds the ASV table only.
# Output: output/m2_microbiome/asv_table.tsv, taxonomy.tsv, track.tsv

library(dada2)

setwd("data/prjna1189573")
dir.create("../../output/m2_microbiome", showWarnings = FALSE, recursive = TRUE)

# --- File lists --------------------------------------------------------------
fwd <- sort(list.files("fastq", pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
rev <- sort(list.files("fastq", pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(fwd) == length(rev), length(fwd) == 31)
sample_names <- sub("_1\\.fastq\\.gz$", "", basename(fwd))

filt_fwd <- file.path("fastq/filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filt_rev <- file.path("fastq/filtered", paste0(sample_names, "_R_filt.fastq.gz"))

# --- Quality-filter and trim (V3-V4, MiSeq 2x300: trunc based on quality) ---
out <- filterAndTrim(fwd, filt_fwd, rev, filt_rev,
                     truncLen = c(280, 260), trimLeft = c(0, 0),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread = 4)
head(out)

# --- Error models -------------------------------------------------------------
errF <- learnErrors(filt_fwd, multithread = 4, randomize = TRUE)
errR <- learnErrors(filt_rev, multithread = 4, randomize = TRUE)

# --- Sample inference ---------------------------------------------------------
dadaF <- dada(filt_fwd, err = errF, multithread = 4)
dadaR <- dada(filt_rev, err = errR, multithread = 4)

# --- Merge pairs --------------------------------------------------------------
merged <- mergePairs(dadaF, filt_fwd, dadaR, filt_rev, verbose = TRUE)

# --- Sequence table & chimera removal ----------------------------------------
seqtab <- makeSequenceTable(merged)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = 4, verbose = TRUE)
cat("seqs:", sum(seqtab), "-> non-chimeric:", sum(seqtab_nochim), "\n")

# --- Track reads --------------------------------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN),
               sapply(merged, getN), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- sample_names
write.table(track, "../../output/m2_microbiome/track.tsv", sep = "\t", quote = FALSE)

# --- Save raw seqtab (DNA sequences as column names) -------------------------
rownames(seqtab_nochim) <- sample_names
saveRDS(seqtab_nochim, "../../output/m2_microbiome/seqtab_nochim.rds")
cat("seqtab saved\n")
