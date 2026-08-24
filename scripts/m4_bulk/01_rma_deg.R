# 01_rma_deg.R — GSE65914: RMA normalization + limma differential expression
# Rosacea (all subtypes) vs healthy controls; plus subtype-stratified
# Output: output/m4_bulk/DEG_all.tsv, DEG_ETR.tsv, DEG_PPR.tsv, DEG_PhR.tsv, expr_matrix.tsv

library(affy)
library(limma)

setwd("data/gse65914")
dir.create("../../output/m4_bulk", showWarnings = FALSE, recursive = TRUE)

# --- Read phenotypes ---------------------------------------------------------
pheno <- read.delim("phenotype.tsv", stringsAsFactors = FALSE)
cel_files <- list.files("CEL", pattern = "\\.CEL\\.gz$", full.names = TRUE)
stopifnot(length(cel_files) == 58)

# Map GSM -> file
gsm <- sub("\\.CEL\\.gz$", "", basename(cel_files))
gsm <- sub("_.*", "", gsm)  # keep GSM id part
# Actually CEL filenames are like GSM1611035_GRDS0027_HV1_1.CEL.gz -> GSM id is prefix
gsm_id <- sub("_.*$", "", basename(cel_files))
pheno <- pheno[match(gsm_id, sub('"', '', pheno$gsm)), ]

# --- RMA ----------------------------------------------------------------------
eset <- justRMA(filenames = cel_files, verbose = TRUE)
exprs <- exprs(eset)
colnames(exprs) <- gsm_id

# --- Annotate probes to symbols ----------------------------------------------
library(hgu133plus2.db)
probe_info <- select(hgu133plus2.db, keys = rownames(exprs),
                     columns = c("SYMBOL"), keytype = "PROBEID")
# Collapse: keep probe with highest mean expression per symbol
probe_info <- probe_info[!is.na(probe_info$SYMBOL), ]
keep <- !duplicated(probe_info$SYMBOL)
# (collapsing by max mean expression is more standard)
probe_mean <- rowMeans(exprs)
probe_info$mean <- probe_mean[probe_info$PROBEID]
probe_info <- probe_info[order(probe_info$mean, decreasing = TRUE), ]
probe_info <- probe_info[!duplicated(probe_info$SYMBOL), ]
exprs_sym <- exprs[probe_info$PROBEID, ]
rownames(exprs_sym) <- probe_info$SYMBOL

write.table(exprs_sym, "../../output/m4_bulk/expr_matrix.tsv", sep = "\t", quote = FALSE)

# --- limma DEG: all rosacea vs HC -------------------------------------------
group_all <- ifelse(pheno$subtype == "HC", "HC", "Rosacea")
design <- model.matrix(~ 0 + factor(group_all))
colnames(design) <- c("HC", "Rosacea")
fit <- lmFit(exprs_sym, design)
contr <- makeContrasts(Rosacea - HC, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
deg_all <- topTable(fit2, number = Inf, sort.by = "none")
deg_all$gene <- rownames(deg_all)
write.table(deg_all, "../../output/m4_bulk/DEG_all.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("DEG all: n with |logFC|>1 & adj.P<0.05 =",
    sum(abs(deg_all$logFC) > 1 & deg_all$adj.P.Val < 0.05), "\n")

# --- Subtype-stratified ------------------------------------------------------
subtype <- pheno$subtype
for (s in c("ETR", "PPR", "PhR")) {
  idx <- subtype %in% c("HC", s)
  sub_expr <- exprs_sym[, idx]
  sub_group <- ifelse(subtype[idx] == "HC", "HC", s)
  d <- model.matrix(~ 0 + factor(sub_group))
  colnames(d) <- c("HC", s)
  f <- lmFit(sub_expr, d)
  cc <- makeContrasts(contrasts = paste0(s, " - HC"), levels = d)
  f2 <- contrasts.fit(f, cc); f2 <- eBayes(f2)
  tt <- topTable(f2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  write.table(tt, paste0("../../output/m4_bulk/DEG_", s, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  cat("DEG", s, ": n |logFC|>1 & adj.P<0.05 =",
      sum(abs(tt$logFC) > 1 & tt$adj.P.Val < 0.05), "\n")
}
cat("done\n")
