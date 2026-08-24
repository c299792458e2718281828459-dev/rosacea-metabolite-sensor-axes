# 03_cibersort_gsea.R — CIBERSORT-style nu-SVR deconvolution (LM22) + GSEA (Hallmark)
# CIBERSORT implementation follows Newman et al. Nat Methods 2015:
#   nu-SVR (linear kernel, nu=0.25 for microarray), non-negative coefficients normalized to 1;
#   permutation P-value from row-shuffled signature matrix (100 permutations).
# LM22 signature: Supplementary Table 1 of Newman 2015 (MOESM207.xls, 547 genes x 22 subsets).
# Output: output/m4_bulk/cibersort_fractions.tsv, cibersort_pvalues.tsv, gsea_results.tsv

library(e1071)

# --- Load data ----------------------------------------------------------------
exprs_sym <- as.matrix(read.delim("output/m4_bulk/expr_matrix.tsv", row.names = 1, check.names = FALSE))
pheno <- read.delim("data/gse65914/phenotype.tsv", stringsAsFactors = FALSE)
pheno$gsm <- sub('"', '', pheno$gsm)
pheno <- pheno[match(colnames(exprs_sym), pheno$gsm), ]
grp <- ifelse(pheno$subtype == "HC", "HC", "Rosacea")

LM22 <- read.delim("data/cibersort/LM22.txt", row.names = 1, check.names = FALSE)
genes_common <- intersect(rownames(exprs_sym), rownames(LM22))
cat("LM22 genes in expression matrix:", length(genes_common), "/ 547\n")
B <- as.matrix(LM22[genes_common, ])
Y <- exprs_sym[genes_common, , drop = FALSE]

# --- nu-SVR deconvolution ------------------------------------------------------
deconvolve <- function(y, B, nu = 0.25) {
  m <- svm(B, y, type = "nu-regression", kernel = "linear", nu = nu, scale = TRUE)
  w <- t(m$coefs) %*% m$SV
  w <- as.numeric(w)
  w[w < 0] <- 0
  s <- sum(w)
  if (s > 0) w <- w / s
  names(w) <- colnames(B)
  w
}

n_perm <- 100
ncores <- min(8, parallel::detectCores())

fit_one_sample <- function(i, Y, B, n_perm) {
  y <- Y[, i]
  w <- deconvolve(y, B)
  rms_obs <- sqrt(mean((B %*% w - y)^2))
  rms_perm <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    Bp <- B[sample(nrow(B)), , drop = FALSE]
    wp <- deconvolve(y, Bp)
    rms_perm[p] <- sqrt(mean((Bp %*% wp - y)^2))
  }
  list(w = w, pval = (1 + sum(rms_perm <= rms_obs)) / (1 + n_perm),
       rmse = rms_obs)
}

set.seed(20260824)
res <- parallel::mclapply(seq_len(ncol(Y)), fit_one_sample,
                          Y = Y, B = B, n_perm = n_perm, mc.cores = ncores)
fractions <- t(sapply(res, function(z) z$w))
colnames(fractions) <- colnames(B); rownames(fractions) <- colnames(Y)
pvals <- sapply(res, function(z) z$pval); names(pvals) <- colnames(Y)
rmses <- sapply(res, function(z) z$rmse)
cat("all samples done\n")

write.table(fractions, "output/m4_bulk/cibersort_fractions.tsv", sep = "\t", quote = FALSE)
write.table(data.frame(sample = names(pvals), p_value = pvals, rmse = rmses),
            "output/m4_bulk/cibersort_pvalues.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# --- Group comparison (Wilcoxon) ----------------------------------------------
res_list <- list()
for (ct in colnames(fractions)) {
  x <- fractions[grp == "Rosacea", ct]
  y2 <- fractions[grp == "HC", ct]
  wt <- wilcox.test(x, y2)
  res_list[[length(res_list) + 1]] <- data.frame(
    cell_type = ct, mean_Rosacea = mean(x), mean_HC = mean(y2),
    median_Rosacea = median(x), median_HC = median(y2),
    p_value = wt$p.value)
}
cell_res <- do.call(rbind, res_list)
cell_res$fdr <- p.adjust(cell_res$p_value, method = "BH")
cell_res <- cell_res[order(cell_res$p_value), ]
write.table(cell_res, "output/m4_bulk/cibersort_group_comparison.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(cell_res[cell_res$p_value < 0.05, c("cell_type", "mean_Rosacea", "mean_HC", "p_value", "fdr")])

# --- GSEA (Hallmark, pre-ranked by t statistic) -------------------------------
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(msigdbr))

deg_all <- read.delim("output/m4_bulk/DEG_all.tsv")
ranks <- setNames(deg_all$t, deg_all$gene)
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

hall <- msigdbr(species = "Homo sapiens", category = "H")
hall <- hall[, c("gs_name", "gene_symbol")]
set.seed(20260824)
gsea <- GSEA(ranks, TERM2GENE = hall, pvalueCutoff = 0.05,
             minGSSize = 15, maxGSSize = 500, seed = TRUE)
gsea_df <- as.data.frame(gsea)
write.table(gsea_df, "output/m4_bulk/gsea_hallmark_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
if (nrow(gsea_df) > 0) {
  print(gsea_df[order(gsea_df$NES, decreasing = TRUE), c("ID", "NES", "p.adjust")][1:12, ])
}
cat("done\n")
