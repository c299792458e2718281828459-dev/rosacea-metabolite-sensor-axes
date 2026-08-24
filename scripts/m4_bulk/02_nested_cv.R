# 02_nested_cv.R — Nested cross-validation for sensor-gene signature
# Design (per reviewer requirements):
#  - Outer CV: repeated stratified 5-fold (5 repeats x 5 folds)
#  - INSIDE each outer training fold: limma DEG (rosacea vs HC) recomputed on training data only
#  - Feature space: genes with metabolite-sensor annotation in MEBOCOST DB (413 genes), fixed a priori
#  - Inner selection: LASSO (glmnet), random forest (Gini top quartile), SVM-RFE; consensus (>=2 votes)
#  - Final model: logistic regression on consensus genes, fit on training fold
#  - Performance: AUC on held-out test fold -> nested CV AUC (mean +- SD)
# Output: output/m4_bulk/nested_cv_results.tsv, feature_stability.tsv

library(limma)
library(glmnet)
library(randomForest)
library(e1071)
library(pROC)
library(parallel)

set.seed(20260824)

# --- Load expression + phenotype --------------------------------------------
exprs_sym <- as.matrix(read.delim("output/m4_bulk/expr_matrix.tsv", row.names = 1, check.names = FALSE))
pheno <- read.delim("data/gse65914/phenotype.tsv", stringsAsFactors = FALSE)
pheno$gsm <- sub('"', '', pheno$gsm)
pheno <- pheno[match(colnames(exprs_sym), pheno$gsm), ]
y <- factor(ifelse(pheno$subtype == "HC", "HC", "Rosacea"), levels = c("HC", "Rosacea"))
stopifnot(length(y) == ncol(exprs_sym))

# --- Candidate feature space (argv[1]: sensor | all) -------------------------
feature_space <- "sensor"
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  feature_space <- commandArgs(trailingOnly = TRUE)[1]
}
if (feature_space == "sensor") {
  sensor_db <- read.delim("data/MEBOCOST-main/data/mebocost_db/human/human_met_sensor_update_Oct21_2025.tsv", check.names = FALSE)
  feature_genes <- intersect(unique(sensor_db$Gene_name), rownames(exprs_sym))
} else {
  feature_genes <- rownames(exprs_sym)
}
cat("feature space:", feature_space, "| n genes:", length(feature_genes), "\n")
X <- t(exprs_sym[feature_genes, , drop = FALSE])   # samples x features (unscaled; scaled per fold inside CV)

# --- Helpers -----------------------------------------------------------------
# Manual stratified k-fold partition (avoids caret dependency)
stratified_folds <- function(y, k) {
  y <- as.character(y)
  groups <- split(seq_along(y), y)
  folds <- vector("list", k)
  for (g in names(groups)) {
    idx <- sample(groups[[g]])
    for (i in seq_along(idx)) folds[[((i - 1) %% k) + 1]] <- c(folds[[((i - 1) %% k) + 1]], idx[i])
  }
  folds
}

run_deg <- function(idx_train) {
  d <- model.matrix(~ y[idx_train])
  fit <- lmFit(exprs_sym[feature_genes, idx_train], d)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  rownames(tt)[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1]
}

select_genes <- function(Xtr, ytr, deg_genes) {
  if (length(deg_genes) < 5) return(character(0))
  Xs <- Xtr[, deg_genes, drop = FALSE]
  # LASSO
  lasso_sel <- tryCatch({
    cvf <- cv.glmnet(Xs, as.numeric(ytr) - 1, family = "binomial", alpha = 1, nfolds = 5)
    co <- coef(cvf, s = "lambda.1se")[-1, , drop = FALSE]
    colnames(Xs)[co != 0]
  }, error = function(e) character(0))
  # Random forest (Gini top quartile)
  rf_sel <- tryCatch({
    rf <- randomForest(Xs, ytr, ntree = 500, importance = TRUE)
    imp <- importance(rf, type = 2)[, 1]
    thr <- quantile(imp, 0.75)
    colnames(Xs)[imp >= thr]
  }, error = function(e) character(0))
  # SVM-RFE
  svm_sel <- tryCatch({
    feats <- colnames(Xs)
    keep <- feats
    while (length(keep) > 20) {
      svmfit <- svm(Xs[, keep, drop = FALSE], ytr, kernel = "linear", scale = FALSE)
      w <- t(svmfit$coefs) %*% svmfit$SV
      w <- w^2
      ranking <- order(w, decreasing = TRUE)
      keep <- keep[ranking][1:max(10, floor(0.7 * length(keep)))]
    }
    keep
  }, error = function(e) character(0))
  list(lasso = lasso_sel, rf = rf_sel, svm = svm_sel)
}

# --- Nested CV ---------------------------------------------------------------
n_repeats <- 5
n_folds <- 5
aucs <- numeric(0)
feature_votes <- setNames(rep(0L, length(feature_genes)), feature_genes)
fold_log <- list()

for (rep in seq_len(n_repeats)) {
  folds <- stratified_folds(y, n_folds)
  for (k in seq_len(n_folds)) {
    idx_test <- folds[[k]]
    idx_train <- setdiff(seq_along(y), idx_test)
    # Per-fold scaling: center/scale computed on training fold only
    mu <- colMeans(X[idx_train, , drop = FALSE])
    sdv <- apply(X[idx_train, , drop = FALSE], 2, sd)
    sdv[sdv < 1e-6] <- 1
    Xtr_all <- scale(X[idx_train, , drop = FALSE], center = mu, scale = sdv)
    Xte_all <- scale(X[idx_test, , drop = FALSE], center = mu, scale = sdv)
    deg_genes <- run_deg(idx_train)
    sel <- select_genes(Xtr_all, y[idx_train], deg_genes)
    votes <- table(unlist(sel))
    consensus <- names(votes)[votes >= 2]
    if (length(consensus) > 0) {
      feature_votes[consensus] <- feature_votes[consensus] + 1L
    }
    auc_val <- NA_real_
    pr_mean_R <- NA_real_
    pr_mean_H <- NA_real_
    if (length(consensus) >= 2) {
      Xtr <- as.data.frame(Xtr_all[, consensus, drop = FALSE])
      Xte <- as.data.frame(Xte_all[, consensus, drop = FALSE])
      Xtr$yy <- as.numeric(y[idx_train]) - 1
      m <- glm(yy ~ ., data = Xtr, family = binomial)
      pr <- predict(m, newdata = Xte, type = "response")
      if (length(unique(y[idx_test])) == 2) {
        auc_val <- as.numeric(auc(roc(y[idx_test], pr, quiet = TRUE)))
      } else {
        auc_val <- mean((pr[y[idx_test] == "Rosacea"] > pr[y[idx_test] == "HC"]))
      }
      pr_mean_R <- mean(pr[y[idx_test] == "Rosacea"])
      pr_mean_H <- mean(pr[y[idx_test] == "HC"])
    }
    aucs <- c(aucs, auc_val)
    fold_log[[length(fold_log) + 1]] <- data.frame(
      rep = rep, fold = k, n_test = length(idx_test),
      n_deg = length(deg_genes), n_consensus = length(consensus),
      auc = auc_val, pr_mean_Rosacea = pr_mean_R, pr_mean_HC = pr_mean_H,
      genes = paste(consensus, collapse = ";"))
    cat(sprintf("rep %d fold %d | DEG %d | consensus %d | AUC %.3f | prR %.3f prH %.3f\n",
                rep, k, length(deg_genes), length(consensus), auc_val, pr_mean_R, pr_mean_H))
  }
}

res <- do.call(rbind, fold_log)
write.table(res, paste0("output/m4_bulk/nested_cv_folds_", feature_space, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stab <- data.frame(gene = names(feature_votes), selected_in_folds = feature_votes)
stab <- stab[stab$selected_in_folds > 0, ]
stab <- stab[order(-stab$selected_in_folds), ]
write.table(stab, paste0("output/m4_bulk/feature_stability_", feature_space, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nNested CV AUC: mean =", round(mean(aucs, na.rm = TRUE), 3),
    "SD =", round(sd(aucs, na.rm = TRUE), 3),
    "over", sum(!is.na(aucs)), "evaluable folds\n")
cat("Top stable features:\n")
print(head(stab, 15))
cat("done\n")
