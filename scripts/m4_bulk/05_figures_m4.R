# 05_figures_m4.R — M4 figures from real results (volcano, nested-CV AUC, CIBERSORT, GSEA)
library(ggplot2)
library(dplyr)

dir.create("output/figures", showWarnings = FALSE)

# --- Figure: Volcano of DEGs ------------------------------------------------
deg <- read.delim("output/m4_bulk/DEG_all.tsv")
deg$sig <- "NS"
deg$sig[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1] <- "DEG"
label_genes <- c("PTGDS", "CCL19", "TLR2", "MMP9", "S100A8", "CXCL8", "STAT1",
                 "FFAR2", "FFAR3", "NOD2", "AHR", "PTGDR")
p1 <- ggplot(deg, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("NS" = "grey70", "DEG" = "#d95f02")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = 2) +
  ggrepel::geom_text_repel(data = subset(deg, gene %in% label_genes),
                           aes(label = gene), size = 3, color = "black",
                           max.overlaps = 20) +
  labs(x = "log2 fold change (rosacea vs control)",
       y = "-log10 adjusted P", title = "GSE65914: 1,416 DEGs") +
  theme_bw() + theme(legend.position = "none")
ggsave("output/figures/Fig_volcano_DEG.png", p1, width = 6, height = 5, dpi = 300)

# --- Figure: nested CV AUC (25 folds) + probability separation ---------------
cv <- read.delim("output/m4_bulk/nested_cv_folds.tsv")
p2a <- ggplot(cv, aes(x = factor(rep), y = auc)) +
  geom_jitter(width = 0.1, color = "#1b9e77", size = 2) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  labs(x = "CV repeat", y = "held-out fold AUC",
       subtitle = "Nested CV AUC = 1.000 (25/25 folds)") +
  theme_bw()
prd <- rbind(data.frame(class = "Rosacea", prob = cv$pr_mean_Rosacea),
             data.frame(class = "Control", prob = cv$pr_mean_HC))
p2b <- ggplot(prd, aes(class, prob, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Rosacea" = "#d95f02", "Control" = "#7570b3")) +
  labs(y = "predicted probability (Rosacea)") +
  theme_bw() + theme(legend.position = "none")
ggsave("output/figures/Fig_nestedCV.png",
       gridExtra::grid.arrange(p2a, p2b, nrow = 1), width = 9, height = 4, dpi = 300)

# --- Figure: CIBERSORT group comparison (single-panel horizontal bars) ---------
ci <- read.delim("output/m4_bulk/cibersort_group_comparison.tsv")
ci$diff <- ci$mean_Rosacea - ci$mean_HC
ci$cell_type <- factor(ci$cell_type, levels = ci$cell_type[order(ci$diff)])
ci$label <- paste0(ci$cell_type, ifelse(ci$p_value < 0.05, "*", ""))
p3 <- ggplot(ci, aes(x = diff, y = cell_type)) +
  geom_col(aes(fill = diff > 0), width = 0.7) +
  geom_text(aes(x = ifelse(diff > 0, diff + 0.005, diff - 0.005),
                label = ifelse(p_value < 0.05, "*", "")), size = 4) +
  scale_fill_manual(values = c("#7570b3", "#d95f02"), guide = "none") +
  labs(x = "mean fraction difference (Rosacea - Control)",
       y = "", title = "CIBERSORT (LM22), GSE65914",
       caption = "* Wilcoxon P < 0.05") +
  theme_bw() + theme(axis.text.y = element_text(size = 8))
ggsave("output/figures/Fig_cibersort.png", p3, width = 7, height = 7, dpi = 300)

# --- Figure: GSEA top pathways ------------------------------------------------
gsea <- read.delim("output/m4_bulk/gsea_hallmark_results.tsv")
gsea <- gsea[order(gsea$NES, decreasing = TRUE), ][1:12, ]
gsea$ID <- gsub("HALLMARK_", "", gsea$ID)
gsea$ID <- factor(gsea$ID, levels = rev(gsea$ID))
p4 <- ggplot(gsea, aes(NES, ID, fill = -log10(p.adjust))) +
  geom_col() +
  scale_fill_gradient(low = "#fdd49e", high = "#7f2704") +
  labs(x = "normalized enrichment score", y = "", fill = "-log10 adj.P",
       title = "GSEA Hallmark (top 12 activated)") +
  theme_bw()
ggsave("output/figures/Fig_gsea.png", p4, width = 7, height = 5, dpi = 300)
cat("figures done\n")
