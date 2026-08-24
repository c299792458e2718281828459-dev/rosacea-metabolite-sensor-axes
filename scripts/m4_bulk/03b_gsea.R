# 03b_gsea.R — GSEA (Hallmark) with local GMT, pre-ranked by t statistic
suppressPackageStartupMessages(library(clusterProfiler))

deg_all <- read.delim("output/m4_bulk/DEG_all.tsv")
ranks <- setNames(deg_all$t, deg_all$gene)
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

gmt <- clusterProfiler::read.gmt("data/hallmark.gmt")
set.seed(20260824)
gsea <- GSEA(ranks, TERM2GENE = gmt, pvalueCutoff = 0.05,
             minGSSize = 15, maxGSSize = 500, seed = TRUE)
gsea_df <- as.data.frame(gsea)
write.table(gsea_df, "output/m4_bulk/gsea_hallmark_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("n enriched sets (p.adjust<0.05):", sum(gsea_df$p.adjust < 0.05), "of", nrow(gsea_df), "\n")
if (nrow(gsea_df) > 0) {
  print(gsea_df[order(gsea_df$NES, decreasing = TRUE),
                c("ID", "NES", "p.adjust")][1:15, ])
}
cat("done\n")
