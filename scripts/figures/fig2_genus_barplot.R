# fig2_genus_barplot.R — Figure 2: skin genus composition (31 samples) + Joura comparison
library(ggplot2)
library(dplyr)
library(tidyr)

dir.create("output/figures", showWarnings = FALSE)
grel <- read.delim("output/m2_microbiome/genus_relabund.tsv", row.names = 1, check.names = FALSE)
grel <- grel * 100  # percentages

# top 12 genera, others pooled
gm <- sort(apply(grel, 2, median), decreasing = TRUE)
top <- names(gm)[1:12]
long <- as.data.frame(grel[, top, drop = FALSE]) %>%
  tibble::rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "genus", values_to = "pct")
long$genus <- factor(long$genus, levels = rev(top))
p <- ggplot(long, aes(sample, pct, fill = genus)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
                               "#e6ab02","#a6761d","#666666","#a6cee3","#b2df8a",
                               "#fb9a99","#cab2d6")) +
  labs(x = "", y = "relative abundance (%)",
       title = "Skin microbiota (16S V3-V4, n = 31 samples, PRJNA1189573)",
       fill = "genus") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.text = element_text(size = 7))
ggsave("output/figures/Fig2_genus_barplot.png", p, width = 11, height = 5.5, dpi = 300)

# comparison panel: our median vs Joura published
pub <- data.frame(
  genus = c("Cutibacterium","Neisseria","Staphylococcus","Corynebacterium",
            "Bacteroides","Faecalibacterium","Prevotella","Blautia",
            "Ruminococcus","Subdoligranulum"),
  published_rosacea = c(35.14,18.65,13.97,10.56,2.94,1.28,0.47,0.83,0.69,0.28))
pub$our_median <- sapply(as.character(pub$genus), function(g) ifelse(g %in% names(gm), gm[g], NA))
publ <- pivot_longer(pub, -genus, names_to = "source", values_to = "pct")
p2 <- ggplot(publ, aes(genus, pct, fill = source)) +
  geom_col(position = position_dodge()) +
  coord_flip() +
  scale_fill_manual(values = c("our_median" = "#1b9e77", "published_rosacea" = "#7570b3"),
                    labels = c("our re-analysis (median)", "published, rosacea (Joura 2024)")) +
  labs(x = "", y = "relative abundance (%)") +
  theme_bw() + theme(legend.position = "top", legend.title = element_blank())
ggsave("output/figures/Fig2b_joura_comparison.png", p2, width = 6.5, height = 4.5, dpi = 300)
cat("Fig2 saved\n")
