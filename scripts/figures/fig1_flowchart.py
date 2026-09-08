#!/usr/bin/env python3
"""fig1_flowchart.py — Figure 1: data-source -> filtering -> integration flow chart.
Layout per R2 round-2 comments: no scRNA-seq module (dropped, controlled access);
circulating metabolomics role shown explicitly as supporting filter; readable filter labels.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(13.5, 10), dpi=300)
ax.set_xlim(0, 13.5); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc="#eef3fb", ec="#1b5a94", fs=8.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                fc=fc, ec=ec, lw=1.2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs)

def arrow(x1, y1, x2, y2, color="#555555", style="-|>", lw=1.2, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=12, lw=lw, color=color, linestyle=ls))

def flabel(x, y, text, color="#444444", fs=7.5):
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=color,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9))

GREEN, ORANGE, PINK, BLUE, GRAY = "#e8f5e9", "#fff3e0", "#fce4ec", "#eef3fb", "#f5f5f5"
GE, OE, PE, BE = "#2e7d32", "#e65100", "#c62828", "#1b5a94"

# ---- Row 1: underlying data sources (top) ----
box(0.3, 8.55, 3.4, 1.15, "Skin 16S rRNA (V3-V4)\nJoura et al. 2024\nSRA: PRJNA1189573\n(31 skin samples)", fc=BLUE, ec=BE)
box(4.0, 8.55, 3.4, 1.15, "Bulk microarray\nBuhl et al. 2015\nGEO: GSE65914\n(58 arrays: 38 R, 20 HC)", fc=BLUE, ec=BE)
box(7.7, 8.55, 3.4, 1.15, "Published metabolomics tables\nLi et al. 2025 (149 serum DAMs)\nZhang et al. 2025\n(plasma, 17 ETR / 17 PPR / 16 HC)", fc=BLUE, ec=BE)
box(11.4, 8.55, 1.9, 1.15, "Published genus stats\nJoura 2024 Table 2\n(10 skin genera,\nWilcoxon P values)", fc=BLUE, ec=BE)

# ---- Row 2: databases / references ----
box(0.3, 7.0, 3.4, 1.05, "SILVA v138.1\nV3-V4 region reference\n(301,898 sequences)")
box(4.0, 7.0, 3.4, 1.05, "MEBOCOST met-sensor DB\n793 pairs / 413 sensor genes\n(release Oct 2025)")
box(7.7, 7.0, 3.4, 1.05, "MiMeDB v2.0\n(29,296 metabolites,\n2,648 microbes)")
box(11.4, 7.0, 1.9, 1.05, "LM22 signature\n(547 genes x 22)\nDGIdb API")

# ---- Row 3: processing modules ----
box(0.3, 5.15, 3.4, 1.3, "M2a DADA2 re-processing\n1,660 ASVs, 228,836 reads\n(median 8,020 / sample)\nIdTaxa genus assignment", fc=GREEN, ec=GE)
arrow(2.0, 8.55, 2.0, 6.45)
box(4.0, 5.15, 3.4, 1.3, "M4a RMA + limma\n1,416 DEGs (all vs HC)\n848 up / 568 down\nETR 1,206 / PPR 1,662 / PhR 1,655", fc=GREEN, ec=GE)
arrow(5.7, 8.55, 5.7, 6.45)
box(7.7, 5.15, 3.4, 1.3, "Supporting metabolomics filter\nintersect candidate list with\nLi/Zhang differential tables\n(5-ALA up in PPR plasma -> TSPO2)", fc=ORANGE, ec=OE)
arrow(9.4, 8.55, 9.4, 6.45)
box(11.4, 5.15, 1.9, 1.3, "Reproducibility check\nof deposited reads\n(Escherichia-dominated;\ncaveat reported)", fc=GRAY, ec="#666666")
arrow(11.4+0.95, 8.55, 11.4+0.95, 6.45)

# ---- Row 4: integration (M2b/M3) ----
box(0.3, 3.15, 3.4, 1.5, "M2b Genus -> metabolite links\n8 differential + 2 dominant genera\nliterature evidence (verified PMIDs)\n+ MiMeDB entry IDs", fc=PINK, ec=PE)
arrow(2.0, 5.15, 2.0, 4.65)
flabel(2.5, 4.9, "genus table")
box(4.0, 3.15, 3.4, 1.5, "M2c Sensor filter\nrestrict to metabolites with\nMEBOCOST-annotated human sensors\n(793 pairs)", fc=PINK, ec=PE)
arrow(5.7, 5.15, 5.7, 4.65)
flabel(6.2, 4.9, "metabolite list")
box(7.7, 3.15, 3.4, 1.5, "M2d Candidate set\n7 candidate metabolites\n(3 SCFAs, heme, succinate,\n5-ALA, adenosylcobalamin)", fc=PINK, ec=PE)
arrow(8.05, 7.0, 4.35, 3.9, color="#666666");  # MiMeDB -> M2b
flabel(5.4, 6.0, "MiMeDB entries")
arrow(9.4, 5.15, 9.4, 4.65)
flabel(9.9, 4.9, "5-ALA filter")
box(11.4, 3.15, 1.9, 1.5, "Prioritized\ncandidate axes\n(Table 2 evidence\nchain)", fc=GREEN, ec=GE)
arrow(9.4+1.6, 3.9, 11.4, 3.9)

# ---- Row 5: downstream analyses (M4b-M5) ----
box(0.3, 1.15, 3.4, 1.5, "M4b Nested CV classifier\n(sensor genes; feature\nselection purpose)\nAUC = 1.000, 25/25 folds", fc=GREEN, ec=GE)
arrow(4.0+1.7, 3.15, 2.0, 2.65, color="#666666")
flabel(3.3, 3.0, "sensor gene list")
box(4.0, 1.15, 3.4, 1.5, "M4c CIBERSORT (LM22)\nimmune composition\n+ M4d GSEA (Hallmark)\nIFN-gamma NES = 3.08", fc=GREEN, ec=GE)
arrow(5.7, 5.15, 5.7, 2.65, color="#666666")
flabel(6.5, 3.0, "DEGs")
box(7.7, 1.15, 3.4, 1.5, "M5 Drug prediction\nDGIdb queries +\nAutoDock Vina docking\n(4 drug-target pairs)", fc=ORANGE, ec=OE)
arrow(7.7+1.7, 3.15, 9.4, 2.65, color="#666666")
flabel(9.2, 3.0, "targets")

# caption
ax.text(0.3, 0.25,
        "Figure 1. Data sources, filtering steps, and information flow. Each underlying data source is shown at the top; each\n"
        "processing step shows the size of the intermediate dataset it produces; dashed gray arrows mark cross-module dependencies.\n"
        "The scRNA-seq module was removed because the controlled-access dataset (GSA: HRA006167) was not granted; cell-type\n"
        "localization is cited from published results. Circulating metabolomics tables serve as a supporting filter for the candidate set.",
        fontsize=8.5, va="bottom", ha="left")

plt.tight_layout()
plt.savefig("output/figures/Fig1_flowchart.png", dpi=300, bbox_inches="tight")
print("Fig1 saved")
