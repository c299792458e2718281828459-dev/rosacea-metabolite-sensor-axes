#!/usr/bin/env python3
"""fig1_flowchart.py — Figure 1: data-source -> filtering -> integration flow chart.
Four vertical chains; cross-column connectors routed through margins (no arrow
crosses any box). scRNA-seq removed (controlled access); metabolomics role shown
as explicit supporting filter.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(16, 11.2), dpi=300)
ax.set_xlim(-0.9, 16.0); ax.set_ylim(0, 11.2); ax.axis("off")

def box(x, y, w, h, text, fc="#eef3fb", ec="#1b5a94", fs=8.5, lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs)

def arrow(x1, y1, x2, y2, color="#555555", style="-|>", lw=1.3, ls="-", conn=None):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=13, lw=lw, color=color, linestyle=ls,
                                 connectionstyle=conn or "arc3,rad=0"))

def flabel(x, y, text, color="#444444", fs=7.5, rot=0):
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=color,
            rotation=rot, bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.95))

GREEN, ORANGE, PINK, BLUE, GRAY = "#e8f5e9", "#fff3e0", "#fce4ec", "#eef3fb", "#f5f5f5"
GE, OE, PE, BE = "#2e7d32", "#e65100", "#c62828", "#1b5a94"

# column x positions
c1, c2, c3, c4 = 0.4, 3.9, 7.5, 11.1
w = 3.2
c1c, c2c, c3c, c4c = c1+w/2, c2+w/2, c3+w/2, c4+w/2
# row y positions
r1, r2, r3, r4, r5 = 8.55, 6.95, 5.05, 3.15, 1.05
h1, h2, h3, h4, h5 = 1.25, 1.15, 1.35, 1.55, 1.55

# ================= column 1: microbiome chain =================
box(c1, r1, w, h1, "Skin 16S rRNA (V3-V4)\nJoura et al. 2024\nSRA: PRJNA1189573\n(31 skin samples)", fc=BLUE, ec=BE)
box(c1, r2, w, h2, "SILVA v138.1 reference\n(V3-V4 region,\n301,898 sequences)")
box(c1, r3, w, h3, "M2a DADA2 re-processing\n1,660 ASVs, 228,836 reads\nmedian 8,020 reads/sample\nReproducibility check:\nEscherichia-dominated (caveat)", fc=GREEN, ec=GE, fs=7.8)
box(c1, r4, w, h4, "M2b Genus -> metabolite links\n10 genera (8 differential +\n2 abundance-dominant)\nliterature evidence (verified\nPMIDs) + MiMeDB entry IDs", fc=PINK, ec=PE, fs=7.8)
arrow(c1c, r1, c1c, r2+h2); arrow(c1c, r2, c1c, r3+h3); arrow(c1c, r3, c1c, r4+h4)

# ================= column 2: host transcriptome chain =================
box(c2, r1, w, h1, "Bulk microarray\nBuhl et al. 2015\nGEO: GSE65914\n(58 arrays: 38 R, 20 HC)", fc=BLUE, ec=BE)
box(c2, r2, w, h2, "M4a RMA + limma\n1,416 DEGs (all vs HC)\n848 up / 568 down\nETR 1,206 / PPR 1,662 / PhR 1,655", fc=GREEN, ec=GE, fs=7.8)
box(c2, r3, w, h3, "M4c CIBERSORT (LM22)\n22 immune subsets, 58 samples\n\nM4d GSEA (Hallmark)\n36 significant gene sets\nIFN-gamma NES = 3.08", fc=GREEN, ec=GE, fs=7.8)
box(c2, r4, w, h4, "M4b Nested CV classifier\nsensor genes; feature-selection\npurpose (not disease prediction)\nAUC = 1.000, 25/25 folds\ninput: GSE65914 expression matrix", fc=GREEN, ec=GE, fs=7.8)
arrow(c2c, r1, c2c, r2+h2); arrow(c2c, r2, c2c, r3+h3)

# ================= column 3: integration chain =================
box(c3, r1, w, h1, "Curated databases\nMEBOCOST (793 pairs, 413 sensors)\nMiMeDB v2.0\nLM22 (547 x 22) + DGIdb API", fc=BLUE, ec=BE, fs=7.8)
box(c3, r2, w, h2, "M2c Sensor filter\nmetabolites with MEBOCOST-\nannotated human sensors\n(793 pairs)", fc=PINK, ec=PE, fs=7.8)
box(c3, r3, w, h3, "M2d Candidate set\n7 candidate metabolites\n(3 SCFAs, heme, succinate,\n5-ALA, adenosylcobalamin)", fc=PINK, ec=PE, fs=7.8)
box(c3, r4, w, h4, "Prioritized candidate axes\n(Table 2 evidence chain)\nsensor DEG status from M4a\nFFAR2/FFAR3/TLR4/NOD2/AHR", fc=GREEN, ec=GE, fs=7.8)
box(c3, r5, w, h5, "M5 Drug prediction + docking\nDGIdb queries\nAutoDock Vina 1.2.7\n4 drug-target pairs", fc=ORANGE, ec=OE, fs=7.8)
arrow(c3c, r1, c3c, r2+h2); arrow(c3c, r2, c3c, r3+h3); arrow(c3c, r3, c3c, r4+h4); arrow(c3c, r4, c3c, r5+h5)

# ================= column 4: published inputs + CV feature space =================
box(c4, r1, w, h1, "Published genus statistics\nJoura 2024, Table 2\n10 skin genera\n(Wilcoxon P values)", fc=BLUE, ec=BE)
box(c4, r2, w, h2, "Published metabolomics tables\nLi et al. 2025 (149 serum DAMs)\nZhang et al. 2025 (plasma,\n17 ETR / 17 PPR / 16 HC)", fc=BLUE, ec=BE, fs=7.8)
box(c4, r3, w, h3, "Supporting metabolomics filter\nintersect candidate list with\nLi/Zhang differential tables\n(5-ALA up in PPR plasma ->\nTSPO2 chain; SCFA absence\nuninformative)", fc=ORANGE, ec=OE, fs=7.6)
arrow(c4c, r2, c4c, r3+h3)
# filter -> candidate set (short horizontal)
arrow(c4, r3+h3/2, c3+w, r3+h3/2)
flabel((c4+c3+w)/2, r3+h3/2+0.42, "5-ALA")

# ================= routed connectors (margin channels only; no box crossings) =================
# 1) published genus stats -> M2b (over the top margin, down the left margin)
arrow(c4c, r1+h1, c4c, 10.75)                       # up from pub-stats top
arrow(c4c, 10.75, -0.45, 10.75)                     # across the top margin
arrow(-0.45, 10.75, -0.45, r4+h4/2)                 # down the left margin
arrow(-0.45, r4+h4/2, c1, r4+h4/2)                  # into M2b left edge
flabel(-0.45, 7.7, "genus differential status (P values)", rot=90, fs=7.2)

# 2) M2b -> M2c (down, across the r4-r5 strip, up the right margin, over the r1-r2 strip, into M2c)
arrow(c1c, r4, c1c, 2.85)                           # down from M2b bottom
arrow(c1c, 2.85, 15.35, 2.85)                       # across the free strip below row 4
arrow(15.35, 2.85, 15.35, 8.42)                     # up the right margin (right of all c4 boxes)
arrow(15.35, 8.42, c3+w, 8.42)                      # across the free strip above row 2
arrow(c3+w, 8.42, c3+w, r2+h2)                      # down into M2c right edge
flabel(15.35, 5.7, "metabolite list\n(49 evidence links)", rot=90, fs=7.2)

# 3) M2c -> M4b (sensor gene list; via the gap between columns 3 and 4, then the r3-r4 strip)
arrow(c3+w, r2+0.3, 10.95, 4.82)                    # down through the c3-c4 gap
arrow(10.95, 4.82, c2c, 4.82)                       # left across the r3-r4 strip
arrow(c2c, 4.82, c2c, r4+h4)                        # into M4b top edge
flabel(8.2, 4.97, "sensor gene list (413)", fs=6.8)

# M4a -> M4c/M4d is vertical (drawn above); M4a also feeds axes evidence via Table 2 text (no arrow)

# caption
ax.text(-0.45, 0.35,
        "Figure 1. Data sources, filtering steps, and information flow of the study. Each underlying data source is shown at the top;\n"
        "each processing step shows the size of the intermediate dataset it produces; labels on routed connectors state what is passed.\n"
        "The scRNA-seq module was removed because the controlled-access dataset (GSA: HRA006167) was not granted; cell-type localization\n"
        "is cited from published results. The published metabolomics tables serve as a supporting filter on the candidate set.",
        fontsize=8.5, va="bottom", ha="left")

plt.tight_layout()
plt.savefig("output/figures/Fig1_flowchart.png", dpi=300, bbox_inches="tight")
print("Fig1 saved")
