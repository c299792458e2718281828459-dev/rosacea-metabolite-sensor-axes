#!/usr/bin/env python3
"""fig1_flowchart.py — Figure 1: data-source -> filtering -> integration flow chart.
Designed AT print size: 7.2 in wide (Frontiers max 180 mm ~ 7.09 in); all fonts
>= 7 pt at final size, so nothing shrinks below legibility when printed.
Four vertical chains; cross-column connectors routed through margins (no arrow
crosses any box).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(7.2, 6.6), dpi=600)
ax.set_xlim(-0.75, 8.15); ax.set_ylim(0, 7.4); ax.axis("off")

def box(x, y, w, h, text, fc="#eef3fb", ec="#1b5a94", fs=7.5, lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                fc=fc, ec=ec, lw=lw))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, linespacing=1.15)

def arrow(x1, y1, x2, y2, color="#555555", lw=0.9):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=9, lw=lw, color=color))

def flabel(x, y, text, color="#444444", fs=7.5, rot=0):
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=color,
            rotation=rot, bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.95))

GREEN, ORANGE, PINK, BLUE = "#e8f5e9", "#fff3e0", "#fce4ec", "#eef3fb"
GE, OE, PE, BE = "#2e7d32", "#e65100", "#c62828", "#1b5a94"

# columns (4), rows (5)
c1, c2, c3, c4 = 0.35, 2.15, 3.95, 5.75
w = 1.62
c1c, c2c, c3c, c4c = c1+w/2, c2+w/2, c3+w/2, c4+w/2
r1, r2, r3, r4, r5 = 6.15, 4.95, 3.6, 2.35, 1.0
h1, h2, h3, h4, h5 = 0.92, 0.92, 1.0, 1.0, 1.0

# ---- column 1: microbiome ----
box(c1, r1, w, h1, "Skin 16S rRNA (V3-V4)\nPRJNA1189573\n(31 skin samples)", fc=BLUE, ec=BE)
box(c1, r2, w, h2, "SILVA v138.1\n(V3-V4 region,\n301,898 seqs)")
box(c1, r3, w, h3, "M2a DADA2\n1,660 ASVs\n228,836 reads\n(median 8,020/sample)", fc=GREEN, ec=GE)
box(c1, r4, w, h4, "M2b Genus to\nmetabolite links\n(10 genera; lit.\nPMIDs + MiMeDB)", fc=PINK, ec=PE)
arrow(c1c, r1, c1c, r2+h2); arrow(c1c, r2, c1c, r3+h3); arrow(c1c, r3, c1c, r4+h4)

# ---- column 2: host transcriptome ----
box(c2, r1, w, h1, "Bulk microarray\nGSE65914\n(58 arrays: 38 R, 20 HC)", fc=BLUE, ec=BE)
box(c2, r2, w, h2, "M4a RMA + limma\n1,416 DEGs\n(848 up / 568 down)", fc=GREEN, ec=GE)
box(c2, r3, w, h3, "M4c CIBERSORT (LM22)\n(22 immune subsets)\nM4d GSEA (Hallmark)\n(36 significant sets)", fc=GREEN, ec=GE)
box(c2, r4, w, h4, "M4b Nested CV\nsensor genes\n(feature selection,\nnot prediction)", fc=GREEN, ec=GE)
arrow(c2c, r1, c2c, r2+h2); arrow(c2c, r2, c2c, r3+h3)

# ---- column 3: integration ----
box(c3, r1, w, h1, "Curated databases\nMEBOCOST (793 pairs)\nMiMeDB v2, LM22,\nDGIdb API", fc=BLUE, ec=BE)
box(c3, r2, w, h2, "M2c Sensor filter\n(metabolites with\nMEBOCOST sensors)", fc=PINK, ec=PE)
box(c3, r3, w, h3, "M2d Candidate set\n(7 metabolites:\n3 SCFAs, heme,\nsuccinate, 5-ALA, B12)", fc=PINK, ec=PE)
box(c3, r4, w, h4, "Prioritized axes\n(Table 2 evidence;\nsensor DEG status\nFFAR2/3, TLR4,\nNOD2, AHR)", fc=GREEN, ec=GE)
box(c3, r5, w, h5, "M5 Drug prediction\n+ docking\n(DGIdb; Vina,\n4 drug-target pairs)", fc=ORANGE, ec=OE)
arrow(c3c, r1, c3c, r2+h2); arrow(c3c, r2, c3c, r3+h3); arrow(c3c, r3, c3c, r4+h4); arrow(c3c, r4, c3c, r5+h5)

# ---- column 4: published inputs ----
box(c4, r1, w, h1, "Published genus stats\nJoura 2024, Table 2\n(10 skin genera,\nWilcoxon P values)", fc=BLUE, ec=BE)
box(c4, r2, w, h2, "Published metabolomics\nLi 2025 (serum),\nZhang 2025 (plasma)", fc=BLUE, ec=BE)
box(c4, r3, w, h3, "Supporting filter\n(candidate x\ncirculating tables)", fc=ORANGE, ec=OE)
arrow(c4c, r2, c4c, r3+h3)
arrow(c4, r3+h3/2, c3+w, r3+h3/2)
flabel((c4+c3+w)/2, r3+h3/2+0.33, "5-ALA;\nsuccinate up", fs=6.8)

# ---- routed connectors (margins; none crosses a box) ----
# published genus stats -> M2b (top margin, left margin)
arrow(c4c, r1+h1, c4c, 7.05)
arrow(c4c, 7.05, -0.5, 7.05)
arrow(-0.5, 7.05, -0.5, r4+h4/2)
arrow(-0.5, r4+h4/2, c1, r4+h4/2)
flabel(-0.5, 5.4, "genus differential\nstatus (P values)", rot=90, fs=7.5)

# M2b -> M2c (down, bottom strip, right margin, top strip)
arrow(c1c, r4, c1c, 0.62)
arrow(c1c, 0.62, 7.75, 0.62)
arrow(7.75, 0.62, 7.75, 5.7)
arrow(7.75, 5.7, c3+w, 5.7)
arrow(c3+w, 5.7, c3+w, r2+h2)
flabel(7.75, 3.2, "metabolite list\n(49 links)", rot=90, fs=7.5)

# M2c -> M4b (sensor gene list; through c3-c4 gap and r3-r4 strip)
arrow(c3+w, r2+0.25, 6.15, 3.05)
arrow(6.15, 3.05, c2c, 3.05)
arrow(c2c, 3.05, c2c, r4+h4)
flabel(5.55, 3.28, "sensor gene list (413)", fs=6.8)

plt.tight_layout()
plt.savefig("output/figures/Fig1_flowchart.png", dpi=600, bbox_inches="tight")
print("Fig1 saved (print-size layout)")
