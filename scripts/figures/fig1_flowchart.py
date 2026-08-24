#!/usr/bin/env python3
"""fig1_flowchart.py — Figure 1: data-source -> filtering -> integration flow chart
Every underlying data source, every filtering step, every intermediate dataset with counts.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(13, 9.5), dpi=300)
ax.set_xlim(0, 13); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc="#eef3fb", ec="#1b5a94", fs=8.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                fc=fc, ec=ec, lw=1.2))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs)

def arrow(x1, y1, x2, y2, color="#555555"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=12, lw=1.2, color=color))

GREEN, ORANGE, PINK, BLUE = "#e8f5e9", "#fff3e0", "#fce4ec", "#eef3fb"
GE, OE, PE, BE = "#2e7d32", "#e65100", "#c62828", "#1b5a94"

# --- Row 1: data sources ---
box(0.2, 8.6, 3.0, 1.1, "scRNA-seq skin\nChen et al. 2024\nGSA: HRA006167\n(131,243 cells, 21 samples)")
box(3.5, 8.6, 3.0, 1.1, "Skin 16S rRNA (V3-V4)\nJoura et al. 2024\nSRA: PRJNA1189573\n(31 skin samples)")
box(6.8, 8.6, 3.0, 1.1, "Bulk microarray\nBuhl et al. 2015\nGEO: GSE65914\n(58 arrays: 38 R, 20 HC)")
box(10.1, 8.6, 2.7, 1.1, "Circulating metabolomics\nLi 2025 (149 DAMs)\nZhang 2025 (ETR/PPR/HC)")

# --- Row 2: databases ---
box(0.2, 7.2, 3.0, 1.0, "MEBOCOST met-sensor DB\n793 pairs / 413 sensor genes / 217 metabolites\n(Oct 2025 release)")
box(3.5, 7.2, 3.0, 1.0, "SILVA v138.1 reference\n(400,300 pure-ACGT seqs)")
box(6.8, 7.2, 3.0, 1.0, "LM22 signature\nNewman et al. 2015\n(547 genes x 22 subsets)")
box(10.1, 7.2, 2.7, 1.0, "DGIdb GraphQL API\n(drug-gene interactions)")

# --- Row 3: first processing steps ---
box(0.2, 5.5, 3.0, 1.2, "M1 Sensor expression atlas\npct.exp + avg.exp per cell type\n(11 annotated cell types)\n[pending controlled-access data]", fc=GREEN, ec=GE)
arrow(1.7, 8.6, 1.7, 6.7)
box(3.5, 5.5, 3.0, 1.2, "M2a DADA2 processing\n1,660 ASVs\n(median reads/sample [X])\nIdTaxa -> genus table", fc=GREEN, ec=GE)
arrow(5.0, 8.6, 5.0, 6.7)
box(6.8, 5.5, 3.0, 1.2, "M4a RMA + limma\n1,416 DEGs (all subtypes vs HC)\nETR 1,206 / PPR 1,662 / PhR 1,655", fc=GREEN, ec=GE)
arrow(8.3, 8.6, 8.3, 6.7)
box(10.1, 5.5, 2.7, 1.2, "DGIdb queries\n27+ drug-gene\ninteractions retrieved", fc=ORANGE, ec=OE)
arrow(11.45, 8.6, 11.45, 6.7)

# --- Row 4: second processing steps ---
box(0.2, 4.0, 3.0, 1.1, "M3 Receiver-side scoring\nsensor expr x cell-type specificity\n+ permutation test (FDR)\n[pending scRNA-seq]", fc=ORANGE, ec=OE)
arrow(1.7, 5.5, 1.7, 5.1)
box(3.5, 4.0, 3.0, 1.1, "M2b Differential genera\n(published Table 2, Joura 2024)\n8 genera P<0.05;\nCutibacterium/Staph dominant n.s.", fc=GREEN, ec=GE)
arrow(5.0, 5.5, 5.0, 5.1)
box(6.8, 4.0, 3.0, 1.1, "M4b Nested-CV classifier\n(sensor genes; DEG per fold;\nLASSO/RF/SVM-RFE consensus)\nAUC = 1.000 (25/25 folds)", fc=GREEN, ec=GE)
arrow(8.3, 5.5, 8.3, 5.1)

# --- Row 5: third processing steps ---
box(0.2, 2.5, 3.0, 1.1, "M2d Metabolomics intersection\n5-ALA up 3.8x in PPR plasma;\nSCFAs not covered by LC-MS\n(supporting evidence only)", fc=GREEN, ec=GE)
arrow(3.2, 4.55, 2.6, 3.6)
box(3.5, 2.5, 3.0, 1.1, "M2c Genus -> metabolite links\n(literature, PMIDs verified)\nx MEBOCOST sensor filter\n= 7 candidate metabolites", fc=GREEN, ec=GE)
arrow(5.0, 4.0, 5.0, 3.6)
box(6.8, 2.5, 3.0, 1.1, "M4c CIBERSORT (LM22)\nM1 mac./Tfh/Treg up;\nresting mast cells down\nGSEA: IFN-g NES 3.08", fc=GREEN, ec=GE)
arrow(8.3, 4.0, 8.3, 3.6)

# --- Triangulation (wide) ---
box(0.2, 1.0, 12.6, 1.2, "MULTI-SOURCE PRIORITIZATION OF CANDIDATE AXES (evidence chain table: per-link data source + statistic)\n1. SCFA axis: depleted Faecalibacterium/Bacteroides/Prevotella/Blautia/Ruminococcus/Subdoligranulum -> acetate/propionate/butyrate -> FFAR2/FFAR3/HCAR2/SLC16A1\n2. Heme axis: Cutibacterium porphyrin degradation (excess heme in rosacea skin) -> heme -> TLR4 (myeloid cells)\n3. Supporting chains: succinate->SUCNR1; 5-ALA->TSPO2; adenosylcobalamin->CUBN/AMN",
     fc=PINK, ec=PE, fs=8)
arrow(1.7, 2.5, 3.5, 2.2)
arrow(5.0, 2.5, 5.4, 2.2)
arrow(8.3, 2.5, 7.9, 2.2)

# --- M5 ---
box(6.8, -0.2, 6.0, 1.0, "M5 Drug prioritization + docking (illustrative)\nRamatroban/setipiprant -> PTGDR2 (CRTH2): -11.2 / -11.8 kcal/mol\nTAK-242 -> TLR4 TIR (4G8A): -4.2 | doxycycline -> MMP9 (1L6J): -3.3",
     fc=ORANGE, ec=OE, fs=8)
arrow(9.8, 1.0, 9.8, 0.8)

# --- Cross-links (non-linear module connections, per Reviewer 2) ---
ax.plot([3.2, 5.2], [7.7, 7.7], ls=(0, (4, 2)), lw=1.1, color="#8e24aa")
ax.text(4.2, 7.85, "met-sensor DB", fontsize=7.5, color="#8e24aa", ha="center")
ax.plot([5.0, 6.8], [4.55, 4.55], ls=(0, (4, 2)), lw=1.1, color="#8e24aa")
ax.text(5.9, 4.7, "candidate metabolites", fontsize=7.5, color="#8e24aa", ha="center")
ax.plot([1.7, 3.5], [6.1, 6.1], ls=(0, (4, 2)), lw=1.1, color="#8e24aa")
ax.text(2.6, 6.25, "sensor-gene filter", fontsize=7.5, color="#8e24aa", ha="center")

ax.text(0.2, 9.85, "Figure 1. Data sources, filtering steps, and information flow. Each intermediate dataset is shown with counts; all intermediate files are provided as supplementary data.",
        fontsize=9, weight="bold")
plt.tight_layout()
plt.savefig("output/figures/Fig1_flowchart.png", bbox_inches="tight", dpi=300)
print("Fig1 saved")
