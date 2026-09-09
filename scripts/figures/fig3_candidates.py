#!/usr/bin/env python3
"""fig3_candidates.py — Figure 3: candidate microbial metabolites -> sensors, per producing genus
Data: output/m2_microbiome/genus_metabolite_sensor_evidence.tsv (all PMIDs verified)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv

rows = list(csv.DictReader(open("output/m2_microbiome/genus_metabolite_sensor_evidence.tsv"), delimiter="\t"))
rows = [r for r in rows if r["metabolite_sensors_MEBOCOST"]]

# direction colors for genera
dir_color = {"up": "#c62828", "down": "#1565c0", "n.s.": "#757575"}
genus_dir = {}
for r in rows:
    genus_dir[r["genus"]] = r["genus_direction_Joura2024"]

# structure: metabolite -> (sensors, genera)
mets = {}
for r in rows:
    m = r["metabolite"]
    mets.setdefault(m, {"sensors": set(), "genera": {}})
    mets[m]["sensors"].add(r["metabolite_sensors_MEBOCOST"])
    mets[m]["genera"].setdefault(r["genus"], r["genus_direction_Joura2024"])

order = sorted(mets, key=lambda m: -len(mets[m]["genera"]))

# metabolomics-filter status per candidate (per R3 comment)
# ▲ concordant circulating signal; ○ not measurable by untargeted LC-MS (volatile SCFA); □ measurable, no concordant signal
filter_status = {
    "Acetic acid": ("\u25cb", "not measurable (volatile SCFA)"),
    "Propionic acid": ("\u25cb", "not measurable (volatile SCFA)"),
    "Butyric acid": ("\u25cb", "not measurable (volatile SCFA)"),
    "Succinic acid": ("\u25b2", "up in rosacea serum (Li 2025, log2FC = 0.50, P = 1e-4)"),
    "5-Aminolevulinic acid": ("\u25b2", "up in PPR plasma (Zhang 2025)"),
    "Heme": ("\u25a1", "no concordant circulating signal"),
    "Adenosylcobalamin (vitamin B12)": ("\u25a1", "no concordant circulating signal"),
}

fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=600)
y = 0
yt = {}
for m in order:
    yt[m] = y
    sym = filter_status.get(m, ("", ""))[0]
    ax.text(-0.55, y, f"{sym}  {m}", ha="right", va="center", fontsize=7.5)
    gen = sorted(mets[m]["genera"], key=lambda g: -genus_dir[g].count("d"))
    x = 0
    for g in gen:
        d = mets[m]["genera"][g]
        ax.scatter(x, y, s=110, color=dir_color[d], zorder=3,
                   edgecolors="white", linewidths=0.5)
        # label's right end centered below the circle (per reviewer)
        ax.text(x + 0.06, y - 0.18, g, ha="right", va="top", fontsize=6.5, rotation=25)
        x += 1.15
    sens = sorted(mets[m]["sensors"])
    import textwrap
    sens_label = "; ".join(sens)
    wrapped = textwrap.wrap(sens_label, width=46)
    xs = x + 0.9
    n = len(wrapped)
    for i, line in enumerate(wrapped):
        off = (i - (n - 1) / 2) * 0.42
        ax.text(xs, y + off, line, fontsize=7, color="#4a148c", va="center", ha="left")
    y += 1.45

ax.set_ylim(-0.9, y - 0.55)
ax.set_xlim(-5.5, 20)
ax.axis("off")
# legend
for lab, c in [("genus up in rosacea", "#c62828"), ("genus down in rosacea", "#1565c0"),
               ("dominant, not significant", "#757575")]:
    pass
import matplotlib.lines as mlines
h1 = mlines.Line2D([], [], marker="o", ls="", color="#c62828", label="genus enriched in rosacea")
h2 = mlines.Line2D([], [], marker="o", ls="", color="#1565c0", label="genus depleted in rosacea")
h3 = mlines.Line2D([], [], marker="o", ls="", color="#757575", label="dominant, not significant (Joura 2024)")
h4 = mlines.Line2D([], [], marker="s", ls="", color="#4a148c", label="sensor (MEBOCOST DB)")
ax.legend(handles=[h1, h2, h3, h4], loc="upper right", fontsize=8, frameon=False)
ax.set_title("Figure 2. Candidate microbial metabolites: producing genera and annotated human sensors", fontsize=10, loc="left")
plt.tight_layout()
plt.savefig("output/figures/Fig3_candidate_metabolites.png", bbox_inches="tight", dpi=300)
print("Fig3 saved")
