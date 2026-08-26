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

fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
y = 0
yt = {}
for m in order:
    yt[m] = y
    ax.text(-0.4, y, m, ha="right", va="center", fontsize=9)
    gen = sorted(mets[m]["genera"], key=lambda g: -genus_dir[g].count("d"))
    x = 0
    for g in gen:
        d = mets[m]["genera"][g]
        ax.scatter(x, y, s=120, color=dir_color[d], zorder=3,
                   edgecolors="white", linewidths=0.5)
        ax.text(x, y - 0.32, g, ha="center", va="top", fontsize=6, rotation=25)
        x += 1
    sens = sorted(mets[m]["sensors"])
    import textwrap
    sens_label = "; ".join(sens)
    wrapped = textwrap.wrap(sens_label, width=44)
    xs = x + 0.9
    n = len(wrapped)
    for i, line in enumerate(wrapped):
        off = (i - (n - 1) / 2) * 0.42
        ax.text(xs, y + off, line, fontsize=7.5, color="#4a148c", va="center", ha="left")
    y += 1

ax.set_ylim(-0.8, y - 0.2)
ax.set_xlim(-6, 32)
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
