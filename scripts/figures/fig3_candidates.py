#!/usr/bin/env python3
"""fig3_candidates.py — Figure 2: candidate microbial metabolites -> genera -> sensors.
Text-flow layout: genera and sensors are wrapped text blocks (arrows encode
differential direction), so overlaps are structurally impossible.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv, textwrap

rows = list(csv.DictReader(open("output/m2_microbiome/genus_metabolite_sensor_evidence.tsv"), delimiter="\t"))
rows = [r for r in rows if r["metabolite_sensors_MEBOCOST"]]
genus_dir = {r["genus"]: r["genus_direction_Joura2024"] for r in rows}
mets = {}
for r in rows:
    mets.setdefault(r["metabolite"], {"sensors": set(), "genera": set()})
    mets[r["metabolite"]]["sensors"].add(r["metabolite_sensors_MEBOCOST"])
    mets[r["metabolite"]]["genera"].add(r["genus"])

filter_status = {
    "Acetic acid": ("○", "not measurable (volatile SCFA)"),
    "Propionic acid": ("○", "not measurable (volatile SCFA)"),
    "Butyric acid": ("○", "not measurable (volatile SCFA)"),
    "Succinic acid": ("▲", "up in rosacea serum (Li 2025, log2FC = 0.50, P = 1e-4)"),
    "5-Aminolevulinic acid": ("▲", "up in PPR plasma (Zhang 2025)"),
    "Heme (via porphyrin degradation)": ("□", "no concordant circulating signal"),
    "Adenosylcobalamin (vitamin B12)": ("□", "no concordant circulating signal"),
}
arrow_sym = {"up": "↑", "down": "↓", "n.s.": "~"}

order = sorted(mets, key=lambda m: -len(mets[m]["genera"]))

# compute layout first
line_h = 0.235
rows_meta = []
for m in order:
    gen = sorted(mets[m]["genera"], key=lambda g: -genus_dir[g].count("d"))
    gen_str = ", ".join(f"{g}{arrow_sym[genus_dir[g]]}" for g in gen)
    sens_str = "; ".join(sorted(mets[m]["sensors"]))
    gen_lines = textwrap.wrap(gen_str, width=50)
    sen_lines = textwrap.wrap(sens_str, width=32)
    n = max(len(gen_lines), len(sen_lines), 1)
    rows_meta.append((m, gen_lines, sen_lines, n))
total = sum(n*line_h + 0.26 for m,_,_,n in rows_meta) + 0.75  # header + legend
fig, ax = plt.subplots(figsize=(7.0, total), dpi=600)
ax.set_xlim(0, 7.0); ax.set_ylim(0, total); ax.axis("off")

# columns (inches): name 0-1.35 | genera 1.5-4.55 | sensors 4.7-7.2
y = total - 0.55
for m, gen_lines, sen_lines, n in rows_meta:
    sym = filter_status.get(m, ("", ""))[0]
    # metabolite name
    ax.text(1.3, y - (n-1)*line_h/2, f"{sym} {m}", ha="right", va="center", fontsize=8)
    # genera
    for i, line in enumerate(gen_lines):
        ax.text(1.45, y - i*line_h, line, ha="left", va="center", fontsize=7.5, color="#333333")
    # sensors
    for i, line in enumerate(sen_lines):
        ax.text(4.6, y - i*line_h, line, ha="left", va="center", fontsize=7.5, color="#4a148c")
    y -= n*line_h + 0.26

# column headers
ax.text(1.3, total-0.18, "Candidate metabolite", ha="right", fontsize=8, style="italic", color="#666666")
ax.text(1.45, total-0.18, "Producing skin genera", ha="left", fontsize=8, style="italic", color="#666666")
ax.text(4.6, total-0.18, "Annotated human sensors (MEBOCOST)", ha="left", fontsize=8, style="italic", color="#666666")
ax.axhline(y + 0.16, color="#999999", lw=0.6)
ax.axvline(1.37, color="#999999", lw=0.6)
ax.axvline(4.52, color="#999999", lw=0.6)

# legend at bottom
leg_y = y - 0.1
ax.text(0.05, leg_y, "Genus direction (Joura 2024): ↑ enriched in rosacea; ↓ depleted; ~ dominant, not significant\n"
        "Metabolomics filter: ▲ concordant circulating signal; ○ not measurable by LC-MS (volatile SCFAs);\n"
        "□ measurable but no concordant circulating signal",
        fontsize=6, va="top", ha="left", color="#444444")
plt.savefig("output/figures/Fig3_candidate_metabolites.png", dpi=600, bbox_inches="tight")
print("Fig2 saved (text-flow layout)")
