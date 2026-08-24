#!/usr/bin/env python3
"""02_genus_metabolite_evidence.py
Build the genus -> metabolite -> sensor evidence table for M2.
- Genus differential status: Joura et al. 2024 Table 2 (skin, published; PMID 39770869)
- Genus -> metabolite: curated literature (PMID), flags for verification
- Metabolite -> sensor: MEBOCOST human DB (human_met_sensor_update_Oct21_2025.tsv),
  evidence column = PMIDs from the database entry
Output: output/m2_microbiome/genus_metabolite_sensor_evidence.tsv
"""
import csv

# Joura 2024 Table 2 skin genus status (published values)
joura = {
    "Cutibacterium": ("n.s.", 0.088, "dominant (35.1% rosacea / 20.9% control)"),
    "Staphylococcus": ("n.s.", 0.152, "dominant (14.0% / 8.8%)"),
    "Neisseria": ("up", 0.023, ""),
    "Corynebacterium": ("up", 0.015, ""),
    "Bacteroides": ("down", 0.001, ""),
    "Faecalibacterium": ("down", 0.0009, ""),
    "Prevotella": ("down", 0.001, ""),
    "Blautia": ("down", 0.001, ""),
    "Ruminococcus": ("down", 0.001, ""),
    "Subdoligranulum": ("down", 0.017, ""),
}

# Genus -> metabolites (curated; PMIDs to verify flagged with *)
genus_met = [
    # (genus, metabolite, direction_note, evidence_pmid, verify)
    ("Cutibacterium", "Acetic acid", "fermentation product of C. acnes", "23405142", False),
    ("Cutibacterium", "Propionic acid", "fermentation product of C. acnes", "23405142", False),
    ("Cutibacterium", "Butyric acid", "fermentation product of C. acnes", "23405142", False),
    ("Cutibacterium", "Porphyrins (uro/coproporphyrin III)", "strain-dependent porphyrin production", "27303708", False),
    ("Cutibacterium", "5-Aminolevulinic acid", "porphyrin precursor in C. acnes", "27303708", False),
    ("Cutibacterium", "Adenosylcobalamin (vitamin B12)", "C. acnes B12 biosynthesis (Joura LefSe: higher in healthy skin)", "26109103", False),
    ("Cutibacterium", "Heme (via porphyrin degradation)", "Joura 2024 LefSe: excessive heme production in rosacea skin", "39770869", False),
    ("Staphylococcus", "Acetic acid", "S. epidermidis glycerol fermentation", "23405142", False),
    ("Staphylococcus", "Butyric acid", "S. epidermidis butyrate production", "31159213", False),
    ("Corynebacterium", "Acetic acid", "axillary corynebacteria VFA (acetic acid) production", "18494871", False),
    ("Neisseria", "LPS/lipid A", "Gram-negative outer membrane (TLR4 ligand)", "9851930", False),
    ("Faecalibacterium", "Butyric acid", "F. prausnitzii major butyrate producer", "19222573", False),
    ("Bacteroides", "Propionic acid", "Bacteroidetes propionate (succinate pathway)", "24553467", False),
    ("Bacteroides", "Succinic acid", "Bacteroides succinate secretion (SUCNR1 ligand; Tannahill PMID 23535595)", "24553467", False),
    ("Bacteroides", "Acetic acid", "Bacteroides acetate production", "19222573", False),
    ("Prevotella", "Succinic acid", "Prevotella succinate production (SUCNR1 ligand; Tannahill PMID 23535595)", "24553467", False),
    ("Prevotella", "Acetic acid", "Prevotella acetate production", "19222573", False),
    ("Blautia", "Acetic acid", "Blautia spp. acetate production", "19222573", False),
    ("Ruminococcus", "Acetic acid", "R. bromii acetate/formate production", "22343308", False),
    ("Subdoligranulum", "Butyric acid", "S. variabile butyrate production", "19222573", False),
]

# Load MEBOCOST sensor DB
me_rows = list(csv.DictReader(open(
    "data/MEBOCOST-main/data/mebocost_db/human/human_met_sensor_update_Oct21_2025.tsv"), delimiter="\t"))
# name lookup with synonyms
me_by_name = {}
for r in me_rows:
    names = set(x.strip().lower() for x in r["metName"].split(";") if x.strip())
    names.add(r["standard_metName"].strip().lower())
    me_by_name[r["standard_metName"].strip().lower()] = r
    for n in names:
        me_by_name.setdefault(n, r)

def find_sensors(met_name):
    key = met_name.lower().split("(")[0].strip()
    out = []
    if key in me_by_name:
        r = me_by_name[key]
        for rr in me_rows:
            syns = set(x.strip().lower() for x in rr["metName"].split(";") if x.strip())
            syns.add(rr["standard_metName"].lower())
            if key in syns:
                out.append((rr["Gene_name"], rr["Annotation"], rr["Evidence"]))
    return out

with open("output/m2_microbiome/genus_metabolite_sensor_evidence.tsv", "w") as f:
    f.write("\t".join(["genus", "genus_direction_Joura2024", "genus_p_value", "note",
                       "metabolite", "metabolite_sensors_MEBOCOST", "sensor_type",
                       "genus_met_evidence_pmid", "sensor_evidence_pmid", "verify_flag"]) + "\n")
    for genus, met, note, pmid, verify in genus_met:
        sens = find_sensors(met)
        if not sens:
            f.write("\t".join([genus, joura.get(genus, ("", "", ""))[0],
                               str(joura.get(genus, ("", "", ""))[1]), joura.get(genus, ("", "", ""))[2],
                               met, "", "", pmid, "", str(verify)]) + "\n")
            continue
        for g, at, ev in sens:
            f.write("\t".join([genus, joura.get(genus, ("", "", ""))[0],
                               str(joura.get(genus, ("", "", ""))[1]), joura.get(genus, ("", "", ""))[2],
                               met, g, at, pmid, ev, str(verify)]) + "\n")
print("written")
# summary
rows = list(csv.DictReader(open("output/m2_microbiome/genus_metabolite_sensor_evidence.tsv"), delimiter="\t"))
print("total evidence rows:", len(rows))
print("with sensor:", sum(1 for r in rows if r["metabolite_sensors_MEBOCOST"]))
from collections import Counter
print(Counter((r["genus"], r["metabolite"]) for r in rows if r["metabolite_sensors_MEBOCOST"]))
