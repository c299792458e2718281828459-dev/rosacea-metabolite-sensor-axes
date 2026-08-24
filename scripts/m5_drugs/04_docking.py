#!/usr/bin/env python3
"""04_docking.py — AutoDock Vina docking of approved drugs to candidate sensor targets.
Pairs (pharmacologically corrected per reviewers):
  - Ramatroban (CRTH2/DP2 = PTGDR2 antagonist + TXA2R antagonist)  -> PTGDR2 (AlphaFold AF-Q9Y5Y4-F1)
  - TAK-242/resatorvid (TLR4 inhibitor, cyclohexene derivative)    -> TLR4 TIR domain (4G8A)
  - Doxycycline (MMP inhibitor; clinical use in rosacea)           -> MMP9 catalytic domain (1L6J)
Note: Vina is a rigid-receptor docking tool; scores are approximate binding
affinities (kcal/mol), presented as in silico illustration only (no MD here).
"""
import subprocess, sys, os, glob

VINA = "tools/vina"
os.makedirs("output/m5_drugs", exist_ok=True)

# (label, receptor_pdbqt, ligand_pdbqt, center xyz, box size xyz, exhaustiveness)
JOBS = [
    ("ramatroban_PTGDR2", "data/docking/pdbs/ptgdr2_rec.pdbqt",
     "data/docking/ligands/ramatroban.pdbqt",
     [16.0, -8.0, -4.0], [30, 30, 30], 32),          # TM bundle center (AF coords)
    ("TAK242_TLR4", "data/docking/pdbs/4g8a_rec.pdbqt",
     "data/docking/ligands/TAK-242.pdbqt",
     [6.0, -5.0, 15.0], [28, 28, 28], 32),          # TIR dimer interface region
    ("doxycycline_MMP9", "data/docking/pdbs/1l6j_rec.pdbqt",
     "data/docking/ligands/doxycycline.pdbqt",
     [34.2, 52.6, 45.1], [22, 22, 22], 32),         # MMP9 catalytic zinc (1L6J)
]

results = []
for label, rec, lig, center, size, exh in JOBS:
    if not os.path.exists(rec) or not os.path.exists(lig):
        print(f"SKIP {label}: missing input"); continue
    out = f"output/m5_drugs/{label}_out.pdbqt"
    cmd = [VINA, "--receptor", rec, "--ligand", lig,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
           "--exhaustiveness", str(exh), "--num_modes", "9", "--out", out,
           "--seed", "20260824"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    print(f"=== {label} ===")
    for line in p.stdout.splitlines():
        if any(k in line for k in ("mode", "1 ", "2 ", "3 ", "4 ", "5 ", "6 ", "7 ", "8 ", "9 ")):
            print(line)
    if p.returncode != 0:
        print("STDERR:", p.stderr[-300:])
    else:
        results.append((label, out))
    print()

print("done")
