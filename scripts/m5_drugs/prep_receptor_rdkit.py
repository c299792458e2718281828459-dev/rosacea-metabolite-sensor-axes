#!/usr/bin/env python3
"""prep_receptor_rdkit.py — minimal PDB -> PDBQT receptor preparation with RDKit.
Adds polar hydrogens, computes Gasteiger charges, assigns AutoDock4 atom types.
Usage: prep_receptor_rdkit.py input.pdb output.pdbqt [chain(s)]
"""
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def ad_type(atom):
    """Rough AutoDock4 type mapping (sufficient for Vina)."""
    a = atom.GetAtomicNum()
    if a == 6:
        return "A" if atom.GetIsAromatic() else "C"
    if a == 7:
        if atom.GetTotalNumHs() > 0:
            return "N"
        return "NA"
    if a == 8:
        return "OA"
    if a == 16:
        return "S"
    if a == 1:
        # H attached to N/O = polar
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() in (7, 8):
                return "HD"
        return "H"
    if a == 15:
        return "P"
    return "C"  # fallback

def prep(pdb_in, pdbqt_out):
    mol = Chem.MolFromPDBFile(pdb_in, removeHs=False, sanitize=False)
    if mol is None:
        raise SystemExit("could not parse PDB")
    mol = Chem.RemoveHs(mol)
    # keep only protein residues (drop waters/hetero)
    keep = []
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        keep.append(atom.GetIdx())
    mol = Chem.RWMol(mol)
    to_remove = [i for i in range(mol.GetNumAtoms()) if i not in keep]
    for i in sorted(to_remove, reverse=True):
        mol.RemoveAtom(i)
    mol = mol.GetMol()
    mol = Chem.AddHs(mol, addCoords=True)
    AllChem.ComputeGasteigerCharges(mol)
    conf = mol.GetConformer()
    with open(pdbqt_out, "w") as f:
        f.write("REMARK prepared by prep_receptor_rdkit.py\n")
        for atom in mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info is None:
                continue
            pos = conf.GetAtomPosition(atom.GetIdx())
            name = atom.GetPDBProp("_Name") if atom.HasProp("_Name") else atom.GetSymbol()
            name = name.ljust(4)[:4]
            resname = (info.GetResidueName() or "UNK").ljust(3)[:3]
            chain = info.GetChainId().strip() or " "
            resnum = info.GetResidueNumber()
            elem = atom.GetSymbol().rjust(2)[:2]
            charge = float(atom.GetProp("_GasteigerCharge"))
            atype = ad_type(atom)
            serial = (atom.GetIdx() % 100000) + 1
            # PDBQT: cols 7-11 serial, 13-16 name, 18-20 resname, 22 chain,
            # 23-26 resseq, 31-54 xyz (8.3), 55-60 occ, 61-66 charge, 78-79 AD type
            f.write("ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" %
                    (serial, name[:4].rjust(4), resname, chain, resnum % 10000,
                     pos.x, pos.y, pos.z, 1.00, charge, atype))
    print("wrote", pdbqt_out)

if __name__ == "__main__":
    prep(sys.argv[1], sys.argv[2])
