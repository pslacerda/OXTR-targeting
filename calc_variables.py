from openbabel import pybel as ob
from glob import glob
from csv import DictWriter
from os.path import basename
from pymol import cmd as pm

def parse_pdbqt_poses(fname):
    
    with open(fname) as file:
        if 'VINA' not in '\n'.join(file.readlines()):
            return
    with open(fname) as file:
        for line in file:
            if line.startswith("MODEL"):
                model = int(line.split()[1])
            elif "VINA RESULT" in line:
                dg = -float(line.split()[3])
            elif "INTER + INTRA" in line:
                inter_intra = float(line.split()[4])
            elif "INTER:" in line:
                inter = float(line.split()[2])
            elif " INTRA:" in line:
                intra = float(line.split()[2])
            elif "UNBOUND" in line:
                unbound = float(line.split()[2])
            elif "active torsions" in line:
                torsions = int(line.split()[1])
            elif line.startswith("ROOT"):
                return {
                    "Model": model,
                    "Dg": dg,
                    "Torsions": torsions,
                    "InterIntra": inter_intra,
                    "Inter": inter,
                    "Intra": intra,
                    "Unbound": unbound,
                }

from yadt.ftmap import load_ftmap, dce, dc, fo
from PyFingerprint.fingerprint import get_fingerprint


VERSION = "fda"


fieldnames = set(
    [
        "Molecule.ChEMBL.ID",
        "Dg",
        "Torsions",
        "S",
        "S0",
        "MD",
        "CD",
        "FO1",
        "FO2",
        "DC",
        "DCE",
        "Model",
        "InterIntra",
        "Inter",
        "Intra",
        "Unbound",
        "NumAtoms"
    ]
)


def parse_fpt(
    smi,
):
    def _inner(method):
        res = {}

        fpts = get_fingerprint  (smi, method)
        for idx, bit in enumerate(fpts.to_numpy().T):
            res[f"{method}_{idx}"] = bit
        return res

    res = {}
    for method in [

        "standard",
        "extended",
        "graph",
        "maccs",
        "pubchem",
        "estate",
        "hybridization",
        "lingo",
        "klekota-roth",
        "shortestpath",
        "cdk-substructure",
        "circular",
        "cdk-atompairs",
        "fp2",
        "spectrophore",
    ]:
        fpt = _inner(method)
        res.update(fpt)
        fieldnames.update(fpt.keys())
    return res


atlas_pdb = f"ftmap/6TPK.pdb"
atlas = load_ftmap(atlas_pdb, origin="ftmap")


feats = []
cnt = 0

parse_fpt("C")

with open(f"fingerprints_{VERSION}.csv", "w", newline="") as file:
    writer = DictWriter(file, fieldnames=sorted(fieldnames), extrasaction="ignore")
    writer.writeheader()
    for idx, poses_pdbqt in enumerate(glob(f"results_{VERSION}/*.pdbqt")):
        pdb = basename(poses_pdbqt).split(".", 1)[0]

        print(cnt, pdb)
        cnt += 1
        mol = next(ob.readfile(filename=poses_pdbqt, format="pdbqt"))
        fpt = parse_fpt(mol.write())

        pm.delete("PEP")
        pm.load(poses_pdbqt, "PEP")

        sel = ' or '.join(hs.selection for hs in atlas)
        S = sum([hs.strength for hs in atlas])
        import numpy as np
        CD = np.average([hs.center_center for hs in atlas])
        MD = np.average([hs.max_dist for hs in atlas])
        model = parse_pdbqt_poses(poses_pdbqt)
        if model is not None:
            writer.writerow({
                "Molecule.ChEMBL.ID": pdb,
                "S": S,
                "MD": MD,
                "CD":  CD,
                'FO1': fo(sel, 'PEP', state2=1, verbose=False),
                'FO2': fo('PEP', sel, verbose=False),
                'DC': dc('PEP',sel, state1=1, verbose=False),
                'DCE': dce('PEP', sel, state1=1, verbose=False),
                'NumAtoms': pm.count_atoms('PEP'),
                **fpt,
                **model
            })
