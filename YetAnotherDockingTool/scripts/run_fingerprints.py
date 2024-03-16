from openbabel import pybel as ob
from glob import glob
from os.path import basename, dirname
from csv import DictWriter
import dask.dataframe as dd

from pymol import cmd as pm

from yadt.ftmap import load_ftmap, dce, dc, fo

def parse_pdbqt_poses(fname):
    models = []
    with open(fname) as file:
        for line in file:
            if line.startswith("MODEL"):
                model = int(line.split()[1])
            elif "VINA RESULT" in line:
                dg = float(line.split()[3])
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
            elif "ROOT" in line:
                models.append({
                    "Model": model,
                    "Dg": dg,
                    "Torsions": torsions,
                })
    return models


fieldnames = set()
def parse_fpt(fname, format="pdbqt"):
    mol = next(ob.readfile(filename=fname, format=format))
    mol.addh()
    bits = {}

    fp = mol.calcfp("FP2")
    for bit in fp.bits:
        bits["FP2_%s" % bit] = 1

    fp = mol.calcfp("MACCS")
    for bit in fp.bits:
        bits["MACCS_%s" % bit] = 1

    fp = mol.calcfp("ECFP10")
    for bit in fp.bits:
        bits["ECFP10_%s" % bit] = 1
    return bits


header = set()
feats = []
for idx, poses_pdbqt in enumerate(glob('tests/Atlas/*_poses.pdbqt')):
    pdb = basename(poses_pdbqt).split('_', 1)[0]
    print(idx, pdb)

    atlas_pdb = f"{dirname(poses_pdbqt)}/{pdb}_atlas.pdb"
    single_pdbqt = f'tests/Atlas/{pdb}_peptide.pdbqt'
    
    mol = next(ob.readfile(filename=single_pdbqt, format="pdbqt"))
    fpt = parse_fpt(single_pdbqt)

    pm.reinitialize()
    pm.load(single_pdbqt, 'PEP_REF')
    pm.load(poses_pdbqt, 'PEP')

    atlas = load_ftmap(atlas_pdb, origin="atlas")
    for model in parse_pdbqt_poses(poses_pdbqt):
        for combination in atlas:
            if combination.kozakov_class in ['D', 'B']:
                feats.append({
                    'Pdb': pdb,
                    'S': combination.strength,
                    'S0': combination.strength0,
                    'MD': combination.max_dist,
                    'CD': combination.center_center,
                    'RMSD': pm.rms_cur('PEP', 'PEP_REF', mobile_state=model['Model']),
                    'FO1': fo(combination.selection, 'PEP', state2=model['Model'], verbose=False),
                    'FO2': fo('PEP', combination.selection, state1=model['Model'], verbose=False),
                    'DC': dc('PEP', combination.selection, state1=model['Model'], verbose=False),
                    'DCE': dce('PEP', combination.selection, state1=model['Model'], verbose=False),
                    **fpt,
                    **model,
                })
                for field in feats[-1]:
                    header.add(field)

print(len(feats))
with open('tests/fingerprints.csv', 'w', newline='') as file:
    writer = DictWriter(file, fieldnames=header)
    writer.writeheader()
    writer.writerows(feats)

