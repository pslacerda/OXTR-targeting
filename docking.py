"""
Run AutoDock Vina

PEDRO SOUSA LACERDA
pslacerda@gmail.com
"""

import sys
import os.path
from vina import Vina
from openbabel import pybel as ob
from importlib.machinery import SourceFileLoader
import csv

config_fname = sys.argv[-2]
config = SourceFileLoader("config", config_fname).load_module()

smiles_fname = sys.argv[-1]

vina = Vina(sf_name=config.FUNCTION, verbosity=config.VERBOSITY, cpu=config.CPU)
vina.set_receptor(os.path.expanduser(config.RECEPTOR))
vina.compute_vina_maps(config.CENTER, config.BOX_SIZE)

reader = csv.DictReader(open(smiles_fname), quotechar='"', delimiter='\t') 
for line in reader:
    try:
        # extract molecule name and SMILES
        smiles = line['SMILES']
        name = line['name']
        if "." in smiles:
            print(
                "Compound "
                + name
                + " has disconnected atoms. ("
                + smiles
                + "). Skipping."
            )
            continue
        else:
            print("Processing compound: " + name + "...")

        out_fname = f"{config.OUTPUT_DIRECTORY}/{name}.pdbqt"

        if os.path.exists(out_fname):
            continue
        # generate a PDBQT conformatio
        try:
            print(name, smiles)
            mol = ob.readstring("smi", smiles)
            if mol.molwt > 1500:
                continue
            mol.addh()
            mol.make3D()
            pdbqt_str = mol.write("pdbqt")

            vina.set_ligand_from_string(pdbqt_str)
            vina.dock(config.EXHAUSTIVENESS)
            vina.write_poses(out_fname, 5, overwrite=True)
        except:
            print("Failed to dock", name)
        if not os.path.exists(out_fname):
            print(f"Failed to dock ")
    except Exception as exc:
        print(exc)
        continue
