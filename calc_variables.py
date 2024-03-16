from openbabel import pybel as ob
from glob import glob
from csv import DictWriter
from os.path import basename
from pymol import cmd as pm
from yadt.ftmap import load_ftmap, dce, dc, fo, load_ftmap
from PyFingerprint.fingerprint import get_fingerprint
import click


def parse_pdbqt_poses(fname):
    with open(fname) as file:
        for line in file:
            if line.startswith("MODEL"):
                model = int(line.split()[1])
            elif "VINA RESULT" in line:
                dg = -float(line.split()[3])
                return {
                    "Model": model,
                    "Dg": dg,
                }


fieldnames = set(
    [
        "Molecule.ChEMBL.ID",
        "Dg",
        # "S",
        # "S0",
        # "MD",
        # "CD",
        # "FO1",
        # "FO2",
        # "DC",
        # "DCE",
        "Model",
        "NumAtoms",
    ]
)


def parse_fpt(smi):
    def _inner(method):
        res = {}

        fpts = get_fingerprint(smi, method)
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


parse_fpt("C")


feats = []
cnt = 0


@cli.command()
@click.option("-p", "--hotspot-program", type=click.Choice(["ftmap", "atlas"]))
@click.argument("hs_file")
@click.argument("output_file")
@click.argument("docking_wildcards", n=-10)
def run_fpts(hs_file, output_file, hotspot_program, docking_wildcards):
    hs = load_ftmap(hs_file, origin=hotspot_program)

    with open(output_file, "w", newline="") as file:
        writer = DictWriter(file, fieldnames=sorted(fieldnames), extrasaction="ignore")
        writer.writeheader()

        for idx, poses_pdbqt in enumerate(flatten_files(docking_wildcards)):
            pdb = basename(poses_pdbqt).split(".", 1)[0]

            print(cnt, pdb)
            cnt += 1
            mol = next(ob.readfile(filename=poses_pdbqt, format="pdbqt"))
            fpt = parse_fpt(mol.write())

            pm.delete("PEP")
            pm.load(poses_pdbqt, "PEP")

            sel = " or ".join(hs.selection for hs in atlas)
            S = sum([hs.strength for hs in atlas])
            import numpy as np

            CD = np.average([hs.center_center for hs in atlas])
            MD = np.average([hs.max_dist for hs in atlas])
            model = parse_pdbqt_poses(poses_pdbqt)
            if model is not None:
                writer.writerow(
                    {
                        "Molecule.ChEMBL.ID": pdb,
                        "S": S,
                        "MD": MD,
                        "CD": CD,
                        "FO1": fo(sel, "PEP", state2=1, verbose=False),
                        "FO2": fo("PEP", sel, verbose=False),
                        "DC": dc("PEP", sel, state1=1, verbose=False),
                        "DCE": dce("PEP", sel, state1=1, verbose=False),
                        "NumAtoms": pm.count_atoms("PEP"),
                        **fpt,
                        **model,
                    }
                )
