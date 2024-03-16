"""
Yet Anther Docking Tool

PEDRO SOUSA LACERDA
pslacerda@gmail.com

Usge:
    dock config.toml
    dock --version

"""

import os.path
from openbabel import pybel as ob
from glob import glob
from subprocess import check_call, CalledProcessError
import click
import os
from contextlib import contextmanager
from openbabel import pybel as ob
from glob import glob
from csv import DictWriter
from os.path import basename
from pymol import cmd as pm
from yadt.ftmap import load_ftmap, dce, dc, fo
import click


@contextmanager
def cd(path):
    wd = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(wd)


def _mgltools(path):
    pythonsh = f"{path}/bin/pythonsh"
    prepare_receptor4 = (
        f"{path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
    )
    prepare_ligand4 = (
        f"{path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    )
    if not (
        os.path.exists(pythonsh)
        and os.path.exists(prepare_receptor4)
        and os.path.exists(prepare_ligand4)
    ):
        raise Exception("Misconfigured MGLTOOLS")
    return pythonsh, prepare_receptor4, prepare_ligand4


def run(args):
    try:
        check_call(args, shell=True)
        return True
    except CalledProcessError:
        return False


@click.group()
def cli():
    pass


@cli.command()
@click.option("-c", "--config", default="docking.cfg")
@click.option("-e", "--engine", default="vina")
@click.argument("output-directory")
@click.argument("prepared-ligands", nargs=-1)
def dock(engine, config, prepared_ligands, output_directory):
    for ligand in prepared_ligands:
        try:
            ofname = os.path.basename(ligand)
            ofname = os.path.splitext(ofname)[0]
            print(ofname)
            ofname = f"{output_directory}/{ofname}.pdbqt"
            if os.path.exists(ofname):
                continue
            run(f"{engine} --config {config} --ligand {ligand} --out {ofname}")
        except Exception as exc:
            print(exc)


@cli.command()
@click.option("--mgltools", envvar="MGLTOOLS")
@click.option("-m", "--method", type=click.Choice(["adt4", "obabel"]))
@click.option("--overwrite/--no-overwrite", default=False)
@click.option("--names-from-filename", is_flag=True, default=False)
@click.option('--mwt', nargs=2, type=float, default=None)
@click.argument("output-directory")
@click.argument("input-wildcard", nargs=-1)
def prepare_compounds(
    mgltools, method, overwrite, output_directory, input_wildcard, names_from_filename, mwt
):
    filenames = flatten_files(input_wildcard)
    for fname in filenames:
        base, ext = os.path.splitext(fname)
        base = os.path.basename(base)
        if ext not in [".smi", ".pdb", ".mol2", ".pdbqt"]:
            raise Exception(f"Unsupported filetype {ext}")
        mols = ob.readfile(ext[1:], fname)
        for idx, mol in enumerate(mols):
            if mwt and not (mwt[0] < mol.molwt < mwt[1]):
                continue
            if ext == ".smi":
                mol.make3D()
            if names_from_filename:
                dst_fname = f"{output_directory}/{base}.pdbqt"
            else:
                dst_fname = f"{output_directory}/{mol.title}.pdbqt"
            if os.path.exists(dst_fname) and overwrite:
                continue
            mol.write("pdbqt", filename=dst_fname, overwrite=overwrite)
            if method == 'obabel':
                mol.write("pdbqt", filename=dst_fname, overwrite=overwrite)    
            elif method == 'adt4':
                dst_fname = os.path.abspath(dst_fname)
                if os.path.exists(dst_fname) and not overwrite:
                    continue
                pythonsh, prepare_receptor4, prepare_ligand4 = _mgltools(mgltools)
                fname = os.path.abspath(fname)
                with cd(os.path.dirname(fname)):
                    fname = os.path.basename(fname)
                    success = run(f"{pythonsh} {prepare_ligand4} -l {fname} -o {dst_fname}")
                    if not success:
                        click.echo(f"Failed on {fname}")
            else:
                raise Exception("unexpected error")
            break
            

@cli.command()
@click.option("--mgltools", envvar="MGLTOOLS")
@click.option("--overwrite/--no-overwrite", default=False)
@click.argument("input")
@click.argument("output")
def prepare_receptor(mgltools, overwrite, input, output):
    pythonsh, prepare_receptor4, prepare_ligand4 = _mgltools(mgltools)
    if not overwrite and os.path.exists(output):
        return
    run(f"{pythonsh} {prepare_receptor4} -r {input} -o {output}")


def flatten_files(patterns):
    patterns = list(patterns)
    filenames = []
    while patterns:
        pat = patterns.pop()
        for filename in glob(pat):
            filenames.append(filename)
    return filenames


fieldnames = set(
    [
        "Name",
        "Dg",
        "S",
        "S0",
        "MD",
        "CD",
        "FO1",
        "FO2",
        "DC",
        "DCE",
        "Model",
        "NumAtoms",
        "Molecular.Weight"
    ]
)


def parse_fpt(smi):
    from PyFingerprint.fingerprint import get_fingerprint
    
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
            
@cli.command()
@click.option('-p', '--hotspot-program', default='ftmap', type=click.Choice(['ftmap', 'atlas']))
@click.argument('hs_file')
@click.argument('output_file')
@click.argument('docking_wildcards', nargs=-1)
def calc_fpts(hs_file, output_file, hotspot_program, docking_wildcards):
    parse_fpt("C")

    cnt = 0
    atlas = load_ftmap(hs_file, origin=hotspot_program)

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
            S0 = next(hs.strength for hs in atlas)
            
            import numpy as np

            CD = np.average([hs.center_center for hs in atlas])
            MD = np.average([hs.max_dist for hs in atlas])

            model = parse_pdbqt_poses(poses_pdbqt)
            if model is not None:
                row = {
                        "Name": pdb,
                        "S": S,
                        "S0": S0,
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
                writer.writerow(row)
                file.flush()

def train_regressor():
    pass

if __name__ == "__main__":
    cli()
