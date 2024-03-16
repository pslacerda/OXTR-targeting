from glob import glob
from os.path import basename, exists
from os import chdir
from subprocess import check_call, CalledProcessError

from yadt.docking import box
from yadt.ftmap import load_ftmap

from pymol import cmd as pm


ADT_PYTHONSH = "/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"
ADT_PREPARE_RECEPTOR = "/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"
ADT_PREPARE_LIGAND = "/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"


SKIP = [
    "2QL7",
    "3FT2",
    "4EXH",
    "1IAU",
    "1KUG",
    "1KUI",
    "1KUK",
    "2J9A",
    "3BOO",
    "3REW",
    "4APH",
    "4IN9",
    "5EFJ",
    "5T5J",
    "5T5L",
    "5T5P",
    "5WEU",
    "6F4R",
    "6F4T",
    "4II9",
    "6RD2",
    "2AOH",
    "1KUI",
    "1KUK",
    "1KUG",
    "4EXH",
    "1IAU",
    "1Y3G",
    "1QBQ",
    "3D3X",
    "3C88",
    "3C89",
    "3C8B",
    "3DDA",
    "3DDB",
    "5T5P",
    "5T5J",
    "5T5L",
    "5AB0",
    "2J9A",
    "5EFJ",
    "3BOO",
    "6F4T",
    "6F4R",
    "4Q4I",
    "4ZHM",
    "2EAX",
    "3AFK0"
]


def run(shell):
    return check_call(shell, shell=True)


chdir("./tests/Atlas")
for idx, atlas_pdb in enumerate(glob("*_atlas.pdb")):
    pdb = basename(atlas_pdb)[:4]
    pept_pdb = f"{pdb}_peptide.pdb"
    prot_pdb = f"{pdb}_protein.pdb"

    print(idx, pdb)
    if pdb in SKIP or not exists(pept_pdb) or not exists(prot_pdb):
        continue

    pm.reinitialize()
    load_ftmap(f"{pdb}_atlas.pdb", origin="atlas")
    try:
        (sx, sy, sz), (cx, cy, cz) = box("*.D.*", 5)
    except:
        continue
    
    try:
        if not exists(pept_pdb + "qt"):
            run(f"{ADT_PYTHONSH} {ADT_PREPARE_LIGAND} -l {pept_pdb} -o {pept_pdb}qt")
        if not exists(prot_pdb + "qt"):
            run(f"{ADT_PYTHONSH} {ADT_PREPARE_RECEPTOR} -r {prot_pdb} -o {prot_pdb}qt")
    except CalledProcessError:
        continue

    out_pdbqt = f"{pdb}_poses.pdbqt"
    out_log = f"{pdb}_log.txt"
    if not exists(out_pdbqt):
        run(
            f"qvina2"
            f" --receptor {prot_pdb}qt"
            f" --ligand {pept_pdb}qt"
            f" --size_x {sx} --size_y {sy} --size_z {sz}"
            f" --center_x {cx} --center_y {cy} --center_z {cz}"
            f" --out {out_pdbqt}"
            f" --log {out_log}"
            f" --seed 1337"
        )
