from subprocess import CalledProcessError, check_call
from glob import iglob, glob
from os.path import basename, exists
from typing import Iterable
from shutil import which

import numpy as np
from pymol import cmd as pm, plugins, CmdException
from openbabel import pybel as ob
from pymol import cgo


__all__ = ["box", "dock"]


@pm.extend
def box(sel, margin, display=0):
    """
    Create a box around a selection.
    USAGE:
        box sel, margin, display=0
    EXAMPLE:
        box resn CU and chain A, 5, display=1
    """
    margin = int(margin)
    box_coords = pm.get_coords(sel)

    if box_coords is None:
        raise Exception(f'Selection "{sel}" doesn\'t evaluate to any atom set.')

    max = np.max(box_coords, axis=0) + margin
    min = np.min(box_coords, axis=0) - margin

    half_size = (max - min) / 2
    center = min + half_size

    size_x, size_y, size_z = half_size * 2
    center_x, center_y, center_z = center

    size_x, size_y, size_z = (
        round(float(size_x), 2),
        round(float(size_y), 2),
        round(float(size_z), 2),
    )
    center_x, center_y, center_z = (
        round(float(center_x), 2),
        round(float(center_y), 2),
        round(float(center_z), 2),
    )
    box = (size_x, size_y, size_z), (center_x, center_y, center_z)
    if bool(int(display)):
        display_box(
            (
                (center_x - size_x / 2, center_x + size_x / 2),
                (center_y - size_y / 2, center_y + size_y / 2),
                (center_z - size_z / 2, center_z + size_z / 2),
            )
        )
    return box


@pm.extend
def dock(
    command,
    program="vina",
    ligands_file="ligands.smi",
    box_sel="not polymer.protein",
    box_margin=4,
    wdir=".",
    target="polymer.protein",
    max_ligands=10,
    group=None,
):
    if command == "screening":
        eng = DockingEngine(program=program, wdir=wdir, group=group)
        eng.set_box(box_sel, box_margin)
        eng.set_target(target)
        eng.dock_from_smiles(ligands_file)
    elif command == "analyze":
        ana = DockingAnalyzer(program=program, wdir=wdir)
        ana.set_box(box_sel, box_margin)
        ana.set_target(target)
        ana.analyze(group, max_ligands)
    else:
        print(f'Unknown command "{command}", must be one of "screening" or "analyze".')


def display_box(box, name="_box"):
    """https://pymolwiki.org/index.php/Autodock_plugin"""

    view = pm.get_view()
    pm.delete(name)
    obj = []
    # build cgo object
    color = [1.0, 1.0, 1.0]
    cylinder_size = 0.2
    for i in range(2):
        for k in range(2):
            for j in range(2):
                if i != 1:
                    obj.append(cgo.CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i + 1], box[1][j], box[2][k]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(cgo.COLOR)
                    obj.extend(color)
                    obj.append(cgo.SPHERE)
                    obj.extend([box[0][i], box[1][j], box[2][k], cylinder_size])

                if j != 1:
                    obj.append(cgo.CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i], box[1][j + 1], box[2][k]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(cgo.COLOR)
                    obj.extend(color)
                    obj.append(cgo.SPHERE)
                    obj.extend([box[0][i], box[1][j + 1], box[2][k], cylinder_size])
                if k != 1:
                    obj.append(cgo.CYLINDER)
                    obj.extend([box[0][i], box[1][j], box[2][k]])
                    obj.extend([box[0][i], box[1][j], box[2][k + 1]])
                    obj.append(cylinder_size)
                    obj.extend(color)
                    obj.extend(color)
                    obj.append(cgo.COLOR)
                    obj.extend(color)
                    obj.append(cgo.SPHERE)
                    obj.extend([box[0][i], box[1][j], box[2][k + 1], cylinder_size])
    pm.load_cgo(obj, name)
    pm.set_view(view)


class DockingError(Exception):
    pass


class InternalDockingError(Exception):
    pass


class DockingEngine:
    def __init__(self, wdir, **options):
        self.wdir = wdir
        if options.get("program") not in ["vina", "qvina2", "qvinaw"]:
            raise DockingError(
                'Option "program" must be "vina" or "qvina2" or "qvinaw".'
            )
        self.program = which(options.pop("program"))
        if self.program is None:
            raise DockingError(f'Program "{self.program}" not available.')
        if "log" in options or "out" in options:
            raise DockingError('Options "--log" or "--out" cannot be set.')
        self.options = options
    
    def run(self, cmd):
        try:
            check_call(cmd)
        except CalledProcessError as exc:
            print(f'Command "{cmd}" failed.')
            raise InternalDockingError from exc
    
    def _get_options(self, name):
        return " ".join(
            "--%s %s" % (k, v)
            for (k, v) in dict(
                log=self.get_log_filename(name),
                out=self.get_poses_filename(name),
                ligand=self.get_prepared_filename(name),
                **self.options,
            ).items()
        )

    def set_box(self, sel, margin=4):
        (
            (self.options["size_x"], self.options["size_y"], self.options["size_z"]),
            (
                self.options["center_x"],
                self.options["center_y"],
                self.options["center_z"],
            ),
        ) = box(sel, margin)

    def set_target(self, sel):
        pythonsh = plugins.preferences["ADT_PYTHONSH"]
        prepare_receptor4 = plugins.preferences["ADT_PREPARE_RECEPTOR4"]

        filename_pdb = f"{self.wdir}/target.pdb"
        filename_pdbqt = f"{filename_pdb}qt"

        try:
            pm.save(filename_pdb, sel)
        except CmdException as exc:
            raise DockingError(f"Failed to save the target selection '{sel}'.") from exc
        self.run(
            f"{pythonsh} {prepare_receptor4} -r {filename_pdb} -o {filename_pdbqt}"
        )
        assert exists(filename_pdbqt)
        self.options["receptor"] = filename_pdbqt

    def get_poses_filename(self, name):
        return f"{self.wdir}/{name}.poses.pdbqt"

    def get_prepared_filename(self, name):
        return f"{self.wdir}/{name}.pdbqt"

    def get_log_filename(self, name):
        return f"{self.wdir}/{name}.log.txt"

    def dock_from_smiles(self, ligands_files):
        if "receptor" not in self.options or "center_x" not in self.options:
            raise DockingError("Must set the receptor and box before running docking.")
        if isinstance(ligands_files, str):
            ligands_files = [ligands_files]
        assert isinstance(ligands_files, Iterable)

        for ligand_file in ligands_files:
            with open(ligand_file) as ligand_file:
                line = ligand_file.readline().strip()
                while line:
                    smiles, name = line.split()
                    prepared_fname = self.get_prepared_filename(name)
                    poses_fname = self.get_poses_filename(name)

                    if exists(poses_fname):
                        continue

                    if not exists(prepared_fname):
                        step1_fname = f"{self.wdir}/{name}.pdb"
                        mol = ob.readstring("smi", smiles)
                        mol.addh()
                        mol.make3D()
                        mol.write("pdb", filename=step1_fname, overwrite=True)
                        pythonsh = plugins.preferences["ADT_PYTHONSH"]
                        prepare_ligand4 = plugins.preferences["ADT_PREPARE_LIGAND4"]
                        try:
                            self.run(
                                f"{pythonsh} {prepare_ligand4} -l {step1_fname} -o {prepared_fname}"
                            )
                        except InternalDockingError as exc:
                            continue
                        assert exists(prepared_fname)

                    program = self.program
                    opts = self._get_options(name)
                    try:
                        self.run(f"{program} {opts}")
                    except InternalDockingError as e:
                        print(f'Failed to dock "{name}".')

                    line = ligand_file.readline().strip()


class DockingAnalyzer:
    def __init__(self, wdir, group=None, **options):
        self.wdir = wdir
        if options.get("program") not in ["vina", "qvina2", "qvinaw"]:
            raise DockingError(
                'Option "program" must be "vina" or "qvina2" or "qvinaw".'
            )
        self.program = which(options.pop("program"))
        if group is None:
            self.group = basename(self.program)
        else:
            self.group = group

        if self.program is None:
            raise DockingError(f'Program "{self.program}" not available.')
        if "log" in options or "out" in options:
            raise DockingError('Options "--log" or "--out" cannot be set.')
        self.options = options

    def set_box(self, sel, margin=4):
        (
            (self.options["size_x"], self.options["size_y"], self.options["size_z"]),
            (
                self.options["center_x"],
                self.options["center_y"],
                self.options["center_z"],
            ),
        ) = box(sel, margin, display=True)

    def set_target(self, sel=None):
        filename_pdbqt = f"{self.wdir}/target.pdbqt"
        if not exists(filename_pdbqt):
            raise DockingError(f'Target file "{filename_pdbqt}" doesn\'t exists.')
        self.options["receptor"] = filename_pdbqt

    def get_poses_filename(self, name):
        return f"{self.wdir}/{name}.poses.pdbqt"

    def get_prepared_filename(self, name):
        return f"{self.wdir}/{name}.pdbqt"

    def get_log_filename(self, name):
        return f"{self.wdir}/{name}.log.txt"

    def analyze(self, group, max_ligands):
        pm.load(self.options["receptor"])

        ligands = []
        for poses_pdbqt in glob(f"{self.wdir}/*.poses.pdbqt"):
            ligand_data = parse_pdbqt_poses(poses_pdbqt)
            ligands.append(ligand_data)
        ligands.sort(key=lambda l: l["Dg"])
        
        max_ligands = int(max_ligands)
        for idx, ligand in enumerate(ligands):
            if idx >= max_ligands:
                break
            dg = ligand["Dg"]
            model_name = basename(ligand['Filename']).rsplit(".", 2)[0]
            model_name = f"{group}.rank.{model_name}_{-dg}"
            pm.group(f"{group}.rank")
            pm.load(poses_pdbqt, model_name)
        pm.group(group, f"{group}.rank")
