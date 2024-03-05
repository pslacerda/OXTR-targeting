import sys
import os.path
from vina import Vina
from importlib.machinery import SourceFileLoader
from glob import glob
from subprocess import run
from shutil import copy


config_fname = sys.argv[-1]
config = SourceFileLoader("config", config_fname).load_module()

vina = Vina(sf_name=config.FUNCTION, verbosity=config.VERBOSITY, cpu=config.CPU)
vina.set_receptor(os.path.expanduser(config.RECEPTOR))
vina.compute_vina_maps(config.CENTER, config.BOX_SIZE)

for idx, fname in enumerate(glob('NuBBE_Files_MOL2/*.mol2')):
    try:
        name = os.path.basename(fname).split('.')[0]
        print(idx, name)
        out_fname = f"{config.OUTPUT_DIRECTORY}/{name}.pdbqt"
        if os.path.exists(out_fname):
            continue
        # generate a PDBQT conformatio
        
        wd = os.getcwd()
        try:
            any = False
            with open(fname) as file:
                for line in file:
                    if line.strip().split()[-1] == 'AC':
                        any = True
            if any:
                continue
            copy(fname, os.path.dirname(out_fname) + f"/{name}.mol2")
            os.chdir(os.path.dirname(out_fname))
            pythonsh = "/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"
            prepare = '/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
            run(f'{pythonsh} {prepare} -l {fname} -o {out_fname}', shell=True)

            
            vina.set_ligand_from_file(out_fname)
            vina.dock(config.EXHAUSTIVENESS)
            vina.write_poses(out_fname, 5, overwrite=True)
        except:
            print("Failed to dock", name)
            os.unlink(out_fname)
        finally:
            os.chdir(wd)
        if not os.path.exists(out_fname):
            print(f"Failed to dock ")
    except Exception as exc:
        print(exc)
        continue
