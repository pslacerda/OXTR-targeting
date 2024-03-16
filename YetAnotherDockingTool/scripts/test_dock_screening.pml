import pymol.plugins
pymol.plugins.preferences['ADT_PYTHONSH'] = '/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/bin/pythonsh'
pymol.plugins.preferences['ADT_PREPARE_LIGAND4'] = '/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py'
pymol.plugins.preferences['ADT_PREPARE_RECEPTOR4'] = '/home/peu/Apps/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py'

from yadt import *
cd tests
load ace_example.pdb
dock screening, ligands_file=test.smi, \
    target=polymer.protein, \
    program=qvinaw, wdir=.
