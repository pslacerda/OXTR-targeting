from yadt import *

cd tests
load_ftmap ace_example.pdb, group=ACE
box *D.*, 5, 1
dock analyze, \
    group=ACE, \
    target=polymer.protein, \
    program=qvinaw, wdir=., \
    max_ligands=4
delete target
