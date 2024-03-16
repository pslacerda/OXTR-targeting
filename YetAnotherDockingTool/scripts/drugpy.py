#!/usr/bin/env python3
from subprocess import check_call
from shutil import which

conda = "conda"
if which(conda) is None:
    conda = "micromamba"

def require(cmd, repo, *pkg):
    try:
        assert which(cmd)
    except:
        check_call([conda, "install", "-y", "-c", repo, *pkg])


require("vina", "bioconda", "autodock-vina")
require("obabel", "conda-forge", "openbabel", "qvina", "scipy", "matplotlib")

try:
    from yadt.docking import *
    from yadt.ftmap import *
except:
    check_call(["pip", "install", "yadt"])
    from yadt.docking import *
    from yadt.ftmap import *

    print("#=============================#")
    print("# yadt::YetAnotherDockingTool #")
    print("#=============================#")
