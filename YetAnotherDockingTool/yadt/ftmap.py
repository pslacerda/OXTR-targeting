from itertools import combinations
from types import SimpleNamespace

import numpy as np
from scipy.spatial import distance_matrix
from pymol import cmd as pm
from matplotlib import pyplot as plt


__all__ = ["load_ftmap", "fo", "dc", "dce"]


def get_kozakov2015_class(combination):
    s0 = combination.clusters[0].strength
    cd = combination.center_center
    md = combination.max_dist
    return (
        "D"
        if (s0 >= 16 and cd < 8 and md >= 10)
        else (
            "Dl"
            if (s0 >= 16 and cd >= 8 and md >= 10)
            else (
                "Ds"
                if (s0 >= 16 and cd < 8 and md < 10)
                else (
                    "B"
                    if (s0 < 16 and cd < 8 and md >= 10)
                    else (
                        "Bl"
                        if (s0 < 16 and cd >= 8 and md >= 10)
                        else "Bs" if (s0 < 16 and cd < 8 and md < 10) else None
                    )
                )
            )
        )
    )


def collect_clusters(prefix="crosscluster"):
    for obj in pm.get_object_list():
        if obj.startswith(f"{prefix}."):
            coords = pm.get_coords(obj)
            try:
                _, _, s, _ = obj.split(".", maxsplit=4)
            except ValueError:
                _, _, s = obj.split(".", maxsplit=3)
            yield SimpleNamespace(selection=obj, strength=int(s), coords=coords)


def collect_kozakov2015(max_length=3, prefix="crosscluster", max_combinations=50):
    idx = 0
    clusters = list(collect_clusters(prefix))
    for length in range(1, max_length + 1):
        for combination in list(combinations(clusters, length)):
            coords = np.concatenate([c.coords for c in combination])
            dist = []
            for cluster1 in combination:
                for cluster2 in combination:
                    if cluster1.selection == cluster2.selection:
                        dist.append(0)
                    else:
                        dist.append(
                            distance_matrix(cluster1.coords, cluster2.coords).min()
                        )
            comb = SimpleNamespace(
                selection=" or ".join(c.selection for c in combination),
                clusters=combination,
                kozakov_class=None,
                strength=sum(c.strength for c in combination),
                strength0=combination[0].strength,
                center_center=np.min(dist),
                max_dist=distance_matrix(coords, coords).max(),
            )
            comb.kozakov_class = get_kozakov2015_class(comb)
            yield comb
            idx += 1
            if idx >= max_combinations:
                break


@pm.extend
def load_ftmap(
    filename, group="FTMap", max_length=3, origin="ftmap", max_combinations=50
):
    """
    Load a FTMap PDB file and classify hotspot ensembles in accordance to
    Kozakov et al. (2015).
    https://doi.org/10.1021/acs.jmedchem.5b00586

    OPTIONS
        filename    mapping PDB file.
        group       optional group name to put objects in.
        max_length  the maximum number of consensus sites to consider.

    EXAMPLES
        load_ftmap ace_example.pdb
        load_ftmap ace_example.pdb, MyProtein
    """
    if origin == "ftmap":
        prefix = "crosscluster"
    elif origin == "atlas":
        prefix = "consensus"
    else:
        raise Exception(f"Unknown origin: {origin}")

    pm.load(filename)

    combinations = list(
        collect_kozakov2015(
            max_length=max_length, prefix=prefix, max_combinations=max_combinations
        )
    )
    idx = 0
    for comb in combinations:
        if comb.kozakov_class:
            new_name = f"{group}.{comb.kozakov_class}.{idx}"
            pm.create(new_name, comb.selection)
            pm.group(f"{group}.{comb.kozakov_class}", new_name)
            comb.selection = new_name
            idx += 1

    pm.group(group, f"{group}.D or {group}.Ds or {group}.B or {group}.Bs")

    pm.set_name("protein", f"{group}.protein")
    pm.group(group, f"{group}.protein")

    for cluster in collect_clusters():
        new_name = f"{group}.cluster.{cluster.strength}_{idx}"
        pm.create(new_name, cluster.selection)
        pm.group(f"{group}.cluster", new_name)
        idx += 1

    pm.delete(f"{prefix}.*")

    pm.show("mesh", f"{group}.D* or {group}.B*")
    pm.color("red", f"{group}.D*")
    pm.color("salmon", f"{group}.B*")

    pm.disable(f"{group}.clusters")

    pm.orient()
    return combinations


@pm.extend
def fo(sel1, sel2, radius=2, state1=1, state2=1, verbose=1):
    """
    Compute the fractional overlap of sel1 respective to sel2.
        FO = Nc/Nt

    Nc is the number of atoms of sel1 in contact with sel2. Nt is the number of atoms
    of sel1.

    Hydrogen atoms are ignored.

    If the contact radius is 0 then the VdW radii will be used.

    The states are for select a single state from a multi-state objects.

    OPTIONS:
        sel1    ligand object.
        sel2    hotspot object.
        radius  the radius so two atoms are in contact (default: 2).
        state1  state of sel1.
        state2  state of sel2.

    EXAMPLES:
        get_fractional_overlap ref_lig, ftmap1234.D.003
        get_fractional_overlap ref_lig, ftmap1234.CS.000_016
    """
    atoms1 = pm.get_coords(f"({sel1}) and not elem H", state=state1)
    atoms2 = pm.get_coords(f"({sel2}) and not elem H", state=state2)
    dist = distance_matrix(atoms1, atoms2) - float(radius) <= 0
    num_contacts = np.sum(np.any(dist, axis=1))
    total_atoms = len(atoms1)
    fo_ = num_contacts / total_atoms
    if int(verbose):
        print(f"FO: {fo_:.2f}")
    return fo_


@pm.extend
def dc(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1):
    """
    Compute the Density Correlation according to:
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3264775/

    USAGE:
        dc(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1)

    sel1 and sel2 are respectively the molecule and hotspot; state1 and state2 are
    the optional corresponding states (default to first state both). The threshold
    distance can be changed with dist (default to 1.25).

    verbose is the standard boolean API option to define verbosity.
    """
    xyz1 = pm.get_coords(f"({sel1}) and not elem H", state1)
    xyz2 = pm.get_coords(f"({sel2}) and not elem H", state2)
    dc_ = (distance_matrix(xyz1, xyz2) < float(dist)).any(1).sum()
    if int(verbose):
        print(f"DC: {dc_:.2f}")
    return dc_


@pm.extend
def dce(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1):
    """
    Compute the Density Correlation Efficiency according to:
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3264775/

    USAGE:
        dce(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1)

    All parameters are the same for dc. The reference selection
    wich atoms are counted is sel1.
    """
    dc_ = dc(sel1, sel2, state1, state2, dist, verbose=False)
    dce_ = dc_ / pm.count_atoms(f"({sel1}) and not elem H")
    if int(verbose):
        print(f"DCE: {dce_:.2f}")
    return dce_
