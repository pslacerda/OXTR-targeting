from pymol import cmd as pm
import numpy as np


@pm.extend
def jaccard_similarity(sel1, sel2):
    s1 = set()
    s2 = set()
    for atom in pm.get_model(sel1).atom:
        s1.add((atom.chain, atom.resi))
    for atom in pm.get_model(sel2).atom:
        s2.add((atom.chain, atom.resi))
    ret = len(s1.intersection(s2)) / len(s1.union(s2))
    print(ret)
    return ret


@pm.extend
def get_box(box_sel, box_margin):
    """
    Create a box (useful for VINA or DOCKTHOR).
    USAGE:
        get_box sel, margin
    EXAMPLE:
        get_box resn CU and chain A, 5

    """
    box_margin = int(box_margin)

    box_coords = pm.get_coords(box_sel)

    max = np.max(box_coords, axis=0) + box_margin
    min = np.min(box_coords, axis=0) - box_margin

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
    print("Size:", size_x, size_y, size_z)
    print("Center:", center_x, center_y, center_z)
