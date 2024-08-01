#!/usr/bin/python3
###############################################################################
# IMPORTS
###############################################################################
import __main__
import sys
import os
import re

__main__.pymol_argv = [
    "pymol",
    "-A3",
]  # https://pymolwiki.org/index.php/Command_Line_Options
path = os.path.dirname(os.path.abspath("__file__"))
sys.path.insert(0, path)
from gen_distribution import *
from get_interactions import get_all_gro
import MDAnalysis as mda
import math

###############################################################################
# PATHS AND DICTIONARIES
###############################################################################
path_dimer = path + "/../SIMULATIONS_" + sys.argv[1] + "/"
path_results = path + "/../RESULTS/"
d_change = {"A": "Q", "B": "R", "C": "S", "D": "T", "E": "U", "F": "V"}

###############################################################################
# FUNCTIONS
###############################################################################
def centroid(selection):
    model = cmd.get_model(selection)
    nAtom = len(model.atom)
    centroid = cpv.get_null()
    for a in model.atom:
        centroid = cpv.add(centroid, a.coord)
    centroid = cpv.scale(centroid, 1.0 / nAtom)
    return centroid


def get_vector(coord_a, coord_b):
    coord_a = np.array(coord_a)
    coord_b = np.array(coord_b)
    vector = (coord_b - coord_a) / np.linalg.norm(coord_b - coord_a)
    return vector


def get_angle_pymol(coord_a, coord_b, v_ref):
    vector_b = get_vector(coord_a, coord_b)
    # angle = np.arccos(np.clip(np.dot(v_ref, vector_b), -1.0, 1.0))
    cos_angle = np.dot(v_ref, vector_b) / (
        np.linalg.norm(v_ref) * np.linalg.norm(vector_b)
    )
    # print (cos_angle)
    angle = np.arccos(float(str(cos_angle)[:-1]))
    angle = np.degrees(angle)
    print(angle)
    return angle


def Euc_dist(p, q):
    """Calculates the distance between two coordinates"""
    dist = math.sqrt((q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2)
    print(dist)
    return dist


###############################################################################
# MAIN PROGRAM
###############################################################################

###############################################################################
# RENAME PDBS
###############################################################################
if __name__ == "__main__":
    dimers = get_files(path_dimer, "")
    print("> Start:" + str(datetime.now()))
    # start_pymol()
    # dimers=['5GLH_EDNRB_41']
    alldistinfo = open(f"{path_results}dimer_dist_" + sys.argv[1] + ".txt", "w")
    allangleinfo = open(f"{path_results}dimer_angles_" + sys.argv[1] + ".txt", "w")
    title = ["DIMER", "FAMILY", "TRJ_0"]
    for i in range(350, 401):
        title.append(f"TRJ_{i}")
    title.append("\n")
    title = "\t".join(title)
    allangleinfo.writelines(title)
    alldistinfo.writelines(title)
    for dimer in dimers:
        print(f"> Getting angles of {dimer}.")
        path_gro_dimer = f"{path_dimer}{dimer}/"
        try:
            grofile = f"2_{dimer}_p_center.gro"
            gro_r = open_file(
                path_gro_dimer, grofile
            )  # Save all pdb info into a variable
        except:
            grofile = f"2_{dimer}_p_md.gro"
            gro_r = open_file(
                path_gro_dimer, grofile
            )  # Save all pdb info into a variable
        wat_ions = ["NA+", "CL-"]
        protein, popc, w, ion, rest = get_all_gro(gro_r, wat_ions)
        atomnum = re.sub("[POPC]", "", popc[0][2])
        atomnum = int(atomnum) - 1
        start2 = 0
        for p in protein:
            act = int(p[0][:-3])
            if act > start2:
                start2 = act
            if act < start2:
                start2 = int(p[2])
                break
        end1 = start2 - 1
        # load_pymol(path_gro_dimer, grofile[:-4], '.gro')
        u = mda.Universe(f"{path_gro_dimer}{grofile}", f"{path_gro_dimer}rotate.xtc")
        info_angle = []
        info_dist = []
        clas = dimer.split("_")
        clas = clas[1]
        info_angle.append(dimer)
        info_angle.append(clas)
        info_dist.append(dimer)
        info_dist.append(clas)
        for i, ts in enumerate(u.trajectory):
            if i >= 350 or i == 0:
                print(f"> Reading trajectory {i} of {len(u.trajectory)}.")
                prot_a = u.atoms[: int(end1)]
                prot_b = u.atoms[int(start2) - 1 : int(atomnum)]
                cenA = prot_a.center_of_geometry()
                cenB = prot_b.center_of_geometry()
                if i == 0:
                    vect_ref = get_vector(cenA, cenB)
                angle = get_angle_pymol(cenA, cenB, vect_ref)
                dist = Euc_dist(cenA, cenB)
                info_angle.append(str(angle))
                info_dist.append(str(dist))
        info_angle.append("\n")
        info_dist.append("\n")
        info_angle = "\t".join(info_angle)
        info_dist = "\t".join(info_dist)
        allangleinfo.writelines(info_angle)
        alldistinfo.writelines(info_dist)
    print("> End:" + str(datetime.now()))
    allangleinfo.close()
    alldistinfo.close()
