#!/usr/bin/python3

###############################################################################
# PATHS AND DICTIONARIES
###############################################################################
import glob
import numpy as np
from pymol import cmd
import pymol
from datetime import datetime
from time import sleep
import __main__
import os
import MDAnalysis as mda
import sys

path = os.path.dirname(os.path.abspath("__file__"))
path_dimer = "/home/adrian/Documents/GitLab/web_dimers/static/PDB/dimers/"
path_dimer_lz = "/home/adrian/Documents/GitLab/web_dimers/static/PDB/dimers_lz/"
path_dimer_dis = f"/home/adrian/Documents/GitLab/web_dimers/static/PDB/dimers_sim/sim_{sys.argv[1]}/"
if sys.argv[1] == 0:
    path_dimer_dis = path_dimer
if sys.argv[1] == -1:
    path_dimer_dis = path_dimer_lz
path_dimer_ref = "/home/adrian/Documents/GitLab/web_dimers/static/PDB/"
path_results = "/home/adrian/Documents/GitLab/gpcr_dimers/RESULTS/"
d_change = {"A": "Q", "B": "R", "C": "S", "D": "T", "E": "U", "F": "V"}

###############################################################################
# IMPORTS
###############################################################################

__main__.pymol_argv = [
    "pymol",
    "-A3",
]  # https://pymolwiki.org/index.php/Command_Line_Options

###############################################################################
# FUNCTIONS
###############################################################################


def get_files(path, ext):
    """Obtain a list that contain the name of all pdb file in a directory"""
    l_files = []
    folder = glob.glob(path + "*" + ext)
    for path_filename in folder:
        path_filename = path_filename.split(
            "/"
        )  # ['', 'home', 'adrian', 'Escritorio', 'TFG', 'dimerInteractions', '2Z73_B_OPSD_TODPA.pse']
        filename = path_filename[-1]
        filename = filename.split(".")  # 4WW3_A_OPSD_TODPA.pse
        filename = filename[0]  # 4WW3_A_OPSD_TODPA
        l_files.append(filename)
    return sorted(l_files)


def get_z_values(files):
    m_pdbs = []
    for f in files:
        name = f.split("/")
        name = name[-1]
        u = mda.Universe(f)
        prot_atom = u.select_atoms("protein")
        coords = prot_atom.positions
        print(name, u, abs(coords[:, 2].max()) + abs(coords[:, 2].min()))
        if abs(coords[:, 2].max()) + abs(coords[:, 2].min()) > 100:
            m_pdbs.append(f)
    return m_pdbs


def open_file(path, filename):
    """Read a file and save it on a variable"""
    f_in = open(path + filename, "r")
    f_r = f_in.readlines()  # All info here
    f_in.close()
    return f_r


def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print("> Created directory: " + path)


def load_pymol(path_pdb, prot, ext):
    # Open reference protein
    cmd.load(path_pdb + prot + ext, quiet=0)
    print("> Loading " + prot + ext)
    sleep(0.5)


def start_pymol():
    # Open pymol
    print("Open PyMol!")
    pymol.finish_launching()


def quit_pymol():
    # Close pymol
    print("> Closing PyMol!")
    pymol.cmd.quit()


def com_pymol(prot, chain):
    # Get center of mass (both chains, angle of reference)
    ref_center = cmd.centerofmass(
        prot + """ and not resname UNX and chain """ + chain
    )  # List of 3 floats
    print("> Center of mass on " + str(ref_center))
    sleep(0.5)
    return ref_center


def pseudoatom_pymol(name, coordinates, reference, indicator):
    # Create pseudoatoms of reference
    if reference == "TRUE":
        cmd.do(
            """pseudoatom """
            + name
            + """_refcenter, pos ="""
            + str(coordinates)
            + """ """
        )
        sleep(0.5)
    else:
        cmd.do(
            """pseudoatom """
            + name
            + """_center"""
            + indicator
            + """, pos ="""
            + str(coordinates)
            + """ """
        )
        sleep(0.5)
    print("> Pseudoatom created on " + str(coordinates))


def perform_align(prot_a, chain_a, prot_b, chain_b):
    # Superposition of two chain of a dimer
    cmd.do(f"select struc, {prot_b} and chain {chain_b}")
    cmd.do(f"select ref, {prot_a}_ref and chain {chain_a}")
    super_rmsd = cmd.super("struc", "ref")
    sleep(0.5)
    align_rmsd = cmd.align("struc", "ref")
    sleep(0.5)
    if super_rmsd[0] > align_rmsd[0]:
        print("> Chains aligned!")
    else:
        super_rmsd = cmd.super("struc", "ref")
        sleep(0.5)
        print("> Chains superposed!")


def color_pymol(color, name):
    # Color these pseudoatoms
    cmd.color(color, name)


def get_vector(coord_a, coord_b):
    coord_a = np.array(coord_a)
    coord_b = np.array(coord_b)
    vector = (coord_b - coord_a) / np.linalg.norm(coord_b - coord_a)
    return vector


def get_angle_pymol(path_pdb, prot_ref, prot, chain_ref, chains, vector_ref, indicator):
    # Get the angle of three points
    name = prot + "_" + str(chains[0] + chains[1])
    cmd.set_name(prot, name)
    prot = prot + "_" + str(chains[0] + chains[1])
    name = prot.split("_")
    perform_align(prot_ref, chain_ref, prot, chains[0])
    center_1 = com_pymol(prot, chains[0])
    pseudoatom_pymol(prot + "_" + chains[0], center_1, "FALSE", indicator)
    sleep(1)
    color_pymol(2, prot + "_" + chains[0] + """_center""" + indicator)
    cmd.do(f'label {prot}_{chains[0]}, "{name[0]}_{name[1]}_{name[2]}"')
    center_2 = com_pymol(prot, chains[1])
    pseudoatom_pymol(prot + "_" + chains[1], center_2, "FALSE", indicator)
    sleep(1)
    color_pymol(2, prot + "_" + chains[1] + """_center""" + indicator)
    cmd.do(f'label {prot}_{chains[1]}, "{name[0]}_{name[1]}_{name[2]}" ')
    vector_nref = get_vector(center_1, center_2)
    angle = np.arccos(np.clip(np.dot(vector_ref, vector_nref), -1.0, 1.0))
    angle = np.degrees(angle)
    sleep(0.5)
    return (angle, center_1, center_2)


###############################################################################
# MAIN PROGRAM
###############################################################################
###############################################################################
# GET DISTRIBUTION INTERFACES
##############################################################################
if __name__ == "__main__":
    # REFERENCE PROTEIN
    print("> Start:" + str(datetime.now()))
    start_pymol()
    prot_ref = "reference_tm"
    load_pymol(path_dimer_ref, prot_ref, ".pdb")
    sleep(0.5)
    # Orient the structure on coords [0,0,0]
    cmd.set_name(prot_ref, prot_ref + "_ref")  # To differ from other objects
    chains_ref = cmd.get_chains(prot_ref)  # AB
    ref_center_a = com_pymol(prot_ref, chains_ref[0])  # A
    ref_center_b = com_pymol(prot_ref, chains_ref[1])  # B
    cmd.do(
        """alter_state 1, ( """ + prot_ref + """_ref), x=x-""" + str(ref_center_a[0])
    )
    cmd.do(
        """alter_state 1, ( """ + prot_ref + """_ref), y=y-""" + str(ref_center_a[1])
    )
    cmd.do(
        """alter_state 1, ( """ + prot_ref + """_ref), z=z-""" + str(ref_center_a[2])
    )
    # Generate ref points (com)
    pseudoatom_pymol(prot_ref + "_" + chains_ref[0], ref_center_a, "TRUE", "")
    color_pymol(4, prot_ref + "_" + chains_ref[0] + """_refcenter""")
    pseudoatom_pymol(prot_ref + "_" + chains_ref[1], ref_center_b, "TRUE", "")
    color_pymol(3, prot_ref + "_" + chains_ref[1] + """_refcenter""")
    vector_ref = get_vector(ref_center_a, ref_center_b)
    #################################################################################
    # PROTEIN
    # Open protein that I want compare
    prots = get_files(path_dimer_dis, ".pdb")
    prot_angle = []
    symmetry = open(f"{path_results}symmetry.txt", "w")
    symmetry.writelines("DIMER\tSYMMETRY\n")
    j = 1
    for prot in prots:
        print("> Protein: " + prot)
        # cord_dimers.append(prot)
        print(f"> Process: {j *100 / len(prots)} %")
        load_pymol(path_dimer_dis, prot, ".pdb")
        sleep(0.5)
        chains_1 = cmd.get_chains(prot)  # AC
        angle_1, center_1a, center_1b = get_angle_pymol(
            path_dimer_dis, prot_ref, prot, chains_ref[0], chains_1, vector_ref, str(1)
        )  # B
        pymol.cmd.do(
            """delete """
            + prot
            + "_"
            + str(chains_1[0] + chains_1[1])
            + "_"
            + chains_1[0]
            + """_center1"""
        )  # Point on the ref. one
        load_pymol(path_dimer_dis, prot, ".pdb")
        sleep(0.5)
        chains_2 = [chains_1[1], chains_1[0]]
        angle_2, center_2a, center_2b = get_angle_pymol(
            path_dimer_dis, prot_ref, prot, chains_ref[0], chains_2, vector_ref, str(2)
        )
        pymol.cmd.do(
            """delete """
            + prot
            + "_"
            + str(chains_2[0] + chains_2[1])
            + "_"
            + chains_2[0]
            + """_center2"""
        )
        cmd.center(prot_ref + "_ref")
        condition = abs(angle_1 - angle_2) < 5
        if condition:
            pymol.cmd.do("""delete """ + prot + "_" + str(chains_2[0] + chains_2[1]))
            pymol.cmd.do(
                """delete """
                + prot
                + "_"
                + str(chains_2[0] + chains_2[1])
                + "_"
                + chains_2[1]
                + """_center2"""
            )
            sleep(0.5)
            pa_1 = (prot + "_" + str(chains_1[0] + chains_1[1]), angle_1)
            if not pa_1 in prot_angle:
                prot_angle.append(pa_1)
                symmetry.writelines(f"{prot}\tYES\n")
        else:
            name1 = (
                prot
                + "_"
                + str(chains_1[0] + chains_1[1])
                + "_"
                + chains_1[1]
                + """_center1"""
            )
            color_pymol(5, name1)
            cmd.set_name(name1, name1 + "_asy")
            name2 = (
                prot
                + "_"
                + str(chains_2[0] + chains_2[1])
                + "_"
                + chains_2[1]
                + """_center2"""
            )
            color_pymol(5, name2)
            cmd.set_name(name2, name2 + "_asy")
            pa_1 = (prot + "_" + str(chains_1[0] + chains_1[1]), angle_1)
            pa_2 = (prot + "_" + str(chains_2[0] + chains_2[1]), angle_2)
            prot_angle.append(pa_1)
            prot_angle.append(pa_2)
            symmetry.writelines(f"{prot}\tNO\n")
        j += 1
        cmd.do("disable all")
    ##cmd.do('enable *center*')
    ##cmd.do('show spheres, *center*')
    ##cmd.do('set label_position, (0,2,0)')

    cmd.save(path_dimer_dis + "distributions.pse")
    cmd.do("delete not *center*")
    cmd.save(path_dimer_dis + "nostructures.pse")
    symmetry.close()

    coordinates = open(f"{path_dimer_dis}coordinates.txt", "w")
    coordinates.writelines(f"DIMER X Y Z\n")
    load_pymol(path_dimer_dis, "nostructures", ".pse")
    objects = cmd.get_object_list("(all)")
    del objects[0:2]
    for obj in objects:
        cord = cmd.get_coords(obj, 1)
        name = obj.split("_")
        name = "_".join(name[0:4])
        cord = str(cord)
        cord = cord.split("]")
        cord = cord[0].split("[")
        cord = " ".join(cord[2:])
        coordinates.writelines(f"{name} {cord}\n")
    coordinates.close()

    print("> End:" + str(datetime.now()))
    quit_pymol()
#
