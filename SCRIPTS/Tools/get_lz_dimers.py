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

###############################################################################
# PATHS
###############################################################################
path_lz_name = path + "/../PDB/dimers_nlz/"
path_lz = path + "/../PDB/dimers_lz/"
path_dimers = path + "/../PDB/dimers/"
path_origin = path + "/../PDB/dimers_origin/"

###############################################################################
# MAIN PROGRAM
###############################################################################
lz_names = os.listdir(path_lz_name)
dimers = os.listdir(path_dimers)
start_pymol()
for d in dimers:
    code = d.split("_")
    code = code[0]
    load_pymol(path_dimer, d, "")
    for lz in lz_names:
        if code in lz:
            load_pymol(path_origin, lz, "")
            lz_name = lz[:-4]
            d_name = d[:-4]
            super_rmsd = cmd.super(d_name, lz_name)
            sleep(0.5)
            align_rmsd = cmd.align(d_name, lz_name)
            sleep(0.5)
            if super_rmsd[0] > align_rmsd[0]:
                print("> Chains aligned!")
                rmsd = align_rmsd[0]
            else:
                super_rmsd = cmd.super(d_name, lz_name)
                rmsd = super_rmsd[0]
                sleep(0.5)
                print("> Chains superposed!")
            print(rmsd)
            if int(rmsd) == 0:
                print(f"> {lz_name} is {d_name}!!!")
                cmd.do(f"save {path_lz}{d}, {lz_name}")
                cmd.reinitialize()
                break
            else:
                cmd.do(f"remove {lz_name}")
                continue
quit_pymol()

