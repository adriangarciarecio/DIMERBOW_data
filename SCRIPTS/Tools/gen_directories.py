#!/usr/bin/python3

###############################################################################
# PATHS AND DICTIONARIES
###############################################################################
import os
import sys
path = os.path.dirname(os.path.abspath('__file__'))
path_simulations = path + '/..' + '/SIMULATIONS_'+ sys.argv[2] +'/'
path_dimer_dis = path + '/..' + '/PDB/dimers/'

###############################################################################
# IMPORTS
###############################################################################
import shutil
import glob
sys.path.insert(0, path)

def open_file (path,filename):
    '''Read a file and save it on a variable'''
    f_in = open(path + filename, "r")
    f_r = f_in.readlines() #All info here
    f_in.close()
    return f_r

def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print ("> Created directory: " + path)

###############################################################################
# MAIN PROGRAM
###############################################################################
create_dir(path_simulations) #This directory will contain all simulations of each dimer
dimers = open_file(path,'/' + sys.argv[1])
for dimer in dimers:
    dimer = dimer.replace("\n", "").replace(".pdb", "")
    path_dimer = path_simulations + dimer + '/'
    create_dir(path_dimer)
    shutil.copy2(f"{path_dimer_dis}{dimer}.pdb", path_dimer)
    print (f'> Loading pdb file {dimer}')
    files_itp = glob.iglob(os.path.join("../GROMACS/MARTINI", "*.itp"))
    print ('> Loading itp files to each directory')
    for file in files_itp:
        if os.path.isfile(file):
            shutil.copy2(file, path_dimer)
    print ('> Loading mdp files to each directory')
    files_mdp = glob.iglob(os.path.join("../GROMACS/MARTINI", "*.mdp"))
    for file in files_mdp:
        if os.path.isfile(file):
            shutil.copy2(file, path_dimer)
