#!/usr/bin/env python3

description = '''Makes a fast and dirty coarse grain-to-atomistic conversion of a protein

'''

import MDAnalysis
import os
import sys
import time


# Requires pulchra: http://cssb.biology.gatech.edu/skolnick/files/PULCHRA/index.html
pulchra_path = '/home/adrian/software/pulchra304/bin/linux/'
path_pdb = '/home/adrian/Documents/bender/SIMULATIONS_5/'
path_output = '/home/adrian/Documents/bender/PDB/dimer_sim/sim_5/'

names = os.listdir(path_pdb)
for name in names:
    file = f'2_{name}_p_center.pdb'
    print (file)
    u = MDAnalysis.Universe(f'{path_pdb}{name}/{file}')
    prot = u.select_atoms('protein and (name BB or name SC1)')
    names = [name.replace('BB', 'CA') for name in prot.names]
    names = [name.replace('SC1', 'CB') for name in names]
    prot.names = names
    prot.write(f'{path_pdb}{name}/{name}_AT.pdb')
    os.system(f'{pulchra_path}pulchra {path_pdb}{name}/{name}_AT.pdb')
    os.system(f'mv {path_pdb}{name}/{name}_AT.rebuilt.pdb {path_output}{name}.pdb')

