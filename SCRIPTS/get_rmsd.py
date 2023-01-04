#!/usr/bin/python3
###############################################################################
# IMPORTS
###############################################################################
import __main__
import sys
import os
import re
__main__.pymol_argv = ['pymol','-A3'] # https://pymolwiki.org/index.php/Command_Line_Options
path = os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(0, path)
from gen_distribution import *
from get_interactions import get_all_gro
import MDAnalysis as mda

###############################################################################
# PATHS AND DICTIONARIES
###############################################################################
path_dimer = path + '/../SIMULATIONS_'+  sys.argv[1] +'/'
path_results = path + '/../RESULTS/'
d_change = {'A':'Q', 'B': 'R', 'C':'S', 'D':'T', 'E':'U', 'F':'V'}

###############################################################################
# FUNCTIONS
###############################################################################


###############################################################################
# MAIN PROGRAM
###############################################################################

###############################################################################
# RENAME PDBS
###############################################################################
if __name__ == "__main__":
    dimers = get_files(path_dimer, '')
    print ('> Start:' + str(datetime.now()))
    #start_pymol()
    #dimers=['2Z73_OPSD_3']
    allrmsdinfo = open(f'{path_results}dimer_rmsd_'+  sys.argv[1] +'.txt', "w")
    for j, dimer in enumerate(dimers):
        title = ['DIMER', 'FAMILY']
        info = []
        info.append(dimer)
        print (f'> Getting rmsd of {dimer}.')
        clas = dimer.split("_")
        clas = clas[1]
        info.append(clas)
        path_xvg_dimer = f'{path_dimer}{dimer}/'
        xvgfile = f'rmsd2.xvg'
        xvg_r = open_file(path_xvg_dimer, xvgfile)#Save all pdb info into a variable
        for i, x in enumerate(xvg_r):
            if not (x.startswith('@') or x.startswith('#')):
                rmsd = x.split(" ")
                rmsd = list(filter(None, rmsd))
                if float(rmsd[0]) > 500000:
                    info.append(rmsd[1][0:-1])
                    if not rmsd[0] in title:
                        title.append(f'Time_{rmsd[0]}')
        if j == 0:
            title.append('\n')
            title = "\t".join(title)
            allrmsdinfo.writelines(title)
        info.append("\n")
        info = "\t".join(info)
        allrmsdinfo.writelines(info)
    print ('> End:' + str(datetime.now()))
    allrmsdinfo.close()
