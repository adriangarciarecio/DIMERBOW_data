# -*- coding: utf-8 -*-
import os
import sys

###############################################################################
# EDIT dim.top
###############################################################################
file_in = open(f'{sys.argv[1]}', 'r')
f_r = file_in.readlines() #All info here
os.system(f'cp {sys.argv[1]} copy_{sys.argv[1]}')
file_out = open(f'{sys.argv[1]}', 'w')
file_out.writelines('#include "martini_v2.2.itp"\n')
file_out.writelines('#include "martini_v2.0_ions.itp"\n')
file_out.writelines('#include "martini_v2.0_POPC_02.itp"\n \n')
for line in f_r:
    if line.startswith('#include') and 'martini' in line:
        f_r.remove(line)
    else:
        file_out.writelines(line)
file_out.writelines('\n') #At the end don't have it and the first molecule is writed at the same line

###############################################################################
# ADD MOLECULES
###############################################################################
file_in = open(f'molecules.txt', 'r')
f_r = file_in.readlines() #All info here
for line in f_r:
    if not line.startswith(';'):
        file_out.writelines(line)

file_out.close()