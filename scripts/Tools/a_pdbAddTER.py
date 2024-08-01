#!/usr/bin/env python3
#
# adds TER records after lipids and water
#

import sys
if len(sys.argv) <= 2:
        print(('usage: python',sys.argv[0],'input.pdb output.pdb'))
        sys.exit(-1)


pdbin = open(sys.argv[1], 'r')
pdbout = open(sys.argv[2], 'w')
reso= 0
lines=pdbin.readlines()
for i, line in enumerate(lines):
    if line.startswith('ATOM'):
        resn=line[17:20]
        resi=int(line[23:26])
        if resi!=reso and resi != reso + 1 and i != 0:
            pdbout.write('TER\n')
        reso = resi
        pdbout.write(line)
pdbout.write('TER\n')

