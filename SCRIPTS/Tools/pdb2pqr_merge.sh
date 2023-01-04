#!/bin/bash -l

# load module to add Hs
module load pdb2pqr-2.1.0

pdb2pqr /people/common/add_water/protein-wat.pdb /people/common/add_water/protein-wat_Hs.pdb --ff=charmm

python2 /people/common/add_water/scripts/merge_prot_wHsv2.py
