import os
import itertools
import glob
import MDAnalysis as mda

path = os.path.dirname(os.path.abspath("__file__"))
pdbfile = "4QIM_SMO_61.pdb"
path_pdb = path + "/../PDB/dimers/"

u = mda.Universe(
    path_pdb + pdbfile
)  # Create a universe with the pdbfile (after select the part that it will analysed):
chains = ["A", "Q"]
prot = "protein and not (name H*) and segid " + chains[0]
sear = (
    "protein and not (name H*) and segid "
    + chains[1]
    + " and around 20 protein and not (name H*) and segid "
    + chains[0]
)
print(u, chains)
sel_prot = u.select_atoms(prot)
print(sel_prot)

protein_a_posits = sel_prot.positions
print(protein_a_posits)
near_atoms = u.select_atoms(sear)
protein_b_posits = near_atoms.positions

print(u, near_atoms)
