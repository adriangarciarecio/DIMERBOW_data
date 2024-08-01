#module load gromacs/2018.5
export GMXLIB=/opt/soft/gromacs-2018.3-source/share/top/amber99sb-ildn.ff
sim=$1
path_out="sim_$sim"
cd ../../SIMULATIONS_$sim/
value="$(ls)"
for f in $value; do
while [ ! -f ../PDB/dimers_sim/"$path_out"/"$f".pdb ] ; do
echo "Converting dimer $f to atomistic model" 
cd $f 
cp ../../SCRIPTS/Tools/backward.py . 
cp ../../SCRIPTS/Tools/initram-v5.sh . 
cp -r ../../SCRIPTS/Tools/Mapping . 
rm *#*
gmx trjconv -s "$f"_p.tpr -f 2_"$f"_p_center.gro -o "$f"_protein.gro << EOF
1
EOF
gmx pdb2gmx -f "$f"_amb.pdb -p topol_at.top -ignh -water tip3p -ff amber99sb-ildn
bash ../../SCRIPTS/Tools/initram-v5.sh -f "$f"_protein.gro -o "$f"_aa_amber.gro -to amber -p topol_at.top
gmx editconf -f "$f"_aa_amber.gro -o "$f"_aa_amber.pdb
echo "Adding chains to dimer $f"
python3 ../../SCRIPTS/Tools/add_chain_letters.py "$path_out" "$f"
cd ..
done
echo "$f converted!"
done 
