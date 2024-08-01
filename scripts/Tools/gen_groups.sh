cd ../../PDB/dimers/
all=$(ls -1 *.pdb | sed -e 's/\..*$//')
cd ../dimers_sim/sim_1/
for f in $all; do
if [ -f $f'.pdb' ]; then
continue
else
echo $f
fi
done 
