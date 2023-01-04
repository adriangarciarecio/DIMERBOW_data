cd ../SIMULATIONS_4/

for f in *; do
cd "$f" 
echo "RMSD of $f"

#GET ROTATE TRAJECTORY ALIGN TO PROT A
gmx trjconv -f center.xtc -s *_p.tpr -n fit.ndx -fit rot+trans  -o rotate.xtc -center << EOF
1 16 0
q     
EOF

cd .. 
done
