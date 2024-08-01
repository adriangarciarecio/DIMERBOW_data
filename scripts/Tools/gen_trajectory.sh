cd ../../SIMULATIONS_$2/

value=$(<"../SCRIPTS/MISC/$1")
for f in $value; do
cd "$f" 
echo "RMSD of $f"

rm fit.ndx 
rm \#*
rm step*.pdb

#Search the atoms numbers of each chain (start-end)
end1=1
echo "RMSD of $f"
tail -n +3 "$f"_p_md.gro | while read LINE; do
res=$(echo $LINE | cut -f1 -d " ")
act=${res:0:-3}
if ((act > end1)); then
end1=$act
fi
if ((act < end1)); then
start2=$(echo $LINE | cut -d " " -f 3)
end1=$((start2-1))
line=$(grep "POPC" "$f"_p_md.gro -m 1)
array=$(echo $line | cut -d " " -f 3)
end2=$((array[0]-1))

echo $res
echo $act
echo $end1
echo $start2
echo $array
echo $end2 

#Make fit index
gmx make_ndx -f "$f"_p_md.gro -o fit.ndx << EOF
a 1-$end1
name 16 chain_a
a $start2-$end2
name 17 chain_b
q
EOF
break
fi
done

##### Wait for output
while [ ! -f 'fit.ndx' ] ; do sleep 5; done
#####

#RMSD 1
gmx rms -f *_p.xtc -s *_e0.tpr -n fit.ndx << EOF
16 17
EOF

#CENTER THE SIMULATION TO MEMBRANE GET BETTER SIMULATIONS
gmx trjconv -f *_p.xtc -s *_p.tpr -n fit.ndx -pbc cluster -o center.xtc -center << EOF
1 13 0
q
EOF
 
gmx trjconv -f center.xtc -s *_p.tpr -n fit.ndx -pbc res -o center.xtc -center << EOF
1 0
q
EOF

GET ROTATE TRAJECTORY ALIGN TO PROT A
gmx trjconv -f center.xtc -s *_p.tpr -n fit.ndx -fit rot+trans -o rotate.xtc -center << EOF
1 16 0
q     
EOF


#LAST FRAME OUT
gmx trjconv -s *_p.tpr -f center.xtc -o "$f"_p_center.gro -dump 600000 << EOF
0
q
EOF

#Convert .gro to .pdb
#gmx editconf -f "$f"_p_center.gro -o "$f"_p_center.pdb

#RMSD 1
gmx rms -f center.xtc -s *_e0.tpr -n fit.ndx -o rmsd2.xvg << EOF
16 17
EOF

cd .. 
done
