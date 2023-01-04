#!/bin/bash -l
#Load the module to fix less atoms
#module load pdb2pqr-2.1.0
#module load gromacs/2018.5
#INDICATE THE SIMULATION
sim=$2
#prepare directories for each protein that we will simulate
ls ../PDB/dimers/ > directories.txt
python3 Tools/gen_directories.py directories.txt $sim
rm directories.txt

cd ../SIMULATIONS_$sim/
##Obtain the list of directories to automatize the process
## for loop starts
value=$(<"../SCRIPTS/MISC/$1")
for f in $value; do
echo "Loading the dimer $f" 
cd $f
rm *.gro
rm *.tpr
echo "Preparing the TER pdb"
grep -i ATOM $f.pdb > $f"_nohet.pdb"

../../SCRIPTS/Tools/a_pdbAddTER.py $f"_nohet.pdb" $f"_ter.pdb"

pdb2pqr $f"_ter.pdb" ${f}_amb.pdb --ff=amber --chain

echo "Loading $f preparation" 

## prepare the protein in Martini
echo "Preparing $f pre-simulation step (1/4)"
#python ../../SCRIPTS/martinize.py -f "$f"_amb.pdb -o $f.top -x $f-prot.pdb -dssp dssp -p backbone -ff elnedyn22
(python ../../SCRIPTS/Tools/martinize.py -f "$f"_amb.pdb -o $f.top -x $f-prot.pdb -dssp dssp -p backbone -elastic -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0 > martini_stdout.log) >& martini_stderr.log

## insert in the lipid bilayer
echo "Preparing $f pre-simulation step (2/4)"
#python2 ../../SCRIPTS/Tools/insane.py -f $f-prot.pdb -l POPC:1 -salt 0.15 -x 15 -y 15 -z 12 -pbc cubic -sol W -o "$f"_ini.gro &> molecules.txt
python2 ../../SCRIPTS/Tools/insane.py -f $f-prot.pdb -l POPC:1 -salt 0.15 -x 15 -y 15 -z 16 -pbc cubic -sol W -o "$f"_ini.gro &> molecules.txt


## retain output of insane and copy it to dim.top
echo "Preparing $f pre-simulation step (3/4)"
python3 ../../SCRIPTS/Tools/write_martini.py $f.top

## create index.ndx
echo "Preparing $f pre-simulation step (4/4)"
gmx make_ndx -f "$f"_ini.gro -o index.ndx << EOF
14|15|16
q
EOF

## Minimization.
## equi6.0 - soft-core minimization
gmx grompp -f e0.mdp -o "$f"_e0.tpr -c "$f"_ini.gro -p $f.top -n index.ndx -r "$f"_ini.gro
gromacs "$f"_e0 hulk

##### Wait for output
while [ ! -f $f'_e0_md.gro' ] ; do sleep 5; done
#####

## equi6.1
gmx grompp -f e1.mdp -o "$f"_e1.tpr -c "$f"_e0_md.gro -p $f.top -n index.ndx -r "$f"_e0_md.gro
gromacs "$f"_e1 hulk

##### Wait for output
while [ ! -f $f'_e1_md.gro' ] ; do sleep 5; done
#####RH1_ADRB2_1_ste.o97276

## equi6.2
gmx grompp -f e2.mdp -o "$f"_e2.tpr -c "$f"_e1_md.gro -p $f.top -n index.ndx -r "$f"_e1_md.gro
gromacs "$f"_e2 hulk

##### Wait for output
while [ ! -f $f'_e2_md.gro' ] ; do sleep 5; done
#####

## equi6.3
gmx grompp -f e3.mdp -o "$f"_e3.tpr -c "$f"_e2_md.gro -p $f.top -n index.ndx -r "$f"_e2_md.gro
gromacs "$f"_e3 hulk

##### Wait for output
while [ ! -f $f'_e3_md.gro' ] ; do sleep 5; done
#####

## equi6.4
gmx grompp -f e4.mdp -o "$f"_e4.tpr -c "$f"_e3_md.gro -p $f.top -n index.ndx -r "$f"_e3_md.gro
gromacs "$f"_e4 hulk

##### Wait for output
while [ ! -f $f'_e4_md.gro' ] ; do sleep 5; done
#####

## equi6.5
gmx grompp -f e5.mdp -o "$f"_e5.tpr -c "$f"_e4_md.gro -p $f.top -n index.ndx -r "$f"_e4_md.gro
gromacs "$f"_e5 hulk

##### Wait for output
while [ ! -f $f'_e5_md.gro' ] ; do sleep 5; done
#####

## equi6.6
gmx grompp -f e6.mdp -o "$f"_e6.tpr -c "$f"_e5_md.gro -p $f.top -n index.ndx -r "$f"_e5_md.gro
gromacs "$f"_e6 hulk

##### Wait for output
while [ ! -f $f'_e6_md.gro' ] ; do sleep 5; done
#####

##step 6.6bis
gmx grompp -f e6b.mdp -o "$f"_e6b.tpr -c "$f"_e6_md.gro -p $f.top -n index.ndx -r "$f"_e6_md.gro
gromacs "$f"_e6b hulk

##### Wait for output
while [ ! -f $f'_e6b_md.gro' ] ; do sleep 5; done
#####

## Production (step 6.7)
gmx grompp -f p.mdp -o "$f"_p.tpr -c "$f"_e6b_md.gro -p $f.top -n index.ndx -r "$f"_e6b_md.gro
gromacs "$f"_p hulk

##### Wait for output
while [ ! -f $f'_p_md.gro' ] ; do sleep 5; done
#####

rm "#"*

#Convert .gro to .pdb
#gmx editconf -f "$f"_p_md.gro -o "$f"_p_md.pdb

#Search the atoms numbers of each chain (start-end)
#end1=1
#echo "RMSD of $f"
#tail -n +3 "$f"_p_md.gro | while read LINE; do 
#res=$(echo $LINE | cut -f1 -d " ")
#act=${res:0:-3} 
#if ((act > end1)); then
#end1=$act
#fi
#if ((act < end1)); then
#start2=$(echo $LINE | cut -d " " -f 3)
#end1=$((start2-1))
#line=$(grep "POPC" "$f"_p_md.gro -m 1)
#array=$(echo $line | cut -d " " -f 3)
#end2=$((array[0]-1))

#Make fit index
#gmx make_ndx -f *_p_md.gro -o fit.ndx << EOF
#a 1-$end1
#name 16 chain_a
#a $start2-$end2
#name 17 chain_b
#q
#EOF
#break 
#fi 
#done

##### Wait for output
#while [ ! -f 'fit.ndx' ] ; do sleep 5; done
#####

#RMSD 1
#gmx rms -f *_p.xtc -s *_e0.tpr -n fit.ndx << EOF
#16 17
#EOF

##### Wait for output
#while [ ! -f 'rmsd.xvg' ] ; do sleep 5; done
#####

#gmx trjconv -f *_p.xtc -s *_p.tpr -n fit.ndx -pbc cluster -o pre_center.xtc -center << EOF
#1 13 0
#q
#EOF

##### Wait for output
#while [ ! -f 'pre_center.xtc' ] ; do sleep 10; done
#####

#gmx trjconv -f pre_center.xtc -s *_p.tpr -n fit.ndx -pbc res -o center.xtc -center << EOF
#1 0
#q
#EOF

##### Wait for output
#while [ ! -f 'center.xtc' ] ; do sleep 10; done
#####

#gmx trjconv -s *_p.tpr -f center.xtc -o "$f"_p_center.gro -dump 600000 << EOF
#0
#q
#EOF

##### Wait for output
#while [ ! -f $f'_p_center.gro' ] ; do sleep 5; done
#####

#Convert .gro to .pdb
#gmx editconf -f "$f"_p_center.gro -o "$f"_p_center.pdb

#RMSD 2
#gmx rms -f center.xtc -s *_e0.tpr -n fit.ndx -o rmsd2.xvg << EOF
#16 17
#EOF

#Return to the previous dir (SIMULATIONS)
cd ..

## for loop ends
done
