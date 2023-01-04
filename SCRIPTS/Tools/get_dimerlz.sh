cd ../../PDB
mkdir dimers_lz
cd dimers 
declare -a dimers=($(ls *.pdb))
cd ../dimers_nlz
declare -a dimers_lz=($(ls *.pdb))
#for d in $dimers_lz; do
#cp ../dimers_origin/$d ../dimers_lz/"${dimers[i]}"
for ((i=0;i<${#dimers_lz[@]};++i)); do
  printf "%s is in %s\n" "${dimers_lz[i]}" "${dimers[i]}"
  cp ../dimers_origin/"${dimers_lz[i]}" ../dimers_lz/"${dimers[i]}"
done



