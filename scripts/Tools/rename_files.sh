for file in *.pdb
do
  mv "$file" "${file/.pdb/_opm.pdb}"
done
