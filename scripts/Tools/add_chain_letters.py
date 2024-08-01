import sys  
# Adds chain letters based on the residue number
# Starts assigning chain A and shifts to B when the numbering restarts
prev_resid = 0
with open(f"{sys.argv[2]}_aa_amber.pdb") as input_pdb:
    with open(f"../../PDB/dimers_sim/{sys.argv[1]}/{sys.argv[2]}.pdb", "w") as output_pdb:
        for line in input_pdb:
            if not line.startswith("ATOM"):
                output_pdb.write(line)
            else:
                resid = int(line[20:26].strip())
                if resid >= prev_resid:
                    output_pdb.write(line[:20] + " A " + line[23:])
                    prev_resid = resid
                else:
                    output_pdb.write(line[:20] + " B " + line[23:])
                    prev_resid = 99999999

