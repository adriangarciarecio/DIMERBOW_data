#!/usr/bin/python3
###############################################################################
# IMPORTS
###############################################################################
import os
import sys
from gen_distribution import create_dir, open_file, get_files, start_pymol, load_pymol, quit_pymol
from pymol import cmd
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.distances
import itertools  # Combinations without repetition
import glob
import collections
from mysql.connector import connect
import requests
from urllib.request import urlretrieve
from bs4 import BeautifulSoup
import re
import pandas as pd
from sqlalchemy import create_engine
from time import sleep
from termcolor import colored
from Tools.SECRETS import *

###############################################################################
# PATHS AND DICTIONARIES
###############################################################################
path = os.path.dirname(os.path.abspath("__file__"))
sys.path.insert(0, path)
path_pse = path + "/../data/PSE/"
path_opm = path + "/../data/OPM/"
path_pdb = path + "/../data/PDB/"
path_misc = path + "/../Misc/"
path_dimer_lz = path + "/../data/PDB/dimers_origin/"
path_dimer = path + "/../data/PDB/dimers_nlz/"
path_dimer_end = path + "/../data/PDB/dimers/"
path_removed = path + "/../data/PDB/dimers_removed/"
d_change = {"A": "Q", "B": "R", "C": "S", "D": "T", "E": "U", "F": "V"}

###############################################################################
# FUNCTIONS
###############################################################################
def get_coord_lz():
    d_coord = {"default": "-1000-0+900-2000"}  # default is to remove -1000-0+900-2000
    try:
        engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')

        print(engine)
    except:
        print("could not connect")
        sys.exit(0)

    # Get all data    
    sql_query = f"SELECT pdb, segments FROM gpcrdb_info WHERE resolution < 90 and type='X-ray';"
    
    data = pd.read_sql(sql_query, engine)
    
    # Create the dict of segments
    for index, row in data.iterrows():
        d_coord[row["pdb"]] = row["segments"]
    
    return d_coord

def rem_lz(lz, lz_info, path_dimer_lz, path_dimer, cutoff_lz):
    codes = lz_info.keys()
    cmd.load(path_dimer_lz + lz + ".pdb", quiet=0)
    name = lz.split("_")
    if name[0] in codes:
        cmd.do(f"select not resi {lz_info[name[0]]} and not het")
        cmd.do("remove sele")
    # Normal case
    else:
        cmd.do(
            f'select resi {lz_info["default"]} and not het'
        )  # default is to remove -1000-0+900-2000
        cmd.do("remove sele")
    sleep(1)
    cmd.do(
        "select not byres (all and not het) around " + cutoff_lz + " + all and not het"
    )
    cmd.do("remove sele")
    cmd.save(path_dimer + lz + ".pdb")
    sleep(1)

def rem_chains(path):
    dic = {}
    try:
        engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')

        print(engine)
    except:
        print("could not connect")
        sys.exit(0)

    # Get all data    
    sql_query = f"SELECT * FROM gpcrdb_info WHERE resolution < 90 and type='X-ray diffraction';"

    data = pd.read_sql(sql_query, engine)
    
    # Get chains from database   
    for index, row in data.iterrows():
        chains = row["preferred_chain"]  
        if str(row["signalling_protein.data.entity1.chain"]) != "None":
            chains = chains + "+" + str(row["signalling_protein.data.entity1.chain"])
        if str(row["signalling_protein.data.entity2.chain"]) != "None":
            chains = chains + "+" + str(row["signalling_protein.data.entity2.chain"])
        if str(row["signalling_protein.data.entity3.chain"]) != "None":
            chains = chains + "+" + str(row["signalling_protein.data.entity3.chain"])
        dic[row["pdb_code"]] = chains     
    return dic

def get_opm(path_opm, pdb):
    link = "https://opm-assets.storage.googleapis.com/pdb/"
    urlretrieve(link + pdb + ".pdb", path_opm + pdb + ".pdb")


def reference_dict(a_pdb, b_pdb, d_ref):
    if a_pdb != b_pdb:
        if not a_pdb in d_ref.keys():  # NO
            if b_pdb in d_ref.keys():  # YES, eliminated
                d_ref[b_pdb] = f"{d_ref[b_pdb]}, {b_pdb}"  # Add himself
                d_ref[a_pdb] = d_ref[b_pdb]  # Append all the list to the not removed
                del d_ref[b_pdb]  # Remove eliminated
            else:  # NO, eliminated (first time on dict)
                d_ref[a_pdb] = b_pdb
        else:  # YES
            if not b_pdb in d_ref[a_pdb]:  # No repeat
                d_ref[a_pdb] = f"{d_ref[a_pdb]}, {b_pdb}"


###############################################################################
# GENERATE ALL PDBS
###############################################################################
cutoff_rmsd = 3
cutoff_lz = str(8)
min_dist = str(10)
option_generate = True
option_nlz = True
option_dist = True
option_super = False
option_check = False

if __name__ == "__main__":
    create_dir(path_pdb)
    create_dir(path_dimer_lz)
    create_dir(path_dimer_end)
    create_dir(path_misc)
    create_dir(path_dimer)
    create_dir(path_opm)
    create_dir(path_removed)

    # Get the all unique code present on pse
    l_psefiles = get_files(path_pse, ".pse")
    l_pdbs_codes = get_files(path_dimer_end, ".pdb")
    l_pdbs = get_files(path_dimer_lz, ".pdb")
    # l_pdbs = get_files(path_pse, ".pse")
    l_pdb_code = []
    for pdb in l_pdbs:
        pdb = pdb.split("_")
        pdb = pdb[0]
        if not pdb in l_pdb_code:
            l_pdb_code.append(pdb)
    # Get the last id:
    for pdb in l_pdbs_codes:
        pdb = pdb.split("_")
        last_id = pdb[2]
    start_pymol()
    cmd.do("set retain_order, 1")
    cmd.do("set pdb_use_ter_records, 1")  # INSERT TER LINE
    cmd.do("set ignore_pdb_segi, 1")
    sleep(1)

    # l_psefiles = ['3ODU_A_CXCR4_HUMAN', '3ODU_B_CXCR4_HUMAN', '3UON_A_ACM2_HUMAN']
    dic_chains = rem_chains(path_misc)  # Get info about real chains of the dimer
    if option_generate == True:
        j = 1
        for psefile in l_psefiles:
            pse = psefile.split("_")
            pse = pse[0]
            if not pse in l_pdb_code:
                print(f"> {j *100 / len(l_psefiles)} %")
                file_okay = glob.glob(f"{path_dimer_lz}{psefile}*.pdb")
                if file_okay == []:
                    print(f"> Getting dimers from {psefile}")
                    load_pymol(path_pse, psefile, ".pse")
                    name = psefile.split("_")
                    if name[0] in dic_chains:
                        cmd.do(
                            f"select * and not chain {dic_chains[name[0]]} and not het"
                        )
                        cmd.remove("sele")
                    l_objects = cmd.get_object_list("(all)")
                    for a, b in itertools.combinations(
                        l_objects, 2
                    ):  # Obtain all dimers NO FILTER!!
                        chain_a = cmd.get_chains(a)
                        chain_b = cmd.get_chains(b)
                        try:#KeyError: 'AAA' 
                            if chain_a[0] == chain_b[0]:
                                cmd.do(
                                    f'alter {b}, chain= "{d_change[chain_b[0]]}"'
                                )  # A --> Q'
                                cmd.save(
                                    f"{path_dimer_lz}{psefile}_{a}_{b}_{chain_a[0]}{d_change[chain_b[0]]}.pdb",
                                    f"{a} + {b}",
                                )
                                cmd.do(f'alter {b}, chain= "{chain_b[0]}"')
                            else:
                                cmd.save(
                                    f"{path_dimer_lz}{psefile}_{a}_{b}_{chain_a[0]}{chain_b[0]}.pdb",
                                    f"{a} + {b}",
                                )
                            sleep(1)
                        except:
                            continue
                    cmd.reinitialize()
                    sleep(1)
            j += 1

    ################################################################################
    ## REMOVE LISOZIM(FILTER I)
    #################################################################################
    # start_pymol()
    lz_files = get_files(path_dimer_lz, ".pdb")
    dimer_files = get_files(path_dimer, ".pdb")
    # lz_files = ['4UG2_A_AA2AR_HUMAN_4UG2A_symB03-10000_AB', '4UG2_B_AA2AR_HUMAN_4UG2B_symA03000000_BA']
    if option_nlz == True:
        lz_info = get_coord_lz()  # Get info about res related with fusion protein
        j = 1
        for lz in lz_files:
            print(f"> {j *100 / len(lz_files)} %")
            if not os.path.exists(f"{path_dimer}{lz}.pdb"):
                print(f"> Removing fusion protein from {lz}")
                rem_lz(
                    lz, lz_info, path_dimer_lz, path_dimer, cutoff_lz
                )  # Open pdb and remove lisozim
                cmd.reinitialize()
            j += 1
            sleep(0.5)
    # quit_pymol()
    #############################################################################
    # MINIMUM DISTANCE (FILTER II)
    #############################################################################
    nlz_files = get_files(path_dimer, ".pdb")
    # nlz_files = ['5F8U_A_ADRB1_MELGA_5F8UA_symB00000100']
    if option_dist == True:
        j = 1
        total = len(nlz_files)
        for nlz in nlz_files:
            print(f"> {j *100 / total} %")
            print(f"> {nlz}")
            try:
                u = mda.Universe(f"{path_dimer}{nlz}.pdb")
            except:
                print(colored("Error open!", "red"))
                os.remove(f"{path_dimer}{nlz}.pdb")
                continue
            name = nlz.split("_")
            g1 = u.select_atoms(f"name CA and segid {name[6][0]}")
            g2 = u.select_atoms(
                f"name CA and segid {name[6][1]} and around {min_dist} protein and name CA and segid {name[1]}"
            )
            dist = mda.analysis.distances.distance_array(g1.positions, g2.positions)
            if dist != []:
                mindist = np.amin(dist)
            else:
                if os.path.exists(f"{path_dimer}{nlz}.pdb"):
                    print(f"> Removing {nlz}")
                    os.remove(f"{path_dimer}{nlz}.pdb")
            j += 1

    ###############################################################################
    # SUPERPOSITION WITH TOPOLOGY (FILTER III)
    ###############################################################################
    nlz_files = get_files(path_dimer, ".pdb")
    # ORIENTATION RIGHT
    # start_pymol()
    # cmd.do('set retain_order, 0') #This setting affects the command super.
    if option_super == True:
        for nlz in nlz_files:
            opm = nlz.split("_")
            opm = opm[0].lower()
            data = requests.get(
                "https://opm-assets.storage.googleapis.com/pdb/" + opm
            ).text
            soup = BeautifulSoup(data, "html.parser")  # Obtain a soup of info.
            if not os.path.exists(f"{path_opm}{opm}.pdb"):
                get_opm(path_opm, opm)
                print(f"> Loading dimer {nlz} for orientation")
                load_pymol(path_dimer, nlz, ".pdb")
                print(f"> Loading opm structure {opm}")
                load_pymol(path_opm, opm, ".pdb")
                super_rmsd = cmd.super(nlz, opm)
                sleep(0.5)
                align_rmsd = cmd.align(nlz, opm)
                sleep(0.5)
                if super_rmsd[0] > align_rmsd[0]:
                    print("> Chains aligned!")
                else:
                    super_rmsd = cmd.super(nlz, opm)
                    sleep(0.5)
                    print("> Chains superposed!")
                # cmd.do('set retain_order, 1')
                cmd.do(f"save {path_dimer}{nlz}.pdb, {nlz}")
                # cmd.do('set retain_order, 0')
                cmd.reinitialize()

        ##GET DICCIONARY OF RESOLUTIONS OF EACH STRUCTURE (CODI PDB)
        try:
            engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')
            print(engine)
            # engine = create_engine(f"mysql+mysqlconnector://adrian:D1m3rB0w!@localhost:3306/lmcdb")
        except:
            print(colored("could not connect", "red"))
            sys.exit(0)

        sql_query = """SELECT pdb, resolution FROM gpcr_pdb
        WHERE resolution < 90 and method='X-ray' and type ='receptor';"""
        data = pd.read_sql(sql_query, engine)
        l_pdbid = list(data.pdb)
        l_res = list(data.resolution)
        d_res = {}
        for i, pdb in enumerate(l_pdbid):
            d_res[pdb] = l_res[i]

        # ##RMSD ANALYZE
        # sleep(1)
        # remove = open(path_misc + "removed.txt", "w")
        # remove.writelines(f"Reference Removed RMSD \n")
        d_ref = {}
        for a, b in itertools.combinations(nlz_files, 2):
            a_pdb = a.split("_")
            a_pdb = a_pdb[0]
            b_pdb = b.split("_")
            b_pdb = b_pdb[0]
            if (
                a[7:12] == b[7:12]
                and os.path.exists(f"{path_dimer}{a}.pdb")
                and os.path.exists(f"{path_dimer}{b}.pdb")
            ):
                load_pymol(path_dimer, a, ".pdb")
                load_pymol(path_dimer, b, ".pdb")
                cmd.do("set retain_order, 0")
                chains_a = cmd.get_chains(a)
                chains_b = cmd.get_chains(b)
                cmd.do("sort")
                # rmsd = cmd.do(f'super {a}, {b}')
                rmsd_a = cmd.align(f"{b} and not het", f"{a} and not het")
                rmsd_s = cmd.super(b, a)
                if round(rmsd_a[0]) < round(rmsd_s[0]):
                    rmsd = rmsd_a
                else:
                    rmsd = rmsd_s
                # print (rmsd)
                if (
                    round(rmsd[0]) <= cutoff_rmsd
                ):  # SELECTED LOOKING DIFERENTS STRUCTURES
                    if (
                        d_res[a[0:4]] < d_res[b[0:4]] or d_res[a[0:4]] == d_res[b[0:4]]
                    ):  # REMOVE THE DIMER WITH MAJOR RESOLUTION
                        print(colored(f"> Removing {b}", "red"))
                        # remove.writelines(f"{a} {b} {rmsd} \n")
                        reference_dict(a_pdb, b_pdb, d_ref)
                        os.remove(f"{path_dimer}{b}.pdb")
                        nlz_files.remove(b)
                    else:
                        print(colored(f"> Removing {a}", "red"))
                        # remove.writelines(f"{a} {b} {rmsd} \n")
                        reference_dict(b_pdb, a_pdb, d_ref)
                        os.remove(f"{path_dimer}{a}.pdb")
                        nlz_files.remove(a)
                else:  # NOW DETECT IF WE SORT THE CHAINS TO INVERSE IT IS THE SAME DIMER
                    # cmd.do('set retain_order, 1')
                    cmd.do(f'alter {b} and chain {chains_b[0]}, chain = "temp"')
                    cmd.do(
                        f'alter {b} and chain {chains_b[1]}, chain = "{chains_b[0]}"'
                    )
                    cmd.do(f'alter {b} and chain temp, chain = "{chains_b[1]}"')
                    cmd.do("sort")
                    rmsd_a = cmd.align(f"{b} and not het", f"{a} and not het")
                    rmsd_s = cmd.super(b, a)
                    if round(rmsd_a[0]) < round(rmsd_s[0]):
                        rmsd = rmsd_a
                    else:
                        rmsd = rmsd_s
                    # print (rmsd)
                    if (
                        round(rmsd[0]) <= cutoff_rmsd
                    ):  # IF IT IS THE SAME REMOVE THE DIMER
                        if (
                            d_res[a[0:4]] < d_res[b[0:4]]
                            or d_res[a[0:4]] == d_res[b[0:4]]
                        ):  # REMOVE THE DIMER WITH MAJOR RESOLUTION
                            print(colored(f"> Removing {b}", "red"))
                            reference_dict(a_pdb, b_pdb, d_ref)
                            os.remove(f"{path_dimer}{b}.pdb")
                            # remove.writelines(f"{a} {b} {rmsd} \n")
                            nlz_files.remove(b)
                        else:
                            print(colored(f"> Removing {a}", "red"))
                            reference_dict(b_pdb, a_pdb, d_ref)
                            os.remove(f"{path_dimer}{a}.pdb")
                            # remove.writelines(f"{a} {b} {rmsd} \n")
                            nlz_files.remove(a)
            else:
                continue
            cmd.reinitialize()
        # REFERENCE PDB CODES TO REMOVED ONES
        reference = open(path_misc + "references.txt", "w")
        reference.writelines(f"Reference PDB_codes \n")
        for key in d_ref.keys():
            reference.writelines(f"{key} {d_ref[key]}\n")
        reference.close()
    # remove.close()

    ###############################################################################
    # RENAME PDBS
    ###############################################################################
    if option_check == True:
        last_id = int(last_id)
        for code in l_pdb_code:
            dimers_id = get_files(path_dimer_end + code + "*", ".pdb")
            dimers = get_files(path_dimer + code + "*", ".pdb")
            i = 0
            possible = False
            possible_dimer = ""
            if dimers_id == [] and not dimers == []:
                for d in dimers:
                    last_id = last_id + 1
                    name = d.split("_")
                    name = name[0] + "_" + name[2]
                    print(colored((code, d, last_id, name), "green"))
                    os.system(
                        f"cp {path_dimer}{d}.pdb {path_dimer_end}{name}_{str(last_id)}.pdb"
                    )
            if dimers == [] and not dimers_id == []:
                for d in dimers_id:
                    name = d.split("_")
                    name = name[0] + "_" + name[2]
                    print(colored((code, d, last_id, name), "red"))
                    os.system(f"mv {path_dimer_end}{d}.pdb {path_removed}{d}.pdb")
                continue
            for d_id in dimers_id:
                load_pymol(path_dimer_end, d_id, ".pdb")
                for i, d in enumerate(dimers):
                    load_pymol(path_dimer, d, ".pdb")
                    rmsd_s = cmd.super(d, d_id)
                    if int(round(rmsd_s[0])) == 0:
                        print(f"> Dimer {d} exists!")
                        continue
                    else:
                        print(f"> Dimer {d} NOT exists!")
                        if possible == True and possible_dimer == d:
                            last_id = last_id + 1
                            name = d.split("_")
                            name = name[0] + "_" + name[2]
                            print(colored((code, d, last_id, name), "green"))
                            os.system(
                                f"cp {path_dimer}{d}.pdb {path_dimer_end}{name}_{str(last_id)}.pdb"
                            )
                        possible = True
                        possible_dimer = d
                    i += 1
                cmd.reinitialize()
    quit_pymol()
##############################################################################
# END

