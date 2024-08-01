#!/usr/bin/env python2

# MAIN IMPORTS
from __future__ import print_function
import pymol
from pymol import cmd
import numpy as np
import sys
import requests
import os
import pandas as pd
import mysql.connector
from sqlalchemy import create_engine
from termcolor import colored
import time 

# VARIABLES
from Tools.SECRETS import *


##############################################################################
# Generates possible GPCR oligomers based on possible symmetry mates

##############################################################################

def get_pdb_opm(pdb, path_opm):
    """Downloads a PDB from OPM"""
    pdb_lower = pdb.lower()
    url = "https://opm-assets.storage.googleapis.com/pdb/" + pdb_lower + ".pdb"
    # print(url)
    req = requests.get(url)
    data = req.text
    file_opmpdb = pdb + "_opm.pdb"
    with open(path_opm + file_opmpdb, "w") as opmpdb:
        opmpdb.write(data)


def cut_pdb_opm(opm_pdb, cut_pdb, sel_chain="A"):
    """Reads a PDB from the OPM an writes a PDB with the TM region"""
    with open(opm_pdb) as input_opm:
        data = input_opm.read()
    if "No such object" in data:
        print("Not in OPM!!!")
        return "no TM"
    ndx_top = data.find("O   DUM")
    ndx_bot = data.find("N   DUM")
    z_top = float(data[ndx_top + 34 : ndx_top + 41])
    z_bot = float(data[ndx_bot + 34 : ndx_bot + 41])
    # take only atoms within the membrane
    lines = data.splitlines()
    out_lines = []
    bot_list, top_list = [], []
    chains = []
    for line in lines:
        ch = line[21:22].strip()
        if not ch in chains and line.startswith("ATOM"):
            chains.append(ch)
    # print(chains)
    # OPM choses a chain which does not need to be A
    if len(chains) == 1 or not sel_chain in chains:
        chain = chains[0]  # use the chain present in OPM
    else:
        chain = sel_chain  # use the selected chain

    for line in lines:
        if line.startswith("ATOM"):
            xcoor = float(line[31:38])
            ycoor = float(line[39:46])
            zcoor = float(line[47:54])
            ch = line[21:22]
            if z_bot < zcoor < z_top and ch == chain:
                out_lines.append(line + "\n")
            if z_bot < zcoor < z_bot + 2 and ch == chain:
                bot_list.append([xcoor, ycoor, zcoor])
            if z_top - 2 < zcoor < z_top and ch == chain:
                top_list.append([xcoor, ycoor, zcoor])
    # calculate center for top and bottom
    bot_list = np.array(bot_list)
    top_list = np.array(top_list)
    # print(top_list)
    if len(bot_list) == 0 or len(top_list) == 0:
        print("Empty list!!!")
        return "no TM"
    bot_point = [np.mean(bot_list[:, i]) for i in range(3)]
    top_point = [np.mean(top_list[:, i]) for i in range(3)]
    # Write output
    with open(cut_pdb, "w") as cut_opm:
        cut_opm.writelines(out_lines)
        cut_opm.write(
            "HETATM {0:4}  O   DUM  {0:4}     {1:7.3f} {2:7.3f} {3:7.3f}\n".format(
                998, top_point[0], top_point[1], top_point[2]
            )
        )
        cut_opm.write(
            "HETATM {0:4}  N   DUM  {0:4}     {1:7.3f} {2:7.3f} {3:7.3f}\n".format(
                999, bot_point[0], bot_point[1], bot_point[2]
            )
        )


def generate_oligomers(pdb, sel_chain, path_cut, uni=""):
    """Generates symmetry mates for a PDB and a selected chain.
    When there are valid oligomers save a PSE session"""
    path_pse = "./data/PSE/"
    # load pdb
    cmd.fetch(pdb)
    time.sleep(1)
    chains = []
    for ch in cmd.get_chains(pdb):
        # print(pdb, " has chain ", ch)
        cmd.copy(pdb + ch, pdb)
        chains.append(ch)
    # print('Chains in PDB', chains)

    cmd.delete(pdb)
    try:
        cmd.load(path_cut + "{0}_cut.pdb".format(pdb))
        time.sleep(1)
    except:
        print("Could not load associated PDB with TM domain only")
        pymol.cmd.quit()
        sys.exit()
    cmd.show_as("cartoon", "all")
    cmd.show("spheres", "resn DUM")

    for ch in chains:
        tmp_chains = chains[:]
        tmp_chains.remove(ch)
        if len(chains) > 1:
            chain_str = "+".join(tmp_chains)  # we create a str for the selection
            remove_str = pdb + ch + " and chain " + chain_str
            cmd.remove(remove_str)

    for ch in chains:
        cmd.super("{0}_cut".format(pdb), pdb + sel_chain)
        time.sleep(1)
        cmd.symexp(
            "sym" + ch, pdb + ch, "({0})".format(pdb + ch), "50"
        )  # cmd.symexp('sym', '4dkl', '(4dkl)', 10)
        time.sleep(1)
    cmd.copy("mv_cut", "{0}_cut".format(pdb))
    
    try:
        o_ref = np.array(cmd.get_atom_coords("resn DUM and name O and {0}_cut".format(pdb)))
    except:
        print("{0}_cut".format(pdb) + " object corrupted!")
        cmd.copy("{0}_cut".format(pdb), "{0}_cut".format(pdb))
        time.sleep(1)
        o_ref = np.array(cmd.get_atom_coords("resn DUM and name O and {0}_cut".format(pdb)))
    n_ref = np.array(cmd.get_atom_coords("resn DUM and name N and {0}_cut".format(pdb)))
    sym_objects = cmd.get_object_list("(sym*)")
    for ch in chains:  # chains[1:]:
        sym_objects.append(pdb + ch)
    sym_objects.sort()
    # print('List of objects', sym_objects)

    v_ref = (o_ref - n_ref) / np.linalg.norm(o_ref - n_ref)

    for sym_object in sym_objects:
        try:
            cmd.super("mv_cut", sym_object)
            time.sleep(1)
        except: # empty selections
            print("removing " + sym_object)
            cmd.delete(sym_object)
            continue
        try: # mv_cut corrupted
            o_mv = np.array(cmd.get_atom_coords("resn DUM and name O and mv_cut"))
        except:
            print("mv_cut object corrupted!")
            cmd.copy("mv_cut", "mv_cut")
            time.sleep(1)
            o_mv = np.array(cmd.get_atom_coords("resn DUM and name O and mv_cut"))
        n_mv = np.array(cmd.get_atom_coords("resn DUM and name N and mv_cut"))
        oo = np.linalg.norm(o_mv - o_ref)
        nn = np.linalg.norm(n_mv - n_ref)
        no = np.linalg.norm(n_mv - o_ref)
        on = np.linalg.norm(o_mv - n_ref)
        v_mv = (o_mv - n_mv) / np.linalg.norm(o_mv - n_mv)
        # http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
        # angle between normal vectors
        angle = np.arccos(np.clip(np.dot(v_ref, v_mv), -1.0, 1.0))
        # print(oo, nn, no, on)
        # different layer
        if abs(no - on) > 10:
            print(sym_object, 'in another layer')
            cmd.delete(sym_object)
        # same layer but far (48 A maximum)
        elif nn > max_dist and oo > max_dist:
            print(sym_object, 'in the same layer but too far')
            print(sym_object, nn, oo)
            cmd.delete(sym_object)
        # not parallel
        elif angle > np.deg2rad(max_angle):  # angle smaller than 30 degrees only
            print(sym_object, 'not parallel to reference', angle)
            cmd.delete(sym_object)
        else:
            print(sym_object, "OK", "angle: ", angle, "dists: ", nn, oo)
            pass
    cmd.delete("{0}_cut".format(pdb))
    cmd.delete("mv_cut")
    try:
        cmd.color("cyan", pdb + sel_chain)
        cmd.show("sticks", "hetatm")
        cmd.show("spheres", "hetatm")
        cmd.set("sphere_scale", "0.2")
        # Save PyMOL session only if there are valid interfaces
        final_objects = cmd.get_object_list("*")
        # print(len(final_objects))
        if len(final_objects) > 1:
            ses_name = pdb + "_" + sel_chain + "_" + uni + ".pse"
            print("Saving a session:", ses_name, "\n")
            cmd.save(path_pse + ses_name, "*")
        else:
            print("No oligos here\n")
    except: # ref chain structures removed only symm structures
        print("No oligos here\n")
        pass
    cmd.delete("*")

def all_in_one(pdb, sel_chain, uni=""):
    """1) Requests the PDB to OPM with get_pdb_opm
       2) Cuts the TM parts with cut_pdb_opm
       3) Generates all dimers and creates .pse sessions with generate_oligomers
       4) Removes PDBs (intermediate)"""
    path_pse = "./data/PSE/"
    path_opm = "./data/OPM/"
    path_cut = "./data/OPM/cut_opm/"
    ses_name = pdb + "_" + sel_chain + "_" + uni + ".pse"
    if not os.path.exists(f"{path_pse}{ses_name}"):
        pdb_lower = pdb.lower()
        if not os.path.exists(f"{path_opm}{pdb}_opm.pdb"):
            get_pdb_opm(pdb, path_opm)  # output is xxx_opm.pdb
        if not os.path.exists(f"{path_cut}{pdb}_cut.pdb"):
            message = cut_pdb_opm(
                path_opm + pdb + "_opm.pdb", path_cut + pdb + "_cut.pdb", sel_chain
            )
        else:
            message = "yes TM"
        if not message == "no TM":
            generate_oligomers(pdb, sel_chain, path_cut, uni)  # RUN PYMOL
            # removes unnecessary files
            fpdb = pdb.lower() + ".pdb"
            if os.path.isfile(fpdb):  # when PDB comes in CIF format
                os.remove(fpdb)
            else:
                fpdb = pdb.lower() + ".cif"
                if os.path.isfile(fpdb):
                    os.remove(fpdb)
            os.remove(path_cut + pdb + "_cut.pdb")
        else:
            print("no TMs, skipped")
        # os.remove(path_opm + pdb + "_opm.pdb")
    else:
        print(f"This file {ses_name} exists!")


def pdb_lmcdb(last_pdb):
    """Returns the list of PDB codes from lmcdb"""

    try:
        engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')

        print(engine)
    except:
        print("could not connect")
        sys.exit(0)

    ##########################################################################
    # USE FOR TESTING
    # sql_query = """SELECT pdbid, chain, uniprot FROM pdb
    # WHERE resolution < 90 and domain like '7tm%' or 'Frizzled' and pdbid >= '5UNH' and
    # expmet='X-RAY DIFFRACTION';"""
    if last_pdb != "NO":
        sql_query = f"SELECT pdb, chain, uniprot_entry FROM gpcr_pdb WHERE resolution < 90 and method='X-ray' and type ='receptor' and pdb >= '{last_pdb}';"
    else:
        sql_query = f"SELECT pdb, Chain, uniprot_entry FROM gpcr_pdb WHERE resolution < 90 and method='X-ray' and type ='receptor';"

    data = pd.read_sql(sql_query, engine)
    pdb_codes = list(data.pdb)
    pdb_chains = list(data.chain)
    pdb_uni = list(data.uniprot_entry)

    return pdb_codes, pdb_chains, pdb_uni


##############################################################################
# MAIN PROGRAM STARTS HERE
max_dist = 48  #  dist in angstroms between monomers
max_angle = 30  # in degrees between normal vectors

# get list of PDBs from the LMCDB
pdbs, chains, unis = pdb_lmcdb(sys.argv[1])

# run PyMOL stuff
pymol.finish_launching(["pymol", "-q"])

for i, u_pdb in enumerate(pdbs):
    pdb = str(u_pdb)
    chns = str(chains[i]).split(", ")
    uni = str(unis[i]).lower()
    for sel_chain in chns:
        print("###", pdb, sel_chain, "###")
        # try:
        all_in_one(pdb, sel_chain, uni)
        # except:
        #     print(colored("### ERROR " + pdb + sel_chain + "###", "red"))
        #     continue
pymol.cmd.quit()
