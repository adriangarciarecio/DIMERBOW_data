#!/usr/bin/python3

###############################################################################
# IMPORTS
###############################################################################
import os
import sys
import shutil
import glob
from gen_distribution import *
path =  os.path.dirname(os.path.abspath('__file__'))
sys.path.insert(0, path)
import math
import requests
import MDAnalysis as mda
from itertools import combinations
from bs4 import BeautifulSoup
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pymol
from pymol import cmd
from operator import itemgetter
import numpy as np
import matplotlib.mlab as mlab
import urllib
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time

###############################################################################
# FUNCTIONS
###############################################################################
def get_all_pdb(f_r, wat_ions):
    '''Obtain protein, ligand and wat_ions in separated lists from a pdb file'''
    letters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    protein, ligand, chains, lipids, wat_ion, rest = [], [], [], [], [], []
    for line in f_r:
        if line.startswith("ATOM") or line.startswith("TER"):
            protein.append(line)
        if line.startswith("HETATM"):
            molec = line[17:20].strip() #All pdb have the same structure and between this lines we can found the ligand name.
            if molec in wat_ions:
                wat_ion.append(line)
            else:
                ligand.append(line)
                if molec not in lipids and molec not in wat_ions: #Gets a list with all name of lipids in pdb file
                    lipids.append(molec)
        if not (line.startswith("CONECT") or line.startswith("END") or line.startswith("NUMMDL") or line.startswith("MODEL")):
            l = line.split(" ")
            l = list(filter(None, l))
            if l[4][0] not in chains and l[4][0] in letters:
                chains.append(l[4][0])
        if not line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
            rest.append(line)
    return protein, ligand, chains, lipids, rest

def get_all_gro(f_r, wat_ions):
    aa = ['TRP', 'TYR','PHE', 'ILE', 'LEU', 'VAL', 'MET', 'ALA', 'PRO', 'CYS','GLY', 'SER','THR', 'ASN', 'GLN', 'HIS', 'ARG','LYS','GLU','ASP']
    protein, popc, w, ion, rest = [] , [] ,[] ,[] ,[]
    for g in f_r:
        g_split = g.split (' ')
        g_split = list(filter(None, g_split))
        if g_split[0][-3:] in aa:
            protein.append(g_split)
        if "POPC" in g_split[0]:
            popc.append(g_split)
        if "W" in g_split[0]:
            w.append(g_split)
        if g_split[0][-3:] in wat_ions:
            ion.append(g_split)
        if not g_split[0][-3:] in aa and not "POPC" in g_split[0] and not "W" in g_split[0] and not g_split[0][-3:] in wat_ions:
            rest.append(g_split)
    return protein, popc, w, ion, rest

def rename_residues(gro_r, path, grofile, ref_r):
    aa = ['TRP', 'TYR','PHE', 'ILE', 'LEU', 'VAL', 'MET', 'ALA', 'PRO', 'CYS','GLY', 'SER','THR', 'ASN', 'GLN', 'HIS', 'ARG','LYS','GLU','ASP']
    f_out = open(f'{path}2_{grofile}', "w")
    f_out.writelines(gro_r[0:2])
    for i, g in enumerate(gro_r):
        g_split = g.split (' ')
        ref_split = ref_r[i].split(' ')
        for j, element in enumerate(g_split):
            if element[-3:] in aa and ref_split[j][-3:] in aa:
                g_split[j] = ref_split[j]
                line = " ".join(g_split)
                f_out.writelines(line)
                break
            if element[-3:] in aa and ref_split[j-1][-3:] in aa:
                g_split[j-1] = ref_split[j-1]
                del g_split[j]
                line = " ".join(g_split)
                f_out.writelines(line)
                break
            if element[-3:] in aa and ref_split[j-2][-3:] in aa:
                g_split[j-2] = ref_split[j-2]
                del g_split[j]
                del g_split[j-1]
                line = " ".join(g_split)
                f_out.writelines(line)
                break
            if element[-3:] in aa and ref_split[j-3][-3:] in aa:
                g_split[j-3] = ref_split[j-3]
                del g_split[j]
                del g_split[j-1]
                del g_split[j-2]
                line = " ".join(g_split)
                f_out.writelines(line)
                break
            if 'POPC' in element:
                f_out.writelines(g)
                break
            if element.endswith('W'):
                f_out.writelines(g)
                break
            if 'ION' in element:
                f_out.writelines(g)
                break
    f_out.writelines(gro_r[-1])

def get_opm_coords(pdb):
    '''Get the coordinates to determine the differents zones on the structure of a protein.'''
    opm_coords = {}#Creates a dictionary for call each coordinates, extracted from opm, according to helix number
    driver.get("https://opm.phar.umich.edu/protein.php?search=" + pdb)
    time.sleep(2)
    all_info = driver.find_element_by_class_name("protein-content")
    table = all_info.find_element_by_class_name("small-text.break")
    coordinates = table.text
    coordinates = coordinates.split(":")
    coords = coordinates[2]
    coords = coords.split(",")
    for i, seg in enumerate(coords):
        c = seg.partition('(')[-1].rpartition(')')[0]
        c = c.replace(" ","").split("-")
        opm_coords[i+1] = c
    return opm_coords

def save_interac_pdb (u, prot, sear, path_out, namefile, coords):
    '''Creates a file that contain all interactions on specfics parameters.'''
    sel_prot = u.select_atoms(prot)
    near_atoms = u.select_atoms(sear)
    f_out = open(path_out + "ints_" + namefile + ".txt", "a")
    inters_all = get_interactions_pdb(sel_prot, near_atoms, namefile, coords)
    inters_all = refine_int_list(inters_all) #Gets the less distance found it between two atoms
    for inter in inters_all:
        f_out.write(inter + "\n")
    f_out.close()

def save_interac_gro (u, prot, sear, path_out, namefile, coords):
    '''Creates a file that contain all interactions on specfics parameters.'''
    sel_prot = prot
    near_atoms = sear
    f_out = open(path_out + "ints_" + namefile + ".txt", "a")
    inters_all = get_interactions_gro(sel_prot, near_atoms, namefile, coords)
    inters_all = refine_int_list(inters_all) #Gets the less distance found it between two atoms
    for inter in inters_all:
        f_out.write(inter + "\n")
    f_out.close()

def refine_int_list(ini_int_list):
    '''Filter the results obtained on get_interactions. A list without repeticions,
    and getting the lowest distance'''
    new_L = []
    temp_L = []
    final_L = []
    for i in ini_int_list:
        inter = i.split()
        prot = inter[0]
        lip = inter[1]
        dist = inter[2]
        r2 = lip.split("_")
        r1 = prot.split("_")
        lip_num = r2[0]+r2[3] #A616 o B616
        resid = r1[0]+ r1[2] + r1[3] #AALA68
        x = resid, lip_num #AALA68, A616
        y = resid, lip_num, dist #AALA68, A616, 9.73
        if x not in new_L:
            new_L.append(x)
            temp_L.append(y)
            final_L.append(i)
        else:
            for k in temp_L:
                if k[0] == y[0] and k[1] == y[1]: #In case we have repetitive interactions
                    if float(dist) < float(k[2]): #Gets the interaction with the shortest distance
                        indice = temp_L.index(k)
                        temp_L.remove(k)
                        del final_L[indice]
                        temp_L.append(y)
                        final_L.append(i)
    return final_L

def get_interactions_pdb(protein_a, protein_b, namefile, coords):
    '''Obtain all interactions between prot a and prot b.
    FORMAT: cha_ident_resprot_pos    cha_ident_reslip_pos_dist --> A_CB_PRO_224 A_C13_1PE_617 8.78'''
    protein_a_posits = protein_a.positions
    protein_b_posits = protein_b.positions
    interactions = []
    for i, prot_a_atom_coord in enumerate(protein_a_posits):
        for j, prot_b_atom_coord in enumerate(protein_b_posits):
            dist = Euc_dist(prot_a_atom_coord, prot_b_atom_coord)
            if dist <= 5:
                coord_z = get_height(prot_a_atom_coord, prot_b_atom_coord)
                if namefile == '5JQH_ADRB2_44':
                    segment_a = get_zone(protein_a[i].resid-1000, coords) #Returns the zone where is located the atom a
                    segment_b = get_zone(protein_b[j].resid-1000, coords)
                else:
                    segment_a = get_zone(protein_a[i].resid, coords) #Returns the zone where is located the atom a
                    segment_b = get_zone(protein_b[j].resid, coords)
                typea, typeb = get_class_pdb(protein_a[i].resname, protein_a[i].name, protein_b[j].resname, protein_b[j].name)
                info_prot_a_atom = protein_a[i].segid + "_" + protein_a[i].name + "_" + protein_a[i].resname + "_" + str(protein_a[i].resid) + "_" + segment_a
                info_prot_b_atom = protein_b[j].segid + "_" + protein_b[j].name + "_" + protein_b[j].resname + "_" + str(protein_b[j].resid) + "_" + segment_b
                info_inter = info_prot_a_atom + " " +  info_prot_b_atom + " " +  ('%.2f' % dist) + " " + coord_z + " " + typea + typeb
                interactions.append(info_inter)
    return interactions

def get_interactions_gro(protein_a, protein_b, namefile, coords):
    '''Obtain all interactions between prot a and prot b.
    FORMAT: cha_ident_resprot_pos    cha_ident_reslip_pos_dist --> A_CB_PRO_224 A_C13_1PE_617 8.78'''
    protein_a_posits = protein_a.positions
    protein_b_posits = protein_b.positions
    interactions = []
    for i, prot_a_atom_coord in enumerate(protein_a_posits):
        for j, prot_b_atom_coord in enumerate(protein_b_posits):
            dist = Euc_dist(prot_a_atom_coord, prot_b_atom_coord)
            if dist <= 5:
                coord_z = get_height(prot_a_atom_coord, prot_b_atom_coord)
                if '5JQH_ADRB2_44' in namefile:
                    segment_a = get_zone(protein_a[i].resid-1000, coords) #Returns the zone where is located the atom a
                    segment_b = get_zone(protein_b[j].resid-1000, coords)
                else:
                    segment_a = get_zone(protein_a[i].resid, coords) #Returns the zone where is located the atom a
                    segment_b = get_zone(protein_b[j].resid, coords)
                typea, typeb = get_class_gro(protein_a[i].resname, protein_a[i].name, protein_b[j].resname, protein_b[j].name)
                info_prot_a_atom = protein_a[i].name + "_" + protein_a[i].resname + "_" + str(protein_a[i].resid) + "_" + segment_a
                info_prot_b_atom = protein_b[j].name + "_" + protein_b[j].resname + "_" + str(protein_b[j].resid) + "_" + segment_b
                info_inter = info_prot_a_atom + " " +  info_prot_b_atom + " " +  ('%.2f' % dist) + " " + coord_z + " " + typea + typeb
                interactions.append(info_inter)
    return interactions

def Euc_dist(p, q):
    '''Calculates the distance between two coordinates'''
    dist = math.sqrt((q[0]-p[0])**2+(q[1]-p[1])**2+(q[2]-p[2])**2)
    return dist

def get_height(coords_a, coords_b):
    '''Return the average between two coordinates z.'''
    coord_za = coords_a[2:3]
    coord_zb = coords_b[2:3]
    coord_z = str((coord_za + coord_zb)/2) #Average
    if coord_z == "[ 0.]": #Simplify that value
        coord_z = "0"
    else:
        coord_z = coord_z.partition('[')[-1].rpartition(']')[0] #Take a part of a str between specific symbols.
    return coord_z

def get_zone(prot_res, coords):
    '''Get the zone where is localized the residue on GPCR.'''
    dic_coords = coords
    coords = []
    for coord in dic_coords:
        c = dic_coords[coord]
        coords.append([c[0],c[1]])
    if prot_res < int(coords[0][0]):
        zone = 'NT'
    if prot_res >= int(coords[0][0]) and prot_res <= int(coords[0][1]):
        zone = 'TM1'
    if prot_res > int(coords[0][1]) and prot_res < int(coords[1][0]):
        zone = 'ICL1'
    if prot_res >= int(coords[1][0]) and prot_res <= int(coords[1][1]):
        zone = 'TM2'
    if prot_res > int(coords[1][1]) and prot_res < int(coords[2][0]):
        zone = 'ECL1'
    if prot_res >= int(coords[2][0]) and prot_res <= int(coords[2][1]):
        zone = 'TM3'
    if prot_res > int(coords[2][1]) and prot_res < int(coords[3][0]):
        zone = 'ICL2'
    if prot_res >= int(coords[3][0]) and prot_res <= int(coords[3][1]):
        zone = 'TM4'
    if prot_res > int(coords[3][1]) and prot_res < int(coords[4][0]):
        zone = 'ECL2'
    if prot_res >= int(coords[4][0]) and prot_res <= int(coords[4][1]):
        zone = 'TM5'
    if prot_res > int(coords[4][1]) and prot_res < int(coords[5][0]):
        zone = 'ICL3'
    if prot_res >= int(coords[5][0]) and prot_res <= int(coords[5][1]):
        zone = 'TM6'
    if prot_res > int(coords[5][1]) and prot_res < int(coords[6][0]):
        zone = 'ECL3'
    if prot_res >= int(coords[6][0]) and prot_res <= int(coords[6][1]):
        zone = 'TM7'
    if prot_res > int(coords[6][1]) and len(coords) == 7 and not prot_res >= 1000:
        zone = 'CT'
    if len(coords) == 8 and prot_res > int(coords[6][1]) and prot_res < int(coords[8][0]):
        zone = 'ICL4'
    if len(coords) == 8 and prot_res >= int(coords[7][0]) and prot_res <= int(coords[8][1]):
        zone = 'TM8'
    if len(coords) == 8 and prot_res > int(coords[7][1]) and not prot_res >= 1000:
        zone = 'CT'
    if prot_res >= 1000:
        zone = 'LZ'
    return zone

def get_class_pdb(resnamea, atomnamea, resnameb, atomnameb):
    aa = ['TRP', 'TYR','PHE', 'ILE', 'LEU', 'VAL', 'MET', 'ALA', 'PRO', 'CYS','GLY', 'SER','THR', 'ASN', 'GLN', 'HIS', 'ARG','LYS','GLU','ASP']
    hidrofobic = ['GLY', 'PRO', 'ALA', 'MET', 'VAL', 'LEU', 'ILE']
    d_gly = {'O':'P', 'N':'P', 'C':'H', 'CA':'H'}
    d_ala = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H'}
    d_val = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG1':'H', 'CG2':'H'}
    d_pro = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'CD':'H'}
    d_leu = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'CD1':'H', 'CD2':'H'}
    d_ile = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG1':'H', 'CG2':'H', 'CD1':'H'}
    d_met = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'CE':'H', 'SD':'H'}

    arom = ['PHE', 'TYR', 'TRP']
    d_phe = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'A', 'CD1':'A', 'CD2':'A', 'CE1':'A', 'CE2':'A', 'CZ':'A'}
    d_tyr = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'A', 'CD1':'A', 'CD2':'A', 'CE1':'A', 'CE2':'A', 'CZ':'A', 'OH':'P'}
    d_trp = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'A', 'CD1':'A', 'NE1':'P', 'CD2':'A', 'CE2':'A', 'CE3':'A', 'CZ3':'A', 'CH2':'A', 'CZ2':'A'}

    charged = ['ASP', 'GLU','LYS','ARG']
    d_arg = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'CD':'H', 'NE':'C', 'CZ':'H', 'NH1':'C', 'NH2':'C'}
    d_lys = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'CD':'H', 'CE':'H', 'NZ':'C'}
    d_asp = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H','CG':'H', 'OD1':'C', 'OD2':'C'}
    d_glu = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H','CG':'H', 'CD':'H', 'OE1':'C', 'OE2':'C'}

    polar = ['CYS', 'GLN', 'ASN', 'THR', 'SER', 'HIS']
    d_cys = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'SG':'H'}
    d_asn = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'OD1':'P', 'ND2':'P'}
    d_gln = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'CD':'H', 'OE1':'P', 'NE2':'P'}
    d_ser = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'OG':'P'}
    d_thr = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG2':'H', 'OG1':'P'}
    d_his = {'O':'P', 'N':'P', 'C':'H', 'CA':'H', 'CB':'H', 'CG':'H', 'ND1':'P', 'CD2':'H', 'NE2':'P', 'CE1':'H'}

    d_all = {'GLY':d_gly, 'PRO':d_pro, 'ALA':d_ala, 'MET':d_met, 'VAL':d_val, 'LEU':d_leu, 'ILE':d_ile, 'PHE':d_phe, 'TYR':d_tyr, 'TRP':d_trp, 'ASP':d_asp, 'GLU':d_glu,'LYS':d_lys,'ARG':d_arg, 'CYS':d_cys, 'GLN':d_gln, 'ASN':d_asn, 'THR':d_thr, 'SER':d_ser, 'HIS':d_his }

    if resnamea in aa and resnameb in aa:
        d_resa = d_all[resnamea]
        d_resb = d_all[resnameb]
        class_typea = d_resa[atomnamea]
        class_typeb = d_resb[atomnameb]
    if resnamea not in aa or resnameb not in aa:
        class_typea = 'X'
        class_typeb = 'X'
    return class_typea, class_typeb

def get_class_gro(resnamea, atomnamea, resnameb, atomnameb):
    aa = ['TRP', 'TYR','PHE', 'ILE', 'LEU', 'VAL', 'MET', 'ALA', 'PRO', 'CYS','GLY', 'SER','THR', 'ASN', 'GLN', 'HIS', 'ARG','LYS','GLU','ASP']
    hidrofobic = ['GLY', 'PRO', 'ALA', 'MET', 'VAL', 'LEU', 'ILE']
    d_gly = {'BB':'P'}
    d_ala = {'BB':'P'}
    d_val = {'BB':'P', 'SC1':'C'}
    d_pro = {'BB':'P', 'SC1':'C'}
    d_leu = {'BB':'P', 'SC1':'C'}
    d_ile = {'BB':'P', 'SC1':'C'}
    d_met = {'BB':'P', 'SC1':'Q'}

    arom = ['PHE', 'TYR', 'TRP']
    d_phe = {'BB':'P', 'SC1':'C', 'SC2':'C', 'SC3':'C'}
    d_tyr = {'BB':'P', 'SC1':'C', 'SC2':'C', 'SC3':'P'}
    d_trp = {'BB':'P', 'SC1':'C', 'SC2':'P', 'SC3':'C', 'SC4':'C'}

    charged = ['ASP', 'GLU','LYS','ARG']
    d_arg = {'BB':'P', 'SC1':'N', 'SC2':'Q'}
    d_lys = {'BB':'P', 'SC1':'C', 'SC2':'Q'}
    d_asp = {'BB':'P', 'SC1':'Q'}
    d_glu = {'BB':'P', 'SC1':'Q'}

    polar = ['CYS', 'GLN', 'ASN', 'THR', 'SER', 'HIS']
    d_cys = {'BB':'P', 'SC1':'N'}
    d_asn = {'BB':'P', 'SC1':'P'}
    d_gln = {'BB':'P', 'SC1':'P'}
    d_ser = {'BB':'P', 'SC1':'P'}
    d_thr = {'BB':'P', 'SC1':'P'}
    d_his = {'BB':'P', 'SC1':'C', 'SC2':'P', 'SC3':'P'}

    d_all = {'GLY':d_gly, 'PRO':d_pro, 'ALA':d_ala, 'MET':d_met, 'VAL':d_val, 'LEU':d_leu, 'ILE':d_ile, 'PHE':d_phe, 'TYR':d_tyr, 'TRP':d_trp, 'ASP':d_asp, 'GLU':d_glu,'LYS':d_lys,'ARG':d_arg, 'CYS':d_cys, 'GLN':d_gln, 'ASN':d_asn, 'THR':d_thr, 'SER':d_ser, 'HIS':d_his }

    if resnamea in aa and resnameb in aa:
        d_resa = d_all[resnamea]
        d_resb = d_all[resnameb]
        class_typea = d_resa[atomnamea]
        class_typeb = d_resb[atomnameb]
    if resnamea not in aa or resnameb not in aa:
        class_typea = 'X'
        class_typeb = 'X'
    return class_typea, class_typeb

def get_numpy_direct(o_files, path_out, ext_out, select):
    '''Return an array with all pairs of aa that concur.'''
    aa = ['TRP', 'TYR','PHE', 'ILE', 'LEU', 'VAL', 'MET', 'ALA', 'PRO', 'CYS','GLY', 'SER','THR', 'ASN', 'GLN', 'HIS', 'ARG','LYS','GLU','ASP']
    d_aa = {'CYS': 9, 'GLN': 14, 'ILE': 3, 'SER': 11, 'VAL': 5, 'GLY': 10, 'PRO': 8, 'LYS': 17, 'ASP': 19, 'THR': 12, 'PHE': 2, 'ALA': 7, 'MET': 6, 'HIS': 15, 'GLU': 18, 'LEU': 4, 'ARG': 16, 'TRP': 0, 'ASN': 13, 'TYR': 1, 'LIG': 20}
    if select == 'direct':
        a = np.zeros((20,20), dtype = np.float)
    if select == 'ligand':
        a = np.zeros((20,1), dtype = np.float)
    for namefile in o_files:
        text_r = open_file(path_out,namefile + ext_out)
        for r in text_r:
            r = r.split(" ")
            aa1 = r[0].split("_")
            name_aa1 = aa1[2]
            aa2 = r[1].split("_")
            name_aa2 = aa2[2]
            distance = r[2]
            interact = r[-1]
            if name_aa1 in aa and name_aa2 in aa and select =='direct':
                cutoff = get_cutoff(distance,interact)
                if distance <= str(cutoff):
                    value = a.item((d_aa[name_aa1], d_aa[name_aa2])) #Return the value with position x, y (determinated by the position of the aa in d_aa)on the array
                    np.put(a[d_aa[name_aa1]], d_aa[name_aa2], value + 1) #Change the value in the array with determinated coordinates.
                    np.put(a[d_aa[name_aa2]], d_aa[name_aa1], value + 1) #For symmetry
            if name_aa1 not in aa or name_aa2 not in aa and select =='ligand':
                if name_aa2 in aa:
                    value = a.item((d_aa[name_aa2],0))
                    np.put(a[d_aa[name_aa2]],0, value + 1)
                if name_aa1 in aa:
                    value = a.item((d_aa[name_aa1],0))
                    np.put(a[d_aa[name_aa1]],0, value + 1)
    #total = np.sum(a)
    #print ('DIRECT'+ select + ': '+ str(total))
    maxim = np.nanmax(a) #Return the maximum of an array or maximum along an axis, ignoring any NaNs.
    minim = np.nanmin(a)
    a = a-minim
    a = a / maxim #Need normalized the array
    return a

def get_numpy_indirect(o_files, path_out, ext_out):
    '''Return an array with all pair of aa that they interaction is mediated by a ligand. '''
    aa = ['TRP', 'TYR','PHE', 'ILE', 'LEU', 'VAL', 'MET', 'ALA', 'PRO', 'CYS','GLY', 'SER','THR', 'ASN', 'GLN', 'HIS', 'ARG','LYS','GLU','ASP']
    d_aa = {'CYS': 9, 'GLN': 14, 'ILE': 3, 'SER': 11, 'VAL': 5, 'GLY': 10, 'PRO': 8, 'LYS': 17, 'ASP': 19, 'THR': 12, 'PHE': 2, 'ALA': 7, 'MET': 6, 'HIS': 15, 'GLU': 18, 'LEU': 4, 'ARG': 16, 'TRP': 0, 'ASN': 13, 'TYR': 1}
    a = np.zeros((20,20), dtype = np.float)
    for namefile in o_files:
        text_r = open_file(path_out,namefile + ext_out)
        l_1chain, l_2chain, l_chain, l_final = [], [], [], []
        for r in text_r:
            r = r.split(" ")
            aa1 = r[0].split("_")
            chain = aa1[0]
            if chain not in l_chain:
                l_chain.append(chain)
        for i in range (len(l_chain)):
            for r in text_r:
                r = r.split(" ")
                aa1 = r[0].split("_")
                name_aa1 = aa1[2]
                aa2 = r[1].split("_")
                name_aa2 = aa2[2]
                z = r[-1].split("\n")
                z = z[0]
                if ((name_aa1 not in aa) or (name_aa2 not in aa)) and l_chain[i] == aa1[0]:
                    if i == 0:
                        l_1chain.append([aa1[0]+ "_" + aa1[1] + "_" + aa1[2] + "_" + aa1[3], aa2[0]+ "_" + aa2[1] + "_" + aa2[2] + "_" + aa2[3], z])
                    if i == 1:
                        l_2chain.append([aa1[0]+ "_" + aa1[1] + "_" + aa1[2] + "_" + aa1[3], aa2[0]+ "_" + aa2[1] + "_" + aa2[2] + "_" + aa2[3], z])
        for i in range(len(l_1chain)):
            for j in range(len(l_2chain)):
                info1 = l_1chain[i][0].split("_") #A_OD2_ASP_60
                info1lip = l_1chain[i][1].split("_") #A_O4_PEG_1204
                info2 = l_2chain[j][0].split("_") #B_OD2_ASP_60
                info2lip = l_2chain[j][1].split("_")#B_O4_PEG_1204
                if info1lip[0]+info1lip[2]+info1lip[3] == info2lip[0]+info2lip[2]+info2lip[3] and info1[2] in aa and info2[2] in aa:
                    l_final.append([info1[0], info1[1], info1[2], info1[3],l_1chain[i][2], info1lip [0], info1lip[1], info1lip[2], info1lip[3], info2[0], info2[1], info2[2], info2[3], info2lip [0], info2lip[1], info2lip[2], info2lip[3], l_2chain[j][2]])
            i += 1
        for l in l_final:
                name_aa1 = l[2]
                name_aa2 = l[11]
                if name_aa1 in aa and name_aa2 in aa:
                    value = a.item((d_aa[name_aa1], d_aa[name_aa2])) #Return the value with position x, y (determinated by the position of the aa in d_aa)on the array
                    np.put(a[d_aa[name_aa1]], d_aa[name_aa2], value + 1) #Change the value in the array with determinated coordinates.
                    np.put(a[d_aa[name_aa2]], d_aa[name_aa1], value + 1) #For symmetry
    #total = np.sum(a)
    #print ('INDIRECT: '+ str(total))
    maxim = np.nanmax(a) #Return the maximum of an array or maximum along an axis, ignoring any NaNs.
    minim = np.nanmin(a)
    a = a-minim
    a = a / maxim#Need normalized the array
    return a

def get_heatmap(a,title, path_fig):
    '''Return an heatmap from a numpy array.'''
    l_aa = ['W', 'Y', 'F', 'I', 'L', 'V', 'M', 'A', 'P', 'C', 'G', 'S', 'T', 'N', 'Q', 'H', 'R', 'K', 'E', 'D', 'X']
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(a, cmap = plt.cm.Reds, alpha = 0.8) #Creates the heatmap
    if a.shape == (20,20):
        fig.set_size_inches(16,13)
    if a.shape == (20,1):
        fig.set_size_inches(0.5, 13)
        lig = ['LIG']
    ax.set_frame_on(False)
    ax.set_yticks(np.arange(a.shape[0]) + 0.5, minor = False) #Put ticks of labels in the middle of the shells.
    ax.set_xticks(np.arange(a.shape[1]) + 0.5, minor = False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    ax.set_xlabel(title)
    ax.set_ylabel("|     Charged     |  |             Polar         |  |                 Hidrofobic                  |  | Aromatic |")
    if a.shape == (20,20):
        ax.set_xticklabels(l_aa, minor = False)
    if a.shape == (20,1):
        ax.set_xticklabels('X', minor = False) #Write the tick labels on l_aa
    ax.set_yticklabels(l_aa, minor = False)
    ax.grid(False)
    cbar = fig.colorbar(heatmap, ticks = [0, 0.5, 1], orientation='horizontal') #Add a legend to know what meaning have each color.
    cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar
    plt.savefig(path_fig + str(datetime.today()) +'_'+ title +'.png', dpi = 100, bbox_inches='tight')
    plt.close()

def get_hist_dist(c_files, path_class, interact, path_fig):
    data = get_data(c_files, path_class, 'distance', interact)
    mean = np.mean(data)
    std = np.std(data)
    title = 'HISTOGRAM OF DISTANCE ' + interact.upper()
    plt.figure()
    n, bins, patches = plt.hist(data, bins = np.arange(2.2, 5.2, 0.1), facecolor='blue', align='mid', rwidth=1, alpha=0.5) #normerd= 1
    plt.xlabel('DISTANCE')
    plt.ylabel('QUANTITY')
    plt.title(title)
    plt.grid(False)
    plt.xticks(np.arange(2.2, 5.2, 0.2))
    plt.xticks(rotation=70)
    plt.yticks(np.arange(0, 20, 5))
    plt.axis('tight')# Tweak spacing to prevent clipping of ylabel
    plt.savefig(path_fig + str(datetime.today()) +'_'+ title +'.png', dpi = 100)

def get_hist_z(c_files, path_class, path_fig):
    data_hh, hh = get_data(c_files, path_class, 'height', 'HH')
    data_ha, ha = get_data(c_files, path_class, 'height', 'HA')
    data_hp, hp = get_data(c_files, path_class, 'height', 'HP')
    data_hc, hc = get_data(c_files, path_class, 'height', 'HC')
    data_aa, aa = get_data(c_files, path_class, 'height', 'AA')
    data_ap, ap = get_data(c_files, path_class, 'height', 'AP')
    data_ac, ac = get_data(c_files, path_class, 'height', 'AC')
    data_pp, pp = get_data(c_files, path_class, 'height', 'PP')
    data_pc, pc = get_data(c_files, path_class, 'height', 'PC')
    data_cc, cc = get_data(c_files, path_class, 'height', 'CC')
    data = [data_hh, data_ha, data_hp, data_hc, data_aa, data_ap, data_ac, data_pp, data_pc, data_cc]
    title = 'HISTOGRAM OF HEIGHT'
    fig, ax = plt.subplots(figsize =(15,10))
    colors = ['#FF0000', '#FF9E00','#FF00EF','#7E0000', '#EFFF00','#09FF00','#786800','#00FFD5','#00786E','#675F5F']
    interacts = ['HH', 'HA', 'HP', 'HC', 'AA ','AP','AC','PP', 'PC', 'CC']
    ax.hist(data, bins = np.arange(-46, 46,2), histtype = 'bar', stacked = True, color = colors, label = interacts, orientation = 'horizontal', rwidth=1, alpha = 0.75, zorder = 10, edgecolor = 'black')
    ax.set(title = None, ylabel = 'HEIGHT', xlabel = 'NUMBER OF RESIDUE-RESIDUE INTERACTIONS')
    ax.grid(False)
    plt.xticks(np.arange(0,210, 10))
    plt.yticks(np.arange(-45, 45, 5))
    box = ax.get_position()# Shrink current axis's height by 10% on the bottom
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))## Put a legend to the right of the current axis
    plt.axhline(y=15, xmin=0, xmax=120, hold=None, c ='black', linewidth = 2, linestyle = '--', zorder = 1)
    plt.axhline(y=-15, xmin=0, xmax=120, hold=None, c ='black', linewidth = 2, linestyle = '--', zorder = 1)
    ax.axis('tight')# Tweak spacing to prevent clipping of ylabel
    plt.savefig(path_fig + str(datetime.today()) +'_'+ title +'.png', dpi = 100, bbox_inches = 'tight')

def get_data(c_files, path_class, type_val, interact):
    polar = ['CYS', 'GLN', 'ASN', 'THR', 'SER', 'HIS']
    arom = ['PHE', 'TYR', 'TRP']
    backborn = ['N','C','O','CA']
    l_interact = [interact, interact[1]+interact[0]]
    data =[]
    total = 0
    bb = 0
    for namefile in c_files:
        text_r = open_file(path_class, namefile +'.txt' )
        for line in text_r:
            line = line[:-1]
            lines = line.split(" ")
            lines = list(filter(None, lines))
            alla = lines[0]
            allb = lines[1]
            alla = alla.split("_")
            allb = allb.split("_")
            resnamea = alla[2]
            resnameb = allb[2]
            atoma = alla[1]
            atomb = allb[1]
            distance = lines[2]
            interact_file = lines[4]
            if interact_file != 'XX':
                #cutoff = get_cutoff(distance,interact)
                #if distance <= str(cutoff):
                    #print (distance, cutoff)
                if lines[4] == l_interact[1] or lines[4] == l_interact[0]:
                    if type_val == 'height':
                        value = lines[3]
                    if type_val == 'distance':
                        value = lines[2]
                    #if resnamea in arom and resnameb in arom and lines[2] > str(4.2):
                        #print namefile
                        #print line
                    #if atoma not in backborn and atomb not in backborn:
                        #bb += 1
                        #print lines[2]
                    data.append(float(value))
                #total += 1
    array_data = np.array(data)
    return array_data, total

def get_cutoff(distance, interact):
    interact = interact[:-1]
    l_interact = [interact, interact[1]+interact[0]]
    if 'HH' in l_interact:
        cutoff = 4.5
    if 'HA' in l_interact:
        cutoff = 5
    if 'HP' in l_interact:
        cutoff = 4
    if 'HC' in l_interact:
        cutoff = 4
    if 'AA' in l_interact:
        cutoff = 5
    if 'AP' in l_interact:
        cutoff = 4
    if 'AC' in l_interact:
        cutoff = 4
    if 'PP' in l_interact:
        cutoff = 4
    if 'PC' in l_interact:
        cutoff = 4
    if 'CC' in l_interact:
        cutoff = 4
    return cutoff

def get_naccess(path_sin, path_pdb, path_na, path):
    p_files = get_files (path_pdb, '.pdb')
    s_files = get_files (path_sin, '.pdb')
    print(os.getcwd())
    for p in p_files:
        call(['/bin/bash', '-i', '-c', 'naccess '+ path_pdb + p + '.pdb -h'])
    for s in s_files:
        call(['/bin/bash', '-i', '-c', 'naccess '+ path_sin + s + '.pdb -h'])
    move_files(path,path_na)

def move_files(path, path_na):
    origen = path
    destino = path_na
    rsa_files = get_files(path, '.rsa')
    for rsa in rsa_files:
        shutil.move(origen + rsa + '.rsa', destino + rsa + '.rsa')
    asa_files = get_files(path, '.asa')
    for asa in asa_files:
        shutil.move(origen + asa + '.asa', destino + asa + '.asa')
    log_files = get_files(path, '.log')
    for log in log_files:
        shutil.move(origen + log + '.log', destino + log + '.log')

def super_acc(path_na):
    rsa_files = get_files(path_na, '.rsa')
    l =[]
    l_final = []
    for rsa in rsa_files:
        l_rsa = rsa.split('_') #2RH1_A_A_ADRB2_HUMAN
        if len(l_rsa[2]) == 2:
            l.append([rsa + '.rsa', rsa[0:8]+rsa[9:] + '.rsa',rsa[0:7]+rsa[8:]+ '.rsa'])
    for group in l:
        fileA = open_file(path_na,group[1])
        fileB = open_file(path_na,group[2])
        filedim = open_file(path_na, group[0])
        for fa in fileA:
            if fa.startswith('TOTAL'):
                fa = fa.split(' ')
                fa = filter(None,fa)
                valueA = float(fa[1])
        for fb in fileB:
            if fb.startswith('TOTAL'):
                fb = fb.split(' ')
                fb = filter(None,fb)
                valueB = float(fb[1])
        for fdi in filedim:
            if fdi.startswith('TOTAL'):
                fdi = fdi.split(' ')
                fdi = filter(None,fdi)
                valuedim = float(fdi[1])
        total = valueA + valueB -valuedim
        if total < 0:
            total = 0.0
        l_final.append(total)
    return l_final

def create_table(lz,nolz):
    data_matrix = ['Name', 'Surface Acces LZ', 'Surface Acces No LZ']
    for l in range(len(lz)):
        if lz[l][0] == nolz[l][0]:
            print ('{0}, {1:.2f}, {2:.2f}'.format(lz[l][0], lz[l][1], nolz[l][1]))

def get_plot(l_1, l_2, path_fig):
    l_1 = np.array(l_1)
    l_2 = np.array(l_2)
    data = [l_1]
    title = "HISTOGRAMA DE SUPERFICIES D'INTERACCIO"
    fig, ax = plt.subplots(figsize =(15,10))
    colors = ['#FF0000','#00FFD5']
    interacts = ['Sf', 'Snf']
    ax.hist(data, bins = np.arange(0, 3600, 100), histtype = 'bar', stacked = True, color = colors[0], label = interacts[0], orientation = 'vertical', rwidth=1, alpha = 0.75)
    ax.set(title = title, ylabel = 'NOMBRE DE DIMERS', xlabel = 'SUPERFICIE ACCESSIBLE')
    ax.grid(False)
    ax.set_ylim([0,30])
    plt.xticks(np.arange(0,3500, 200))
    plt.yticks(np.arange(0, 30, 5))
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))## Put a legend to the right of the current axis
    plt.savefig(path_fig + str(datetime.today()) +'_'+ title +'.png', dpi = 100, bbox_inches = 'tight')

###############################################################################
# FILES
###############################################################################
path_pdb =  path + '/../PDB/dimers/'
path_lz = path + '/../PDB/dimers_lz/'
path_gro = path + '/../SIMULATIONS_'+ sys.argv[1] +'/'
path_opm = path + "/../PDB/OPM/"
path_out1 =  path + "/../RESULTS/Interactions_pdb/"
path_outlz = path + "/../RESULTS/Interactions_pdb_lz/"
path_out2 =  path + "/../RESULTS/Interactions_gro_"+ sys.argv[1] +"/"
path_fig = path + '/../RESULTS/Pictures/'
path_na = path + '/../RESULTS/Naccess/'
path_pdb_count = path + "/../RESULTS/COUNTS_pdb/"
path_pdblz_count = path + "/../RESULTS/COUNTS_pdb_lz/"
path_gro_count = path + "/../RESULTS/COUNTS_gro_"+ sys.argv[1] +"/"

#create_dir(path_opm)
create_dir(path + '/../RESULTS')
create_dir(path_out1)
create_dir(path_outlz)
create_dir(path_out2)
create_dir(path_fig)
create_dir(path_na)
create_dir(path_pdb_count)
create_dir(path_pdblz_count)
create_dir(path_gro_count)

pdb_files = get_files (path_pdb, ".pdb") #List of all infile names
#pdb_files = [pdb_files[13], pdb_files[18], pdb_files[22], pdb_files[29], pdb_files[43], pdb_files[44], pdb_files[45], pdb_files[46], pdb_files[52], pdb_files[53], pdb_files[54]]
o1_files = get_files (path_out1, ".txt") #List of all outfile names
o2_files = get_files (path_out2, ".txt") #List of all outfile names

###############################################################################
# GET ATOM SELECTION AND INTERACTIONS
###############################################################################
if __name__ == "__main__":
    selection = input("> Select the type of data that you want analyze (pdb/gro): ")
    i = 1
    driver = webdriver.Firefox()
    pdb_fail = ["3OE8", "4JKV", "5LWE", "5TUD"]
    if selection == 'pdb':
        lz = input("> Choose without fusion protein (0) or with fusion protein (1): ")
        if str(lz) == "1": 
            path_out1 = path_outlz
            path_pdb = path_lz
            path_pdb_count = path_pdblz_count
    for namefile in pdb_files:
        if selection == 'pdb':
            print (f'> {i *100 / len(pdb_files)} %')
            file_okay = glob.glob(f'{path_out1}/ints_{namefile}*.txt')
            if file_okay == []:
                pdbfile = namefile + ".pdb"
                print (f'> Reading {pdbfile}')
                pdb_r = open_file(path_pdb,pdbfile)#Save all pdb info into a variable
                wat_ions = ["HOH", "NA", "CL", "ZN", "SO4", "PO4"] #List of ligands that we don't want
                protein, ligand, chains, lipids, rest = get_all_pdb(pdb_r, wat_ions)#Separate the pdb file into protein and ligands, and get the name of all lipids
                pdb = namefile.split("_")
                print (f'> Getting OPM coordinates of {namefile}')
                #Obtain coords to difference the extra, trans, intra zone
                coords = get_opm_coords(pdb[0].lower())  
               
                u = mda.Universe(path_pdb + pdbfile)#Create a universe with the pdbfile (after select the part that it will analysed):
                ##INTERACCION PROT --> PROT:
                print ('> Getting interactions prot-prot!')
                chains = sorted(chains)
                prot_a = "protein and not (name H*) and segid "+ chains[0]
                prot_b = "protein and not (name H*) and segid "+ chains[1] + " and around 10 protein and not (name H*) and segid "+ chains[0]
                if pdb[0] == "4QIM" or (pdb[0] in pdb_fail and lz =="1"):
                    prot_b = "protein and not (name H*) and segid "+ chains[1] + " and around 11 protein and not (name H*) and segid "+ chains[0]
                try:
                    save_interac_pdb(u, prot_a, prot_b, path_out1, namefile, coords)
                except:
                    prot_b = "protein and segid "+ chains[1] + " and around 10 protein and segid "+ chains[0]
                    if pdb[0] == "4QIM" or (pdb[0] in pdb_fail and lz =="1"):
                        prot_b = "protein and segid "+ chains[1] + " and around 11 protein and segid "+ chains[0]
                    save_interac_pdb(u, prot_a, prot_b, path_out1, namefile, coords)
                #INTERACCIONS PROT --> LIP:
                print ('> Getting interactions prot-lip!')
                if lipids != []: #If the file don't have lipids, as from here don't continue
                    prot = "protein and not (name H*)"
                    lip = "resname {}".format(" + ".join(lipids)) + " and around 10 protein"
                    if pdb[0] == "4QIM" or (pdb[0] in pdb_fail and lz =="1"):
                        lip = "resname {}".format(" + ".join(lipids)) + " and around 11 protein"
                    save_interac_pdb(u, prot, lip, path_out1, namefile, coords) 
        ##################################################################################################################################
        if selection == 'gro':
            print (f'> {i *100 / len(pdb_files)} %')
            file_okay = glob.glob(f'{path_out2}/ints_{namefile}*.txt')
            if file_okay == []:
                path_gro_dimer = f'{path_gro}{namefile}/'
                grofiles = [f'{namefile}_ini.gro', f'{namefile}_p_md.gro']
                #grofiles = [f'{namefile}_p_md.gro']
                for grofile in grofiles:
                    wat_ions = ["NA+", "CL-"]
                    if 'p_md.gro' in grofile and os.path.exists(path_gro_dimer + f'{namefile}_p_center.gro'):
                        grofile = f'{namefile}_p_center.gro' #That dimers that have a bad last frame (prova.sh)
                    print (f'> Reading {grofile}')
                    gro_r = open_file(path_gro_dimer, grofile)#Save all pdb info into a variable
                    protein, popc, w, ion, rest = get_all_gro(gro_r, wat_ions)#Separate the pdb file into protein and ligands, and get the name of all lipids
                    if '_ini.gro' in grofile:
                        ref_r = gro_r
                        print (f'> Converting {grofile} to pdb')
                        os.system(f"gmx editconf -f {path_gro_dimer}{namefile}_ini.gro -o {path_gro_dimer}{namefile}_ini.pdb")
                    if '_p_md.gro' in grofile:
                        print (f'> Rename {grofile} to 2_{grofile}')
                        rename_residues(gro_r, path_gro_dimer, grofile, ref_r)
                        print (f'> Converting 2_{grofile} to pdb')
                        os.system(f"gmx editconf -f {path_gro_dimer}2_{namefile}_p_md.gro -o {path_gro_dimer}2_{namefile}_p_md.pdb")
                    if '_p_center.gro' in grofile:
                        print (f'> Rename {grofile} to 2_{grofile}')
                        rename_residues(gro_r, path_gro_dimer, grofile, ref_r)
                        print (f'> Converting 2_{grofile} to pdb')
                        os.system(f"gmx editconf -f {path_gro_dimer}2_{namefile}_p_center.gro -o {path_gro_dimer}2_{namefile}_p_center.pdb")
                    try:#Obtain coords to difference the extra, trans, intra zone
                        coords = get_opm_coords(namefile)
                        print (f'> Getting OPM coordinates of {namefile}')
                    except: #In case that code name of protein not found it on opm, get the representative one.
                        pdb = namefile.split("_")
                        print (f'> Getting OPM coordinates of {namefile}')
                        #Obtain coords to difference the extra, trans, intra zone
                        coords = get_opm_coords(pdb[0].lower())  
                    if '_ini.gro' in grofile:
                        u = mda.Universe( f'{path_gro_dimer}{namefile}_ini.pdb')#Create a universe with the pdbfile (after select the part that it will analysed):
                    if '_p_md.gro' in grofile:
                        u = mda.Universe( f'{path_gro_dimer}2_{namefile}_p_md.pdb')#Create a universe with the pdbfile (after select the part that it will analysed):
                    if '_p_center.gro' in grofile:
                        u = mda.Universe( f'{path_gro_dimer}2_{namefile}_p_center.pdb')
                    #resnum = re.sub('[POPC]', '', popc[0][0])
                    #resnum = int(resnum)-1
                    atomnum = re.sub('[POPC]', '', popc[0][2])
                    atomnum = int(atomnum)-1
                    start2 = 0
                    for p in protein:
                        act = int(p[0][:-3])
                        if act > start2:
                            start2 = act
                        if act < start2:
                            start2 = int(p[2])
                            break
                    end1 = start2-1
                    ##INTERACCION PROT --> PROT:
                    print ('> Getting interactions prot-prot!')
                    prot_a = u.atoms[:int(end1)]
                    prot_b = u.atoms[int(start2)-1:int(atomnum)]
                    save_interac_gro(u, prot_a, prot_b, path_out2, grofile[:-4], coords)

        ########################################################################
        if selection != 'pdb' and selection != 'gro':
            print ('> Choose a type of data: pdb/gro!')
            print ('> Restarting!')
            driver.close()
            sys.exit()
        i += 1
    driver.close()
    if selection == 'pdb':
        ###############################################################################
        # GET HEATMAPS
        ###############################################################################
        ###############################################################################
        # DIRECT INTERACTIONS
        ###############################################################################
        select = [ 'direct', 'ligand'] #'direct' or 'ligand'
        print ('> Getting direct interactions')
        for sel in select:
            direct = get_numpy_direct(o1_files, path_out1, '.txt', sel)
            if sel == 'ligand':
                title = "| HEATMAP LIGAND INTERACTIONS | "
            if sel == 'direct':
                title = "| HEATMAP DIRECT INTERACTIONS | "
            get_heatmap(direct, title, path_fig)
        ###############################################################################
        # INDIRECT INTERACTIONS
        ###############################################################################
        indirect = get_numpy_indirect(o1_files, path_out1, '.txt')
        print ('> Getting indirect interactions')
        title = "| HEATMAP INDIRECT INTERACTIONS | "
        get_heatmap(indirect, title, path_fig)
        ###############################################################################
        # TYPES OF INTERACTIONS AND DISTRIBUTION ON THE MEMBRANE
        ###############################################################################
        interacts = ['HH', 'HA', 'HP', 'HC', 'AA','AP','AC','PP', 'PC', 'CC']
        print ('> Getting the type of interactions')
        #interacts = ['CC']
        for interact in interacts:
            get_hist_dist(o1_files, path_out1, interact, path_fig)
        get_hist_z(o1_files, path_out1, path_fig)
        ###############################################################################
        # NACCES SURFACE
        ###############################################################################
        #get_naccess(path_sinlz, path_pdblz, path_nalz, path)
        #na = super_acc(path_na)
        ###create_table(lz,nolz)
        #print ('> Getting the surface analysis')
        #get_plot(lz, nolz, path_fig)
        #move_files(path_pdblz, path_nalz)

        ###############################################################################
        # RESULTS COUNT
        ###############################################################################
        ###############################################################################
        # EACH DIMER
        ###############################################################################
        inter = get_files (path_out1, '.txt')
        print ('> Getting the count analysis for each dimer')
        for i in inter:
            d_counts={'TM7\tTM5': 0, 'TM7\tTM3': 0, 'TM4\tTM6': 0, 'NT\tTM2': 0, 'ICL2\tTM1': 0, 'CT\tTM6': 0, 'CT\tTM3': 0, 'ECL1\tTM6': 0, 'ICL1\tCT': 0, 'ICL1\tICL1': 0, 'CT\tCT': 0, 'NT\tTM1': 0, 'TM1\tTM1': 0, 'ICL2\tICL3': 0, 'ICL2\tICL2': 0, 'TM5\tTM5': 0, 'TM5\tICL3': 0, 'ICL3\tICL3': 0, 'TM4\tTM5': 0, 'TM4\tTM4': 0, 'ECL2\tECL2': 0, 'ECL2\tTM5': 0, 'ECL2\tECL3': 0, 'TM5\tECL3': 0, 'NT\tECL3': 0, 'NT\tTM7': 0, 'NT\tTM6': 0, 'TM1\tTM6': 0, 'TM1\tECL3': 0, 'TM1\tTM7': 0, 'TM1\tTM5': 0, 'TM1\tICL3': 0, 'ICL1\tICL3': 0, 'ICL2\tICL1': 0, 'TM4\tICL1': 0, 'TM4\tTM1': 0, 'TM4\tTM2': 0, 'TM4\tTM3': 0, 'TM4\tECL1': 0, 'ECL2\tECL1': 0, 'TM5\tTM2': 0, 'TM5\tECL1': 0, 'ECL1\tTM1': 0, 'CT\tICL3': 0, 'CT\tTM5': 0, 'TM3\tECL2': 0, 'TM4\tNT': 0, 'NT\tNT': 0, 'NT\tECL2': 0, 'NT\tECL1': 0, 'NT\tTM3': 0, 'TM1\tTM2': 0, 'ECL1\tECL1': 0, 'TM7\tCT': 0, 'TM2\tTM2': 0, 'TM2\tECL1': 0, 'ECL2\tTM6': 0, 'TM5\tTM6': 0, 'ICL3\tTM6': 0, 'TM6\tTM6': 0, 'ICL3\tTM7': 0, 'TM6\tTM7': 0, 'TM6\tECL3': 0, 'ECL3\tECL3': 0, 'NT\tTM5': 0, 'TM1\tECL2': 0, 'TM7\tTM4': 0, 'CT\tTM4': 0, 'CT\tICL2': 0, 'TM1\tICL1': 0, 'TM1\tTM3': 0, 'TM2\tTM7': 0, 'ECL1\tTM7': 0, 'ECL1\tECL3': 0, 'TM2\tTM3': 0, 'ECL1\tTM3': 0, 'TM5\tICL2': 0, 'ICL3\tTM4': 0, 'ICL3\tTM3': 0, 'TM7\tECL3': 0, 'TM7\tECL2': 0, 'TM7\tICL2': 0, 'ECL3\tTM3': 0, 'ECL2\tTM2': 0, 'TM1\tCT': 0, 'TM3\tTM3': 0, 'ICL1\tTM5': 0, 'ICL2\tTM4': 0, 'TM4\tECL2': 0, 'TM3\tICL2': 0, 'TM7\tTM7': 0}
            total = 0
            print (f'> Opening {i} file.')
            info = open_file (path_out1,i+'.txt')
            print (f'> Counting...')
            for j in info:
                line = j.split(' ')
                if line[-1][:2] != 'XX':
                    info_a = line[0].split('_')
                    info_b = line[1].split('_')
                    if not info_a[-1] + '\t' + info_b[-1] in d_counts.keys() and not info_b[-1] + '\t' + info_a[-1] in d_counts.keys():
                        d_counts [info_a[-1] + '\t' + info_b[-1]] = 1
                    else:
                        try:
                            d_counts[info_a[-1] + '\t' + info_b[-1]] = d_counts[info_a[-1] + '\t' + info_b[-1]] + 1
                        except:
                            d_counts[info_b[-1] + '\t' + info_a[-1]] = d_counts[info_b[-1] + '\t' + info_a[-1]] + 1
                    total += 1
            count_inter = open(f'{path_pdb_count}count_inter_{i}.txt', "w")
            count_inter.writelines(f'{total}\n') #All info here
            d_counts = [(k,d_counts[k]) for k in sorted(d_counts, key=d_counts.get,reverse =True)]
            for k, v in d_counts:
                if v != 0:
                    count_inter.writelines(f'{k}\t{v}\t{v*100/total}\n') #All info here
                else:
                    count_inter.writelines(f'{k}\t{v}\t0\n') #All info here
            count_inter.close()
        ################################################################################
        ## TOTAL
        ################################################################################
        inter = get_files (path_out1, '.txt')
        print ('> Getting the total count analysis')
        d_counts={'TM7\tTM5': 0, 'TM7\tTM3': 0, 'TM4\tTM6': 0, 'NT\tTM2': 0, 'ICL2\tTM1': 0, 'CT\tTM6': 0, 'CT\tTM3': 0, 'ECL1\tTM6': 0, 'ICL1\tCT': 0, 'ICL1\tICL1': 0, 'CT\tCT': 0, 'NT\tTM1': 0, 'TM1\tTM1': 0, 'ICL2\tICL3': 0, 'ICL2\tICL2': 0, 'TM5\tTM5': 0, 'TM5\tICL3': 0, 'ICL3\tICL3': 0, 'TM4\tTM5': 0, 'TM4\tTM4': 0, 'ECL2\tECL2': 0, 'ECL2\tTM5': 0, 'ECL2\tECL3': 0, 'TM5\tECL3': 0, 'NT\tECL3': 0, 'NT\tTM7': 0, 'NT\tTM6': 0, 'TM1\tTM6': 0, 'TM1\tECL3': 0, 'TM1\tTM7': 0, 'TM1\tTM5': 0, 'TM1\tICL3': 0, 'ICL1\tICL3': 0, 'ICL2\tICL1': 0, 'TM4\tICL1': 0, 'TM4\tTM1': 0, 'TM4\tTM2': 0, 'TM4\tTM3': 0, 'TM4\tECL1': 0, 'ECL2\tECL1': 0, 'TM5\tTM2': 0, 'TM5\tECL1': 0, 'ECL1\tTM1': 0, 'CT\tICL3': 0, 'CT\tTM5': 0, 'TM3\tECL2': 0, 'TM4\tNT': 0, 'NT\tNT': 0, 'NT\tECL2': 0, 'NT\tECL1': 0, 'NT\tTM3': 0, 'TM1\tTM2': 0, 'ECL1\tECL1': 0, 'TM7\tCT': 0, 'TM2\tTM2': 0, 'TM2\tECL1': 0, 'ECL2\tTM6': 0, 'TM5\tTM6': 0, 'ICL3\tTM6': 0, 'TM6\tTM6': 0, 'ICL3\tTM7': 0, 'TM6\tTM7': 0, 'TM6\tECL3': 0, 'ECL3\tECL3': 0, 'NT\tTM5': 0, 'TM1\tECL2': 0, 'TM7\tTM4': 0, 'CT\tTM4': 0, 'CT\tICL2': 0, 'TM1\tICL1': 0, 'TM1\tTM3': 0, 'TM2\tTM7': 0, 'ECL1\tTM7': 0, 'ECL1\tECL3': 0, 'TM2\tTM3': 0, 'ECL1\tTM3': 0, 'TM5\tICL2': 0, 'ICL3\tTM4': 0, 'ICL3\tTM3': 0, 'TM7\tECL3': 0, 'TM7\tECL2': 0, 'TM7\tICL2': 0, 'ECL3\tTM3': 0, 'ECL2\tTM2': 0, 'TM1\tCT': 0, 'TM3\tTM3': 0, 'ICL1\tTM5': 0, 'ICL2\tTM4': 0, 'TM4\tECL2': 0, 'TM3\tICL2': 0, 'TM7\tTM7': 0}
        total = 0
        for i in inter:
            print (f'> Opening {i} file.')
            info = open_file (path_out1,i+'.txt')
            print (f'> Counting...')
            for j in info:
                line = j.split(' ')
                if line[-1][:2] != 'XX':
                    info_a = line[0].split('_')
                    info_b = line[1].split('_')
                    if not info_a[-1] + '\t' + info_b[-1] in d_counts.keys() and not info_b[-1] + '\t' + info_a[-1] in d_counts.keys():
                        d_counts [info_a[-1] + '\t' + info_b[-1]] = 1
                    else:
                        try:
                            d_counts[info_a[-1] + '\t' + info_b[-1]] = d_counts[info_a[-1] + '\t' + info_b[-1]] + 1
                        except:
                            d_counts[info_b[-1] + '\t' + info_a[-1]] = d_counts[info_b[-1] + '\t' + info_a[-1]] + 1
                    total += 1
        count_inter = open(f'{path_pdb_count}count_inter.txt', "w")
        count_inter.writelines(f'{total}\n') #All info here
        d_counts = [(k, d_counts[k]) for k in sorted(d_counts, key=d_counts.get,reverse =True)]
        for k, v in d_counts:
            if v != 0:
                count_inter.writelines(f'{k}\t{v}\t{v*100/total}\n') #All info here
            else:
                count_inter.writelines(f'{k}\t{v}\t0\n') #All info here
        count_inter.close()
    ###########################################################################
    if selection == 'gro':
        ###############################################################################
        # RESULTS COUNT
        ###############################################################################
        ###############################################################################
        # EACH DIMER
        ###############################################################################
        inter = get_files (path_out2, '.txt')
        print ('> Getting the count analysis for each dimer')
        for i in inter:
            d_counts={'TM7\tTM5': 0, 'TM7\tTM3': 0, 'TM4\tTM6': 0, 'NT\tTM2': 0, 'ICL2\tTM1': 0, 'CT\tTM6': 0, 'CT\tTM3': 0, 'ECL1\tTM6': 0, 'ICL1\tCT': 0, 'ICL1\tICL1': 0, 'CT\tCT': 0, 'NT\tTM1': 0, 'TM1\tTM1': 0, 'ICL2\tICL3': 0, 'ICL2\tICL2': 0, 'TM5\tTM5': 0, 'TM5\tICL3': 0, 'ICL3\tICL3': 0, 'TM4\tTM5': 0, 'TM4\tTM4': 0, 'ECL2\tECL2': 0, 'ECL2\tTM5': 0, 'ECL2\tECL3': 0, 'TM5\tECL3': 0, 'NT\tECL3': 0, 'NT\tTM7': 0, 'NT\tTM6': 0, 'TM1\tTM6': 0, 'TM1\tECL3': 0, 'TM1\tTM7': 0, 'TM1\tTM5': 0, 'TM1\tICL3': 0, 'ICL1\tICL3': 0, 'ICL2\tICL1': 0, 'TM4\tICL1': 0, 'TM4\tTM1': 0, 'TM4\tTM2': 0, 'TM4\tTM3': 0, 'TM4\tECL1': 0, 'ECL2\tECL1': 0, 'TM5\tTM2': 0, 'TM5\tECL1': 0, 'ECL1\tTM1': 0, 'CT\tICL3': 0, 'CT\tTM5': 0, 'TM3\tECL2': 0, 'TM4\tNT': 0, 'NT\tNT': 0, 'NT\tECL2': 0, 'NT\tECL1': 0, 'NT\tTM3': 0, 'TM1\tTM2': 0, 'ECL1\tECL1': 0, 'TM7\tCT': 0, 'TM2\tTM2': 0, 'TM2\tECL1': 0, 'ECL2\tTM6': 0, 'TM5\tTM6': 0, 'ICL3\tTM6': 0, 'TM6\tTM6': 0, 'ICL3\tTM7': 0, 'TM6\tTM7': 0, 'TM6\tECL3': 0, 'ECL3\tECL3': 0, 'NT\tTM5': 0, 'TM1\tECL2': 0, 'TM7\tTM4': 0, 'CT\tTM4': 0, 'CT\tICL2': 0, 'TM1\tICL1': 0, 'TM1\tTM3': 0, 'TM2\tTM7': 0, 'ECL1\tTM7': 0, 'ECL1\tECL3': 0, 'TM2\tTM3': 0, 'ECL1\tTM3': 0, 'TM5\tICL2': 0, 'ICL3\tTM4': 0, 'ICL3\tTM3': 0, 'TM7\tECL3': 0, 'TM7\tECL2': 0, 'TM7\tICL2': 0, 'ECL3\tTM3': 0, 'ECL2\tTM2': 0, 'TM1\tCT': 0, 'TM3\tTM3': 0, 'ICL1\tTM5': 0, 'ICL2\tTM4': 0, 'TM4\tECL2': 0, 'TM3\tICL2': 0, 'TM7\tTM7': 0}
            total = 0
            print (f'> Opening {i} file.')
            info = open_file (path_out2,i+'.txt')
            print (f'> Counting...')
            for j in info:
                line = j.split(' ')
                if line[-1][:2] != 'XX':
                    info_a = line[0].split('_')
                    info_b = line[1].split('_')
                    if not info_a[-1] + '\t' + info_b[-1] in d_counts.keys() and not info_b[-1] + '\t' + info_a[-1] in d_counts.keys():
                        d_counts [info_a[-1] + '\t' + info_b[-1]] = 1
                    else:
                        try:
                            d_counts[info_a[-1] + '\t' + info_b[-1]] = d_counts[info_a[-1] + '\t' + info_b[-1]] + 1
                        except:
                            d_counts[info_b[-1] + '\t' + info_a[-1]] = d_counts[info_b[-1] + '\t' + info_a[-1]] + 1
                    total += 1
            count_inter = open(f'{path_gro_count}count_inter_{i}.txt', "w")
            count_inter.writelines(f'{total}\n') #All info here
            d_counts = [(k,d_counts[k]) for k in sorted(d_counts, key=d_counts.get,reverse =True)]
            for k, v in d_counts:
                if v != 0:
                    count_inter.writelines(f'{k}\t{v}\t{v*100/total}\n') #All info here
                else:
                    count_inter.writelines(f'{k}\t{v}\t0\n') #All info here
            count_inter.close()
        ###############################################################################
        # TOTAL
        ###############################################################################
        inter = get_files (path_out2, '.txt')
        print ('> Getting the total count analysis')
        d_counts_ini = {'TM7\tTM5': 0, 'TM7\tTM3': 0, 'TM4\tTM6': 0, 'NT\tTM2': 0, 'ICL2\tTM1': 0, 'CT\tTM6': 0, 'CT\tTM3': 0, 'ECL1\tTM6': 0, 'ICL1\tCT': 0, 'ICL1\tICL1': 0, 'CT\tCT': 0, 'NT\tTM1': 0, 'TM1\tTM1': 0, 'ICL2\tICL3': 0, 'ICL2\tICL2': 0, 'TM5\tTM5': 0, 'TM5\tICL3': 0, 'ICL3\tICL3': 0, 'TM4\tTM5': 0, 'TM4\tTM4': 0, 'ECL2\tECL2': 0, 'ECL2\tTM5': 0, 'ECL2\tECL3': 0, 'TM5\tECL3': 0, 'NT\tECL3': 0, 'NT\tTM7': 0, 'NT\tTM6': 0, 'TM1\tTM6': 0, 'TM1\tECL3': 0, 'TM1\tTM7': 0, 'TM1\tTM5': 0, 'TM1\tICL3': 0, 'ICL1\tICL3': 0, 'ICL2\tICL1': 0, 'TM4\tICL1': 0, 'TM4\tTM1': 0, 'TM4\tTM2': 0, 'TM4\tTM3': 0, 'TM4\tECL1': 0, 'ECL2\tECL1': 0, 'TM5\tTM2': 0, 'TM5\tECL1': 0, 'ECL1\tTM1': 0, 'CT\tICL3': 0, 'CT\tTM5': 0, 'TM3\tECL2': 0, 'TM4\tNT': 0, 'NT\tNT': 0, 'NT\tECL2': 0, 'NT\tECL1': 0, 'NT\tTM3': 0, 'TM1\tTM2': 0, 'ECL1\tECL1': 0, 'TM7\tCT': 0, 'TM2\tTM2': 0, 'TM2\tECL1': 0, 'ECL2\tTM6': 0, 'TM5\tTM6': 0, 'ICL3\tTM6': 0, 'TM6\tTM6': 0, 'ICL3\tTM7': 0, 'TM6\tTM7': 0, 'TM6\tECL3': 0, 'ECL3\tECL3': 0, 'NT\tTM5': 0, 'TM1\tECL2': 0, 'TM7\tTM4': 0, 'CT\tTM4': 0, 'CT\tICL2': 0, 'TM1\tICL1': 0, 'TM1\tTM3': 0, 'TM2\tTM7': 0, 'ECL1\tTM7': 0, 'ECL1\tECL3': 0, 'TM2\tTM3': 0, 'ECL1\tTM3': 0, 'TM5\tICL2': 0, 'ICL3\tTM4': 0, 'ICL3\tTM3': 0, 'TM7\tECL3': 0, 'TM7\tECL2': 0, 'TM7\tICL2': 0, 'ECL3\tTM3': 0, 'ECL2\tTM2': 0, 'TM1\tCT': 0, 'TM3\tTM3': 0, 'ICL1\tTM5': 0, 'ICL2\tTM4': 0, 'TM4\tECL2': 0, 'TM3\tICL2': 0, 'TM7\tTM7': 0}
        d_counts_final = {'TM7\tTM5': 0, 'TM7\tTM3': 0, 'TM4\tTM6': 0, 'NT\tTM2': 0, 'ICL2\tTM1': 0, 'CT\tTM6': 0, 'CT\tTM3': 0, 'ECL1\tTM6': 0, 'ICL1\tCT': 0, 'ICL1\tICL1': 0, 'CT\tCT': 0, 'NT\tTM1': 0, 'TM1\tTM1': 0, 'ICL2\tICL3': 0, 'ICL2\tICL2': 0, 'TM5\tTM5': 0, 'TM5\tICL3': 0, 'ICL3\tICL3': 0, 'TM4\tTM5': 0, 'TM4\tTM4': 0, 'ECL2\tECL2': 0, 'ECL2\tTM5': 0, 'ECL2\tECL3': 0, 'TM5\tECL3': 0, 'NT\tECL3': 0, 'NT\tTM7': 0, 'NT\tTM6': 0, 'TM1\tTM6': 0, 'TM1\tECL3': 0, 'TM1\tTM7': 0, 'TM1\tTM5': 0, 'TM1\tICL3': 0, 'ICL1\tICL3': 0, 'ICL2\tICL1': 0, 'TM4\tICL1': 0, 'TM4\tTM1': 0, 'TM4\tTM2': 0, 'TM4\tTM3': 0, 'TM4\tECL1': 0, 'ECL2\tECL1': 0, 'TM5\tTM2': 0, 'TM5\tECL1': 0, 'ECL1\tTM1': 0, 'CT\tICL3': 0, 'CT\tTM5': 0, 'TM3\tECL2': 0, 'TM4\tNT': 0, 'NT\tNT': 0, 'NT\tECL2': 0, 'NT\tECL1': 0, 'NT\tTM3': 0, 'TM1\tTM2': 0, 'ECL1\tECL1': 0, 'TM7\tCT': 0, 'TM2\tTM2': 0, 'TM2\tECL1': 0, 'ECL2\tTM6': 0, 'TM5\tTM6': 0, 'ICL3\tTM6': 0, 'TM6\tTM6': 0, 'ICL3\tTM7': 0, 'TM6\tTM7': 0, 'TM6\tECL3': 0, 'ECL3\tECL3': 0, 'NT\tTM5': 0, 'TM1\tECL2': 0, 'TM7\tTM4': 0, 'CT\tTM4': 0, 'CT\tICL2': 0, 'TM1\tICL1': 0, 'TM1\tTM3': 0, 'TM2\tTM7': 0, 'ECL1\tTM7': 0, 'ECL1\tECL3': 0, 'TM2\tTM3': 0, 'ECL1\tTM3': 0, 'TM5\tICL2': 0, 'ICL3\tTM4': 0, 'ICL3\tTM3': 0, 'TM7\tECL3': 0, 'TM7\tECL2': 0, 'TM7\tICL2': 0, 'ECL3\tTM3': 0, 'ECL2\tTM2': 0, 'TM1\tCT': 0, 'TM3\tTM3': 0, 'ICL1\tTM5': 0, 'ICL2\tTM4': 0, 'TM4\tECL2': 0, 'TM3\tICL2': 0, 'TM7\tTM7': 0}
        total_1 = 0
        total_2 = 0
        for i in inter:
            print (f'> Opening {i} file.')
            info = open_file (path_out2,i+'.txt')
            print (f'> Counting...')
            if 'ini' in i:
                for j in info:
                    line = j.split(' ')
                    if line[-1][:2] != 'XX':
                        info_a = line[0].split('_')
                        info_b = line[1].split('_')
                        if not info_a[-1] + '\t' + info_b[-1] in d_counts_ini.keys() and not info_b[-1] + '\t' + info_a[-1] in d_counts_ini.keys():
                            d_counts_ini [info_a[-1] + '\t' + info_b[-1]] = 1
                        else:
                            try:
                                d_counts_ini[info_a[-1] + '\t' + info_b[-1]] = d_counts_ini[info_a[-1] + '\t' + info_b[-1]] + 1
                            except:
                                d_counts_ini[info_b[-1] + '\t' + info_a[-1]] = d_counts_ini[info_b[-1] + '\t' + info_a[-1]] + 1
                        total_1 += 1
            if 'ini' not in i:
                for j in info:
                    line = j.split(' ')
                    if line[-1][:2] != 'XX':
                        info_a = line[0].split('_')
                        info_b = line[1].split('_')
                        if not info_a[-1] + '\t' + info_b[-1] in d_counts_final.keys() and not info_b[-1] + '\t' + info_a[-1] in d_counts_final.keys():
                            d_counts_final [info_a[-1] + '\t' + info_b[-1]] = 1
                        else:
                            try:
                                d_counts_final[info_a[-1] + '\t' + info_b[-1]] = d_counts_final[info_a[-1] + '\t' + info_b[-1]] + 1
                            except:
                                d_counts_final[info_b[-1] + '\t' + info_a[-1]] = d_counts_final[info_b[-1] + '\t' + info_a[-1]] + 1
                        total_2 += 1

        count_inter_ini = open(f'{path_gro_count}count_inter_ini.txt', "w")
        count_inter_ini.writelines(f'{total_1}\n') #All info here
        d_counts_ini = [(k, d_counts_ini[k]) for k in sorted(d_counts_ini, key=d_counts_ini.get,reverse =True)]
        for k, v in d_counts_ini:
            if v != 0:
                count_inter_ini.writelines(f'{k}\t{v}\t{v*100/total_1}\n') #All info here
            else:
                count_inter_ini.writelines(f'{k}\t{v}\t0\n') #All info here
        count_inter_ini.close()

        count_inter_end = open(f'{path_gro_count}count_inter_end.txt', "w")
        count_inter_end.writelines(f'{total_2}\n') #All info here
        d_counts_final = [(k, d_counts_final[k]) for k in sorted(d_counts_final, key=d_counts_final.get,reverse =True)]
        for k, v in d_counts_final:
            if v != 0:
                count_inter_end.writelines(f'{k}\t{v}\t{v*100/total_2}\n') #All info here
            else:
                count_inter_end.writelines(f'{k}\t{v}\t0\n') #All info here
        count_inter_end.close()


