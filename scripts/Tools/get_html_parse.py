#!/usr/bin/python3

###############################################################################
# IMPORTS
###############################################################################
import os
import sys
import glob
import shutil
import requests
import re

path = os.path.dirname(os.path.abspath("__file__"))
# sys.path.insert(0, path + "/..")
# print(path)
# from get_interactions import get_opm_coords, get_files
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time

###############################################################################


def get_opm_coords(pdb):
    """Get the coordinates to determine the differents zones on the structure of a protein."""
    opm_coords = (
        ""
    )  # Creates a dictionary for call each coordinates, extracted from opm, according to helix number
    driver.get("https://opm.phar.umich.edu/protein.php?search=" + pdb)
    time.sleep(2)
    all_info = driver.find_element_by_class_name("protein-content")
    table = all_info.find_element_by_class_name("small-text.break")
    coordinates = table.text
    coordinates = coordinates.split(":")
    coords = coordinates[2]
    coords = coords.split(",")
    for seg in coords:
        c = seg.partition("(")[-1].rpartition(")")[0]
        c = c.replace(" ", "")
        opm_coords = opm_coords + " " + c
    return opm_coords


def get_files(path, ext):
    """Obtain a list that contain the name of all pdb file in a directory"""
    l_files = []
    folder = glob.glob(path + "*" + ext)
    for path_filename in folder:
        path_filename = path_filename.split(
            "/"
        )  # ['', 'home', 'adrian', 'Escritorio', 'TFG', 'dimerInteractions', '2Z73_B_OPSD_TODPA.pse']
        filename = path_filename[-1]
        filename = filename.split(".")  # 4WW3_A_OPSD_TODPA.pse
        filename = filename[0]  # 4WW3_A_OPSD_TODPA
        l_files.append(filename)
    return sorted(l_files)


###############################################################################
# FILES
###############################################################################

###############################################################################
driver = webdriver.Firefox()
driver.get("https://gnomad.broadinstitute.org/variant/1-55516888-G-GA")
time.sleep(2)
all_info = driver.find_element_by_class_name("List-sc-19igk48-0.hztfx")
text = all_info.text
text = text.split(")")
print(text[0][7:])
driver.close()

