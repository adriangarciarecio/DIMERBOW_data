#!/usr/bin/env python3


from bs4 import BeautifulSoup
import numpy as np
from urllib.request import urlopen
import os

def get_all_gpcrdb_structures():
    '''Returns the list of PDB codes from GPCRDB.
    The first row contains field_names'''
    resp = urlopen('http://gpcrdb.org/structure')
    soup = BeautifulSoup(resp.read(), "html5lib")
    table = soup.find('table')
    rows = table.find_all('tr')
    results = []
    for row in rows:
        table_headers = row.find_all('th')
        if table_headers:
            results.append([headers.get_text() for headers in table_headers])

        table_data = row.find_all('td')
        if table_data:
            results.append([data.get_text().strip() for data in table_data])
    # some fields with multiple entries (ligands ...) contain '\n'

    for i, pdb in enumerate(results):
        for j, column in enumerate(pdb):
            if '\n' in column:
                #print(column)
                results[i][j] = results[i][j].replace('\n', ',').replace(' ', '').replace(',,', ', ')
    # field_names = results.pop(0)  # first row has titles
    field_names = results[1]
    field_names.pop(0)
    field_names.pop(-1)
    print(field_names)
    results.pop(0)
    results.pop(1)

    for i in range(len(results)):
        results[i].pop(0)  # remove empty column (selection)
        results[i].pop(-1) # remove the 

    fields = np.array(field_names)
    # fields[i,:] --> item i
    # fields[:,f] -->  field f
    fields_gpcrdb = np.vstack([fields, results])
    return fields_gpcrdb


fields_gpcrdb = get_all_gpcrdb_structures()

#for i in range(len(fields_gpcrdb)):
#    print(fields_gpcrdb[i])


from mysql.connector import connect
import sys
# Connect
try:
    cnx = connect(host='localhost', user='root', password='Cair0!',
        database='dimerbow')
except:
    print('could not connect')
    sys.exit(0)

update_fields = ['ligand', 'ligand_type']
# a dictionary field_names --> position in column
fnames2position = {}
for i, name in enumerate(fields_gpcrdb[0]):
    fnames2position[name] = i


for update_field in update_fields:
    lmcdb2gpcrdb = {'ligand': 'X-ray ligand', 'ligand_type':'X-ray lig. function' }
    # Get pdbids with empty field (current_field)
    cursor = cnx.cursor()
    sql_query = "select pdbid, uniprot, resolution, resgroup from pdb where {} is NULL;".format(update_field)
    cursor.execute(sql_query)
    pdbs_lmcdb = []
    for fields_lmcdb in cursor:
        pdbs_lmcdb.append(fields_lmcdb[0])
    # Update current fields
    for pdb_lmcdb in pdbs_lmcdb:
        for i in range(len(fields_gpcrdb)):
            # IN CASE OF PROBLEMS it is possible that they have changed the format
            pdb_gpcrdb = fields_gpcrdb[i][4]   # OLD fields_gpcrdb[i,:][3]
            if pdb_gpcrdb == pdb_lmcdb:
                cur = cnx.cursor()
                new_value = fields_gpcrdb[i][fnames2position[lmcdb2gpcrdb[update_field]]]
                if 'pubchem' in new_value:
                    new_value = new_value.replace(', pubchem', '') # added to remove the pubchem link in field ligand_name
                #if "'" in new_value:
                #    new_value = new_value.replace("'", '') # added to remove the pubchem link in field ligand_name
                #UPDATE pdb set ligand = 'TAK-875' where pdbid = '4PHU';
                sql_query = "UPDATE pdb set {} = '{}' where pdbid = '{}'".format(update_field, new_value, pdb_lmcdb)
                print(sql_query)
                cur.execute(sql_query)
ask_commit = input('Commit changes? (y/n)')
if ask_commit == 'y':
    cnx.commit()
