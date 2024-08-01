# MAIN IMPORT
import requests
import pandas as pd
import json
import os 
from time import sleep
import sys

path = os.path.dirname(os.path.abspath("__file__"))
sys.path.append(os.path.dirname(path))

# IMPORTS
from SECRETS import *
from gen_distribution import start_pymol, quit_pymol, load_pymol
from pymol import cmd
from pypdb import *


# DATABASE IMPORTS
from sqlalchemy import create_engine


def update_gpcrdb_info():
    
    # Get data from GPCRdb 
    print("- Getting information from GPCRdb")
    strucgpcrdb=requests.get('https://gpcrdb.org/services/structure/').json()
    gpcrdb_info = pd.DataFrame.from_dict(pd.json_normalize(strucgpcrdb), orient='columns')
    gpcrdb_info = gpcrdb_info.sort_values(by=['pdb_code'])
    gpcrdb_info["ligands"] = gpcrdb_info["ligands"].apply(json.dumps)
    print(gpcrdb_info)
    
    # Run the connection engine to the database
    print("- Engine on the connection to MemProtDb database")
    engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')
    
    # Store the info in the database
    print("- Update gpcrdb_info table")
    gpcrdb_info.to_sql(name='gpcrdb_info', if_exists='replace', con=engine)

def update_pdb_info():  
    
    # PFAM GPCRs
    pfamdict = {
    "classA": "PF00001",
    "classB": "PF00002",
    "classC": "PF00003",
    "frizzled": "PF01534",}
    
    # Run the connection engine to the database
    print("- Engine on the connection to MemProtDb database")
    engine = create_engine(f'postgresql+psycopg2://{DB_USER}:{DB_PASSWORD}@localhost/{DB_NAME}')
    
    # Get list of pdbs
    sql_query = f"SELECT pdb_code, preferred_chain, protein FROM gpcrdb_info WHERE resolution < 90 and type='X-ray diffraction';"
    data = pd.read_sql(sql_query, engine)
    
    pdbids = data["pdb_code"].values
    
    # test
    # pdbids = ["8C9W", "8GNG"]
        
    all_info, l_seg_pdb, l_error, l_seg = [], [], [], []
    clas = ""
    j = 1
    for pdb in pdbids:
        print(f"> {j *100 / len(pdbids)} %") 
        print(f"- Getting information for {pdb}")
        url_pdb=f"https://data.rcsb.org/rest/v1/core/entry/{pdb.lower()}"
        response = requests.get(url_pdb)
        data = json.loads(response.text)
        try:
            res = data["rcsb_entry_info"]["diffrn_resolution_high"]["value"]
        except:
            res = ""
        meth = data["rcsb_entry_info"]["experimental_method"]
        ids = data["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
        for id in ids:
            url = f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb.lower()}/{id}'
            response = requests.get(url)
            data = json.loads(response.text)
            #chain = data["entity_poly"]["pdbx_strand_id"]# la chain que correspon al polymer entity  
            chain_pdb = data["rcsb_polymer_entity_container_identifiers"]["asym_ids"]
            chain_auth = data["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]
            for i, ch in enumerate(chain_pdb):
                mapping = f'https://data.rcsb.org/rest/v1/core/polymer_entity_instance/{pdb.lower()}/{ch}' 
                response = requests.get(mapping)
                data_map = json.loads(response.text)
                auth_map = data_map["rcsb_polymer_entity_instance_container_identifiers"]["auth_to_entity_poly_seq_mapping"]# ["852","853","854","855","856","
                try:
                    l_rcsb = data["rcsb_polymer_entity_align"]
                    for rcsb in l_rcsb: 
                        # uniprot
                        acc=rcsb['reference_database_accession']
                        url_uni = f'https://rest.uniprot.org/uniprotkb/{acc}.txt'
                        response = requests.get(url_uni)
                        data_uni = response.text
                        # Get the Uniprot entry
                        search = re.search('ID   ()\w+', data_uni)    
                        entry = str(search.group()).replace("ID", "").replace(" ", "").lower()
                        # extreiem el pfam 
                        search = re.search("Pfam; \w+", data_uni)    
                        pfam = str(search.group()).replace("Pfam;", "").replace(" ", "")
                        if pfam in pfamdict.values():
                            clas = "receptor"
                        else:
                            clas = "fusion protein"
                        #Segment region
                        l_rscb_fea = data["rcsb_polymer_entity_feature"] 
                        for rscb_fea in l_rscb_fea:
                            # el receptor comença un segment a ref_beg_seq_id 
                            l_reg = rcsb['aligned_regions']
                            for reg in l_reg:
                                res_ini_pdb = reg['entity_beg_seq_id'] - 1  # 1 --: 0 
                                res_end_pdb = res_ini_pdb + reg['length'] - 1 
                                res_ini_auth = auth_map[res_ini_pdb] 
                                res_end_auth = auth_map[res_end_pdb] 
                                l_seg.append(f"{res_ini_pdb}-{res_end_pdb}")
                                l_seg_pdb.append(f"{res_ini_auth}-{res_end_auth}")
                            if not [pdb,res,meth,id,acc,entry,pfam,chain_auth[i][0],",".join(l_seg_pdb), clas] in all_info:
                                if clas != "":
                                    all_info.append([pdb,res,meth,id,acc,entry,pfam,chain_auth[i][0],",".join(l_seg_pdb), clas])
                                elif clas == "" and pfam == "":
                                    all_info.append([pdb,res,meth,id,acc,entry,pfam,chain_auth[i][0],",".join(l_seg_pdb), clas])
                            l_info, l_seg_pdb, l_seg = [], [], []
                            clas = ""
                except Exception as e:
                    if not [pdb,res,meth,id,"","","",chain_auth[i][0],"","", ""] in all_info:
                        all_info.append([pdb,res,meth,id,"","","",chain_auth[i][0],"","", ""])
                    l_info, l_seg_pdb, l_seg = [], [], []
                    clas = ""
                    if pdb not in l_error:
                        l_error.append(pdb) 
        j += 1
        
    print(l_error)
    gpcr_pdb = pd.DataFrame(all_info, columns=['pdb', 'resolution', 'method', 'id', 'uniprot_accesion', 'uniprot_entry', 'pfam', 'chain', 'segments', 'type', 'what'])
    gpcr_pdb = gpcr_pdb.drop(['what'],axis = 1)
    print(gpcr_pdb)
    
    # Store the info in the database
    print("- Update gpcrdb_info table")
    gpcr_pdb.to_sql(name='gpcr_pdb', if_exists='replace', con=engine)
    
# MAIN
update_gpcrdb_info()
update_pdb_info()
