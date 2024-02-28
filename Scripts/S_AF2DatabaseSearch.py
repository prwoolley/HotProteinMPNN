import ijson
import requests, sys, json
import ipywidgets as wgt
from urllib.request import urlopen
import os
import pandas as pd
import gemmi

def Search3DBeacons(ID):
  WEBSITE_API = "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/"

  r = ijson.parse(urlopen(f"{WEBSITE_API}{ID}.json"))
  structures = list(ijson.items(r, "structures.item", use_float=True))
  return structures


# Function to download a file from a given URL and save it to the Google Drive
def download_file(url):
  os.chdir('/stor/work/Ellington/ProteinMPNN/HotProtein/S/AF2/CIF/')
  subprocess.run(["wget", url])  # Use a list of arguments for clarity

def convert_to_pdb(url):
    name = url.split('/')[-1]
    try:
        cif = gemmi.read_structure('/stor/work/Ellington/ProteinMPNN/HotProtein/S/AF2/CIF/'+name)
        pdb_str = cif.make_minimal_pdb()
        with open('/stor/work/Ellington/ProteinMPNN/HotProtein/S/AF2/PDB/'+name[:-4]+'.pdb', 'w') as f:
            f.write(pdb_str)
    except:
        # Add name to text file
        with open('/stor/work/Ellington/ProteinMPNN/HotProtein/S/AF2/failed.txt', 'a') as f:
            f.write(name+'\n')
    
    

uniprot_ids = os.listdir('/stor/work/Ellington/ProteinMPNN/HotProtein/S/S_target/')
uniprot_ids = [x.split('.')[0] for x in uniprot_ids]

af_results = pd.DataFrame(columns=['uniprot_id','af_model_url'])

for i in range(len(uniprot_ids)):
    try:
        search_result = Search3DBeacons(uniprot_ids[i])
    except requests.HTTPError as e:
        af_results.loc[len(af_results.index)] = [uniprot_ids[i],None] 
        continue
    except Exception as e:
        af_results.loc[len(af_results.index)] = [uniprot_ids[i],None] 
        continue

    for j in search_result:
        if j['summary']['provider'] == 'AlphaFold DB':
            af_results.loc[len(af_results.index)] = [uniprot_ids[i],j['summary']['model_url']] 
    
    if i%100 == 0:
        af_results.to_csv('/stor/work/Ellington/ProteinMPNN/HotProtein/S/S_AF2_model_urls.csv',sep=',',index=False)   
af_results.to_csv('/stor/work/Ellington/ProteinMPNN/HotProtein/S/S_AF2_model_urls.csv',sep=',',index=False)    

af_results = pd.read_csv('/stor/work/Ellington/ProteinMPNN/HotProtein/S/S_AF2_model_urls.csv')
af_results = af_results.dropna()
af_results['af_model_url'].apply(lambda url: download_file(url) if url else None)
af_results.dropna()['af_model_url'].apply(lambda url: convert_to_pdb(url))
