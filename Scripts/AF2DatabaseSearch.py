import argparse
import ijson
import requests
import subprocess
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
def download_file(url,cif_output_dir):
  os.chdir(cif_output_dir)
  subprocess.run(["wget", url])  # Use a list of arguments for clarity

def convert_to_pdb(url,cif_output_dir,pdb_output_dir):
    name = url.split('/')[-1]
    try:
        cif = gemmi.read_structure(cif_output_dir+name)
        pdb_str = cif.make_minimal_pdb()
        with open(pdb_output_dir+name[:-4]+'.pdb', 'w') as f:
            f.write(pdb_str)
    except:
        # Add name to text file
        with open('/stor/work/Ellington/ProteinMPNN/HotProtein/S/AF2/failed.txt', 'a') as f:
            f.write(name+'\n')
    
    
def main(args):
    uniprot_ids = os.listdir(args.uniprot_fasta_dir)
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
            af_results.to_csv(args.csv_output,sep=',',index=False)   
    af_results.to_csv(args.csv_output,sep=',',index=False)    

    af_results = pd.read_csv(args.csv_output)
    af_results = af_results.dropna()
    af_results['af_model_url'].apply(lambda url: download_file(url) if url else None)
    af_results.dropna()['af_model_url'].apply(lambda url: convert_to_pdb(url))


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument("--uniprot_fasta_dir", type=str, default="my_path/UniProt", help="path for directory with UniProt fasta files")
    argparser.add_argument("--pdb_output_dir", type=str, default="my_path/AF2/PDB", help="path for where to dump PDB files")
    argparser.add_argument("--cif_output_dir", type=str, default="my_path/AF2/CIF", help="path for where to dump CIF intermediate files")
    argparser.add_argument("--csv_output", type=str, default="my_path/AF2/AF2_model_urls.csv", help="path for where to save CSV with UniProt IDs and AlphaFold2 model URLs")
    
    args = argparser.parse_args()
    main(args)
