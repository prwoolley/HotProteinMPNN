import argparse
import requests
import pandas as pd


def main(args):
    with open(args.uniprot_txt, 'r') as file:
        uniprots = [line.strip() for line in file]    

    # Base URL for UniProt REST API
    base_url = "https://rest.uniprot.org/uniprotkb/"

    # Retrieve and save Fasta files
    for accession in uniprots:
        url = f"{base_url}{accession}.fasta"
        try:
            response = requests.get(url)
            response.raise_for_status()  # Raise an exception for error status codes
            with open(f"{args.fasta_output_dir}/{accession}.fasta", "w") as f:
                f.write(response.text)
            print(f"Saved Fasta file for {accession}")
        except requests.exceptions.RequestException as e:
            print(f"Error retrieving {accession}: {e}")
        except Exception as e:
            print(f"Error retrieving {accession}: {e}")

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument("--uniprot_txt", type=str, default="my_path/UniProt.txt", help="path for text file with UniProt IDs")
    argparser.add_argument("--fasta_output_dir", type=str, default="my_path/FASTA", help="path for where to dump FASTA files")
    
    args = argparser.parse_args()
    main(args)