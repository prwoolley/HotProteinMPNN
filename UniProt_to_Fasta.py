import requests
import pandas as pd

mpnn_uniprot = pd.read_csv('/stor/work/Ellington/ProteinMPNN/HotProtein/ProteinMPNN_TrainingData_Unique_UniProtAccessions.txt', header=None)
mpnn_fasta_directory = '/stor/work/Ellington/ProteinMPNN/HotProtein/ProteinMPNN_TrainingData_Fasta'

# Base URL for UniProt REST API
base_url = "https://rest.uniprot.org/uniprotkb/"

# Retrieve and save Fasta files
for accession in mpnn_uniprot[0]:
    url = f"{base_url}{accession}.fasta"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for error status codes
        with open(f"{mpnn_fasta_directory}/{accession}.fasta", "w") as f:
            f.write(response.text)
        print(f"Saved Fasta file for {accession}")
    except requests.exceptions.RequestException as e:
        print(f"Error retrieving {accession}: {e}")
    except Exception as e:
        print(f"Error retrieving {accession}: {e}")