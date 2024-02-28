import pandas as pd
import pickle as pkl
import os
from bioservices import UniProt

proteinmpnn_data = pd.read_csv('/stor/work/Ellington/ProteinMPNN/training/ProteinMPNN_training_data/pdb_2021aug02/list.csv')
proteinmpnn_data = proteinmpnn_data['CHAINID'].str.split('_').str[0].unique()
proteinmpnn_data = [x.upper() for x in proteinmpnn_data]
accessions = []

uniprot = UniProt()
for pdb in proteinmpnn_data:
    if len(accessions) % 100 == 0:
        with open('/stor/work/Ellington/ProteinMPNN/HotProtein/ProteinMPNN_TrainingData_UniProtAccessions.txt', 'w') as f:
            for item in accessions:
                f.write("%s\n" % item)
    try:
        results = uniprot.mapping(fr="PDB", to="UniProtKB",query=pdb)  # ACC for UniProtKB AC (identifier)
        if results['results']!=[]:
            accession = results['results'][0]['to']['primaryAccession']
            accessions.append(accession)
        else:
            continue
    except Exception:
        continue
with open('/stor/work/Ellington/ProteinMPNN/HotProtein/ProteinMPNN_TrainingData_UniProtAccessions.txt', 'w') as f:
    for item in accessions:
        f.write("%s\n" % item)