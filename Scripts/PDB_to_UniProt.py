import argparse
from bioservices import UniProt



def main(args):
    with open(args.pdb_txt, 'r') as file:
        pdbs = [line.strip() for line in file]
    pdbs = [x.upper() for x in pdbs]
    accessions = []
    uniprot = UniProt()
    for pdb in pdbs:
        if len(accessions) % 100 == 0: # Save checkpoint every 100 accessions
            with open(args.csv_output, 'w') as f:
                for item in accessions:
                    f.write("%s\n" % item)
        try:
            results = uniprot.mapping(fr="PDB", to="UniProtKB",query=pdb)  # ACC for UniProtKB AC (identifier)
            if results['results'] != []:
                accession = results['results'][0]['to']['primaryAccession']
                accessions.append(f"{accession},{pdb}")
            else:
                continue
        except Exception:
            continue
    with open(args.csv_output, 'w') as f:
        for item in accessions:
            f.write("%s\n" % item)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument("--pdb_txt", type=str, default="my_path/pdbs.txt", help="path for text file with PDB IDs")
    argparser.add_argument("--csv_output", type=str, default="my_path/pdb_uniprot.csv", help="path for outputting PDB IDs and mapped UniProt accessions")
    
    args = argparser.parse_args()
    main(args)
