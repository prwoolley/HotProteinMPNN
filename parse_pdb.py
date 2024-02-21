import torch
import numpy as np
import Bio.PDB as bp
import os

# Location of Alphafold2 predicted structures
#pdbs = '/stor/work/Ellington/ProteinMPNN/HotProtein/s2c5/d2_clean/'
pdbs = '/stor/work/Ellington/ProteinMPNN/HotProtein/S/AF2/PDB/'


RES_NAMES = [
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL'
]

RES_NAMES_1 = 'ARNDCQEGHILKMFPSTWYV'

to1letter = {aaa:a for a,aaa in zip(RES_NAMES_1,RES_NAMES)}
to3letter = {a:aaa for a,aaa in zip(RES_NAMES_1,RES_NAMES)}

ATOM_NAMES = [
    ("N", "CA", "C", "O", "CB"), # ala
    ("N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"), # arg
    ("N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"), # asn
    ("N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"), # asp
    ("N", "CA", "C", "O", "CB", "SG"), # cys
    ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"), # gln
    ("N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"), # glu
    ("N", "CA", "C", "O"), # gly
    ("N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"), # his
    ("N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"), # ile
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"), # leu
    ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"), # lys
    ("N", "CA", "C", "O", "CB", "CG", "SD", "CE"), # met
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"), # phe
    ("N", "CA", "C", "O", "CB", "CG", "CD"), # pro
    ("N", "CA", "C", "O", "CB", "OG"), # ser
    ("N", "CA", "C", "O", "CB", "OG1", "CG2"), # thr
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE2", "CE3", "NE1", "CZ2", "CZ3", "CH2"), # trp
    ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"), # tyr
    ("N", "CA", "C", "O", "CB", "CG1", "CG2") # val
]
        
idx2ra = {(RES_NAMES_1[i],j):(RES_NAMES[i],a) for i in range(20) for j,a in enumerate(ATOM_NAMES[i])}

aa2idx = {(r,a):i for r,atoms in zip(RES_NAMES,ATOM_NAMES) 
          for i,a in enumerate(atoms)}
aa2idx.update({(r,'OXT'):3 for r in RES_NAMES})

ATOM_NAMES = ["N", "CA", "C", "O", "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE", "CE1", "CE2"]

aa = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}


def parse_af2_pdb(fname, pdb_id):
    parser = bp.PDBParser()
    structure = parser.get_structure('X', fname)
    
    chain = list(structure.get_chains())[0]
    residues = list(chain.get_residues())

    seq = "".join([aa[r.resname] for r in residues])
    
    n_res = len(residues)
    xyz = torch.full((n_res, len(ATOM_NAMES), 3), np.nan)
    #mask = torch.zeros((n_res, len(ATOM_NAMES)), dtype=torch.long)
    n_res = len(residues)
    #n_atom = max([len(list(r.get_atoms())) for r in residues])
    mask = torch.zeros(n_res, len(ATOM_NAMES))
    bfac = torch.full((n_res, len(ATOM_NAMES)), np.nan) 
    occ = torch.zeros(n_res, len(ATOM_NAMES))
    
    #atoms_per_res = [len(list(r.get_atoms())) for r in residues]
    #start = 0
    #for i, n_atoms in enumerate(atoms_per_res):
    #    end = start + n_atoms
    #    mask[i, :n_atoms] = 1
    #    start = end
        

    # Fill in values if atom exists
    for i, r in enumerate(residues):
        for j, atom in enumerate(r.get_atoms()):
            if atom.name in ATOM_NAMES:
                k = ATOM_NAMES.index(atom.name)
                xyz[i,k,:] = torch.tensor(atom.coord)
                mask[i,k] = 1
                bfac[i,k] = atom.bfactor
                occ[i,k] = atom.occupancy
            
    data = {'seq': seq, 
            'xyz': xyz,
            'mask': mask,
            'bfac': bfac,
            'occ': occ}

    metadata = {'id': pdb_id,
                'chains': ['A'],
                'seq': [[seq]]}
    
    return data, metadata

def prep_files(items):
    for pdb in items:
        fname = pdbs+pdb
        if os.path.exists(fname) and os.path.isfile(fname):
            pdb_name = pdb.split('-')[1]
            data, metadata = parse_af2_pdb(fname,pdb_name)
            torch.save(data,'/stor/work/Ellington/ProteinMPNN/HotProtein/S/PDB_pt/'+pdb_name+'_A.pt')
            torch.save(data,'/stor/work/Ellington/ProteinMPNN/HotProtein/Models/training_data/pdb/'+pdb_name+'_A.pt')
            torch.save(metadata,'/stor/work/Ellington/ProteinMPNN/HotProtein/S/PDB_pt/'+pdb_name+'.pt')
            torch.save(metadata,'/stor/work/Ellington/ProteinMPNN/HotProtein/Models/training_data/pdb/'+pdb_name+'.pt')
            #with open('/stor/work/Ellington/ProteinMPNN/HotProtein/s2c2_prepped/unrelaxed_model_1_ptm/s2c2.unrelaxed_model_1_ptm.fa', 'a') as file:
            #    file.write('>'+pdb+'\n'+data['seq']+'\n')
                

items = os.listdir(pdbs)
prep_files(items)


