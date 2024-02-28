import argparse
import glob
import csv
import numpy as np
import torch
from dateutil import parser 
from utils_mmCIF import worker_init_fn, get_pdbs, loader_pdb, PDB_dataset, StructureDataset, StructureLoader
from model_utils_mmCIF import featurize, loss_nll, ProteinMPNN

# This function is pulled from utils_Baker.py but modified for testing
def build_training_clusters(params):
    test_ids = set([int(l) for l in open(params['TEST']).readlines()])

    # read & clean list.csv
    with open(params['LIST'], 'r') as f:
        reader = csv.reader(f)
        next(reader)
        rows = [[r[0],r[3],int(r[4])] for r in reader
                if float(r[2])<=params['RESCUT'] and
                parser.parse(r[1])<=parser.parse(params['DATCUT'])]
    
    ## we want to remove proteins that are listed in list.csv that aren't actually in the dataset
    files = glob.glob('../../training/ProteinMPNN_training_data//pdb_2021aug02/pdb/*/*.pt')
    files = [x.split('/')[-1].split('.pt')[0] for x in files]
    files = set(files) # set for faster membership checking

    # compile test sets
    test = {}
    for r in rows:
        if r[0] in files:
            if r[2] in test_ids:
                if r[2] in test.keys():
                    test[r[2]].append(r[:2])
                else:
                    test[r[2]] = [r[:2]]
     
    return test


def main(args):
    classes = ['class_0','class_1','class_2','class_3','class_4','class_0_1_2','class_3_4']
    

    device = torch.device("cuda:1" if (torch.cuda.is_available()) else "cpu")
    print(device)

    data_path = '../../training/ProteinMPNN_training_data/pdb_2021aug02'
    params = {
        "LIST"    : f"{data_path}/list.csv",
        "TEST"    : f"{data_path}/test_clusters.txt",
        "DIR"     : f"{data_path}",
        "DATCUT"  : "2030-Jan-01",
        "RESCUT"  : 3.5, #resolution cutoff for PDBs
        "HOMO"    : 0.70 #min seq.id. to detect homo chains
        }

    LOAD_PARAM = {'batch_size': 1,
                    'shuffle': True,
                    'pin_memory':False,
                    'num_workers': 4,
                    }

    # loading testing data
    test = build_training_clusters(params) 
    test_set = PDB_dataset(list(test.keys()), loader_pdb, test, params)
    test_loader = torch.utils.data.DataLoader(test_set, worker_init_fn=worker_init_fn, **LOAD_PARAM)
    pdb_dict_test = get_pdbs(test_loader, 1, args.max_protein_length, args.num_examples_per_epoch)
    dataset_test = StructureDataset(pdb_dict_test, truncate=None, max_length=args.max_protein_length)
    loader_test = StructureLoader(dataset_test, batch_size=args.batch_size)  

    # Once the test data is loaded, we will iterate through the models
    model = ProteinMPNN(node_features=args.hidden_dim, 
                edge_features=args.hidden_dim, 
                hidden_dim=args.hidden_dim, 
                num_encoder_layers=args.num_encoder_layers, 
                num_decoder_layers=args.num_encoder_layers, 
                k_neighbors=args.num_neighbors, 
                dropout=args.dropout, 
                augment_eps=args.backbone_noise)
    model.to(device)

    #Iterating through each model
    for i in classes:
        model_list = sorted(glob.glob('../Models/'+i+'/model/model_weights/*'))
        logfile = '../Models/'+i+'/model/BakerLab_test_log.txt'
        for j in range(0, len(model_list)):
            PATH = model_list[j]
            print(model_list[j])
            
            checkpoint = torch.load(PATH)
            epoch = checkpoint['epoch'] #write epoch from the checkpoint
            model.load_state_dict(checkpoint['model_state_dict'])
            print("Loaded model from epoch: ", epoch)

            model.eval()
            with torch.no_grad():
                test_sum, test_weights = 0., 0.
                test_acc = 0.
                for _, batch in enumerate(loader_test):
                    X, S, mask, lengths, chain_M, residue_idx, mask_self, chain_encoding_all = featurize(batch, device)
                    log_probs = model(X, S, mask, chain_M, residue_idx, chain_encoding_all)
                    mask_for_loss = mask * chain_M
                    loss, loss_av, true_false = loss_nll(S, log_probs, mask_for_loss)

                    test_sum += torch.sum(loss * mask_for_loss).cpu().data.numpy()
                    test_acc += torch.sum(true_false * mask_for_loss).cpu().data.numpy()
                    test_weights += torch.sum(mask_for_loss).cpu().data.numpy()

            test_loss = test_sum / test_weights
            test_accuracy = test_acc / test_weights
            test_perplexity = np.exp(test_loss)

            test_perplexity_ = np.format_float_positional(np.float32(test_perplexity), unique=False, precision=3)
            test_accuracy_ = np.format_float_positional(np.float32(test_accuracy), unique=False, precision=3)
            with open(logfile, 'a') as f:
                    f.write(f'epoch: {epoch}, test: {test_perplexity_}, test_acc: {test_accuracy_}\n')
            print(f'epoch: {epoch}, test: {test_perplexity_}, test_acc: {test_accuracy_}')

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument("--path_for_test_data", type=str, default="my_path/pdb_test", help="path for loading test data")
    argparser.add_argument("--path_for_outputs", type=str, default="./exp_020", help="path for logs and model weights")
    argparser.add_argument("--batch_size", type=int, default=10000, help="number of tokens for one batch")
    argparser.add_argument("--max_protein_length", type=int, default=10000, help="maximum length of the protein complex")
    argparser.add_argument("--hidden_dim", type=int, default=128, help="hidden model dimension")
    argparser.add_argument("--num_encoder_layers", type=int, default=3, help="number of encoder layers")
    argparser.add_argument("--num_decoder_layers", type=int, default=3, help="number of decoder layers")
    argparser.add_argument("--num_neighbors", type=int, default=48, help="number of neighbors for the sparse graph")
    argparser.add_argument("--dropout", type=float, default=0.1, help="dropout level; 0.0 means no dropout")
    argparser.add_argument("--rescut", type=float, default=3.5, help="PDB resolution cutoff")
    argparser.add_argument("--debug", type=bool, default=False, help="minimal data loading for debugging")
    argparser.add_argument("--backbone_noise", type=float, default=0.2, help="amount of noise added to backbone during training")
    argparser.add_argument("--num_examples_per_epoch", type=int, default=1000000, help="number of training example to load for one epoch")
    args = argparser.parse_args()
    main(args)