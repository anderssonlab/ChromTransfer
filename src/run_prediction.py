import sys
import os
import argparse
from custom_functions import *



parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bed', type=str,
    help="input file in bed format")
parser.add_argument('-f', '--fasta', type=str,
    help="save file in fasta format")
parser.add_argument('-r', '--genome', type=str,
    help="genome in fasta format")
parser.add_argument('-c', '--csv', type=str,
    help="genome in fasta format")
parser.add_argument('-g','--gpu', default=-1, type=int,
    help='the GPU number, -1 indicates CPU')
parser.add_argument('-o','--output', default='output_folder', type=str,
    help='prefix of the output folder where the outputs are saved')
args = parser.parse_args()

gpu = args.gpu

if gpu > -1:
	device = '/gpu:%i' % gpu
	os.environ['CUDA_VISIBLE_DEVICES']=str(args.gpu)
else:
	device = '/cpu:0'

#os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

if args.bed == None:
    parser.print_help()
    sys.stderr.write("Please specify the bed file!\n")
    sys.exit(1)


## from bed to one hot encoding:


one_hot_seq, name_seq = parse_sequence(args.bed, args.genome, args.fasta)

    
all_preds = []
all_y = []
all_l = []
all_id = []

#X = np.asarray(X).astype(np.float32)


open_chromatin = prediction_open_chromatin(one_hot_seq)
a549_chromatin, a549_promoter, a549_enhancer = prediction_cell_line(one_hot_seq, 'A549')
hct116_chromatin, hct116_promoter, hct116_enhancer = prediction_cell_line(one_hot_seq, 'HCT116')
hepg2_chromatin, hepg2_promoter, hepg2_enhancer = prediction_cell_line(one_hot_seq, 'HEPG2')
k562_chromatin, k562_promoter, k562_enhancer = prediction_cell_line(one_hot_seq, 'K562')


print(open_chromatin, a549_chromatin, hct116_chromatin, hepg2_chromatin, k562_chromatin,
      a549_promoter, hct116_promoter, hepg2_promoter, k562_promoter,
      a549_enhancer, hct116_enhancer, hepg2_enhancer, k562_enhancer)
'''
all_preds.append(test_type_avg)
all_y.append(y_val)
all_l.append(id_val)

#prediction
all_preds_conc = np.concatenate(all_preds, axis=0)
#label
all_y_conc = np.concatenate(all_y,axis=0)
#score of the label predicted
all_l_conc = np.concatenate(all_l,axis=0)

for i, c, p  in zip(all_y_conc, all_preds_conc, all_l_conc):
    print(str(i)+','+str(c[0])+','+str(p))#+','+str(c[1])) #+','+str(c[2]))
'''    
    


