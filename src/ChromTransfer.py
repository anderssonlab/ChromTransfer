import sys
import os
import argparse
from custom_functions import *
from keras.preprocessing.sequence import pad_sequences
import numpy as np
import pandas as pd

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




os.makedirs(args.output, exist_ok=True)

## from bed to one hot encoding:

one_hot_seq, name_seq = parse_sequence(args.bed, args.genome, args.fasta)

one_hot_seq = pad_sequences(one_hot_seq, padding="post")


## running the transfer learning method
print('open chromatin prediction')

open_chromatin = prediction_open_chromatin(one_hot_seq)

## running the cell line methods

print('cell line chromatin')

a549_chromatin = prediction_cell_line1(one_hot_seq, 'A549', 'rDHS')
hct116_chromatin = prediction_cell_line1(one_hot_seq, 'HCT116', 'rDHS')
hepg2_chromatin = prediction_cell_line1(one_hot_seq, 'HEPG2', 'rDHS')
k562_chromatin = prediction_cell_line1(one_hot_seq, 'K562', 'rDHS')
gm12878_chromatin = prediction_cell_line1(one_hot_seq, 'GM12878', 'rDHS')
mcf7_chromatin = prediction_cell_line1(one_hot_seq, 'MCF7', 'rDHS')


open_chromatin = np.concatenate(open_chromatin, axis=0)

a549_chromatin = np.concatenate(a549_chromatin, axis=0)

hct116_chromatin = np.concatenate(hct116_chromatin, axis=0)

hepg2_chromatin = np.concatenate(hepg2_chromatin, axis=0)

k562_chromatin = np.concatenate(k562_chromatin, axis=0)

gm12878_chromatin = np.concatenate(gm12878_chromatin, axis=0)

mcf7_chromatin = np.concatenate(mcf7_chromatin, axis=0)


print('creating final output file')

final_pred_array = np.array([open_chromatin,  a549_chromatin, hct116_chromatin, 
                             hepg2_chromatin, k562_chromatin, gm12878_chromatin, mcf7_chromatin])

tmp_df = pd.DataFrame(final_pred_array, columns = name_seq, index = ['TL_Open_Closed', 'A549', 'HCT116',  'HEPG2',  'K562',  'GM12878','MCF7' ]).T

final_df = tmp_df[['TL_Open_Closed', 'A549', 'HCT116', 'HEPG2', 'K562', 'GM12878', 'MCF7',
                   ]].round(3)

final_df.to_csv(args.output+'output.csv')
os.remove(args.fasta)


