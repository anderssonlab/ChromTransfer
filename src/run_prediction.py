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

print('cell line promoter')

a549_prom = prediction_cell_line1(one_hot_seq, 'A549', 'cage')
hct116_prom = prediction_cell_line1(one_hot_seq, 'HCT116', 'cage')
hepg2_prom = prediction_cell_line1(one_hot_seq, 'HEPG2', 'cage')
k562_prom = prediction_cell_line1(one_hot_seq, 'K562', 'cage')

print('cell line enhancer')

a549_enha = prediction_cell_line1(one_hot_seq, 'A549', 'starr')
hct116_enha = prediction_cell_line1(one_hot_seq, 'HCT116', 'starr')
hepg2_enha = prediction_cell_line1(one_hot_seq, 'HEPG2', 'starr')
k562_enha = prediction_cell_line1(one_hot_seq, 'K562', 'starr')



open_chromatin = np.concatenate(open_chromatin, axis=0)

a549_chromatin = np.concatenate(a549_chromatin, axis=0)
a549_prom = np.concatenate(a549_prom, axis=0)
a549_enha = np.concatenate(a549_enha, axis=0)

hct116_chromatin = np.concatenate(hct116_chromatin, axis=0)
hct116_prom = np.concatenate(hct116_prom, axis=0)
hct116_enha = np.concatenate(hct116_enha, axis=0)

hepg2_chromatin = np.concatenate(hepg2_chromatin, axis=0)
hepg2_prom = np.concatenate(hepg2_prom, axis=0)
hepg2_enha = np.concatenate(hepg2_enha, axis=0)

k562_chromatin = np.concatenate(k562_chromatin, axis=0)
k562_prom = np.concatenate(k562_prom, axis=0)
k562_enha = np.concatenate(k562_enha, axis=0)


print('creating final output file')

final_pred_array = np.array([open_chromatin, 
                             a549_chromatin, a549_prom, a549_enha,
                             hct116_chromatin, hct116_prom, hct116_enha,
                             hepg2_chromatin, hepg2_prom, hepg2_enha,
                             k562_chromatin, k562_prom, k562_enha])

tmp_df = pd.DataFrame(final_pred_array, columns = name_seq, index = ['TL_Open_Closed', 
                                                                     'A549_rDHS', 'A549_prom','A549_enh',
                                                                     'HCT116_rDHS', 'HCT116_prom','HCT116_enh',
                                                                     'HEPG2_rDHS', 'HEPG2_prom','HEPG2_enh',
                                                                     'K562_rDHS', 'K562_prom','K562_enh']).T

final_df = tmp_df[['TL_Open_Closed', 
                   'A549_rDHS', 'HCT116_rDHS', 'HEPG2_rDHS', 'K562_rDHS',
                   'A549_prom', 'HCT116_prom', 'HEPG2_prom', 'K562_prom',
                   'A549_enh', 'HCT116_enh', 'HEPG2_enh', 'K562_enh']].round(3)

final_df.to_csv(args.output+'output.csv')
os.remove(args.fasta)


