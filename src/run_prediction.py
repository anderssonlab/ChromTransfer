import sys
import os
import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import tensorflow 
import keras
from keras.models import model_from_json, load_model
import argparse
from custom_functions import *



parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bed', type=str,
    help="input file in bed format")
parser.add_argument('-f', '--fasta', type=str,
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

dataframe_info = extract_fasta(args.fasta, args.bed, args.csv)
file_info = final_dataframe(args.csv)
one_hot_seq, name_seq = parse_sequence(file_info)



all_preds = []
all_y = []
all_l = []
all_id = []

#X = np.asarray(X).astype(np.float32)

def prediction_open_chromatin(one_hot_seq):
    n_seq = one_hot_seq.shape[0]
    _chromatin_ =  np.zeros((n_seq, 1))

    #prediction_ = list()

    model_tl = "models/tensorflow1/tl_model_tf1/model_compiled.h5"
    chromatin_tl = load_model(model_tl, compile=False)
    pred_chromatin_tl = chromatin_tl.predict(one_hot_seq, batch_size=1)
    _chromatin_ += pred_chromatin_tl
    
    return _chromatin_
    
    
def prediction_cell_line(one_hot_seq, cell_line):
    n_seq = one_hot_seq.shape[0]
    type_chromatin =  np.zeros((n_seq, 1))
    type_enhancer =  np.zeros((n_seq, 1))
    type_promoter =  np.zeros((n_seq, 1))

    #prediction_ = list()

    partitions_interval = np.arange(3)
    inner_partitions_interval = partitions_interval[partitions_interval]

    for val_partition in inner_partitions_interval:
        #model_file = "%spartition_%i/model_compiled.h5" % (logdir, val_partition)
        model_chromatin = 'models/tensorflow1/rDHS/'+ cell_line +'/partition_%i/model_compiled.h5' % (val_partition)
        model_promoter = 'models/tensorflow1/cage/'+ cell_line +'/partition_%i/model_compiled.h5' % (val_partition)
        model_enhancer = 'models/tensorflow1/starr/'+ cell_line +'/partition_%i/model_compiled.h5' % (val_partition)
        #print(model_file)

        chromatin = load_model(model_chromatin, compile=False)
        promoter = load_model(model_promoter, compile=False)
        enhancer = load_model(model_enhancer, compile=False)
        
        pred_chromatin = chromatin.predict(one_hot_seq, batch_size=1)
        pred_promoter = promoter.predict(one_hot_seq, batch_size=1)
        pred_enhancer = enhancer.predict(one_hot_seq, batch_size=1)
        
        type_chromatin += pred_chromatin
        type_promoter += pred_promoter
        type_enhancer += pred_enhancer


    type_chromatin_avg = type_chromatin/3.0
    type_promoter_avg = type_promoter/3.0
    type_enhancer_avg = type_enhancer/3.0


    return type_chromatin_avg, type_promoter_avg, type_enhancer_avg

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
    


