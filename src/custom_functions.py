from Bio import SeqIO
import collections
import numpy as np
import os
import pandas as pd
import pybedtools
from pybedtools import BedTool
import seaborn as sns
import sys
import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter(action='ignore', category=FutureWarning)
import tensorflow 
import keras
from keras.models import model_from_json, load_model




def encode_sequence(sequence, onehot_array):
    code = list("ACGT")
    seq_len = len(sequence)
    n_feat = len(code)
    arr_2d = np.zeros((seq_len, n_feat))
    for i, nucleotide in enumerate(sequence):
        if nucleotide not in code:
            arr_2d[i] = np.zeros(n_feat)
        else:
            ind = code.index(nucleotide)

        arr_2d[i] = onehot_array[ind]

    return arr_2d


def parse_sequence(bed_file, genome_file, fasta_file):
    
    bed_info = pybedtools.BedTool(bed_file)
    bed_info = bed_info.sequence(fi=genome_file).save_seqs(fn=fasta_file)
    
    oneHot = np.identity(4)
    
    seqs = []
    ids = []
    
    for info_ in SeqIO.parse(fasta_file, 'fasta'):
        ids.append(info_.id)
        seqs.append(encode_sequence(info_.seq.upper(), oneHot))
            
    seq_array = np.array(seqs)
    ids_array = np.array(ids)
    #print(ids_array, seq_array)
    return seq_array, ids_array



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

    partitions_interval = np.arange(3)
    inner_partitions_interval = partitions_interval[partitions_interval]

    for val_partition in inner_partitions_interval:
        print(val_partition)
        
        model_chromatin = 'models/tensorflow1/rDHS/'+ cell_line +'/partition_%i/model_compiled.h5' % (val_partition)
        model_promoter = 'models/tensorflow1/cage/'+ cell_line +'/partition_%i/model_compiled.h5' % (val_partition)
        model_enhancer = 'models/tensorflow1/starr/'+ cell_line +'/partition_%i/model_compiled.h5' % (val_partition)

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
