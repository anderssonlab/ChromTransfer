from Bio import SeqIO
import collections
import numpy as np
import os
import pandas as pd
import pybedtools
from pybedtools import BedTool
import sys
import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter(action='ignore', category=FutureWarning)
import tensorflow 
import keras
from tensorflow.keras.models import model_from_json, load_model




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
        info_seq = info_.seq.upper().replace('N','')
        ids.append(info_.id)
        seqs.append(encode_sequence(info_seq), oneHot))
        
    seq_array = np.array(seqs)
    ids_array = np.array(ids)
    return seq_array, ids_array



def prediction_open_chromatin(one_hot_seq):
    
    n_seq = one_hot_seq.shape[0]
    _chromatin_ =  np.zeros((n_seq, 1))
    model_tl = "models/tensorflow1/tl_model_tf1/model_compiled.h5"
    chromatin_tl = load_model(model_tl, compile=False)
    pred_chromatin_tl = chromatin_tl.predict(one_hot_seq, batch_size=16)
    _chromatin_ += pred_chromatin_tl
    return _chromatin_
    
    
def prediction_cell_line1(one_hot_seq, cell_line, tech_):
    
    n_seq = one_hot_seq.shape[0]
    type_cell_ =  np.zeros((n_seq, 1))
    model_cell = 'models/tensorflow1/'+tech_+'/'+ cell_line +'/model_compiled.h5' 
    cell_ = load_model(model_cell, compile=False)
    pred_cell = cell_.predict(one_hot_seq, batch_size=16)
    type_cell_ += pred_cell
    return type_cell_
        
