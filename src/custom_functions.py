import collections
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
import os
import seaborn as sns
import sys



def extract_fasta(fasta_file, bed_file, csv_file):
    bed_fasta = 'bedtools getfasta -fi ' + fasta_file + ' -bed ' + bed_file + ' -fo ' + csv_file  +' -tab -fullHeader'
    os.system(bed_fasta)
    #print(bed_fasta)
    
    
def final_dataframe(csv_file):
    fasta_df = pd.read_csv(csv_file, sep='\t')
    fasta_df.columns = ['Region', 'sequence']
    #df_final_for_npz = dataframe_info.merge(fasta_df, on='Region', how='inner')
    return fasta_df


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


def parse_sequence(file_info):
    
    oneHot = np.identity(4)
    
    seqs = []
    df = file_info 
    ids_array = np.array((df.Region))
    df.sequence = df.sequence.str.upper()
    
    for nuc_seq in df.sequence:
        seqs.append(encode_sequence(nuc_seq, oneHot))
    seq_array = np.array(seqs)
    #print(ids_array, seq_array[0])
    return seq_array, ids_array


