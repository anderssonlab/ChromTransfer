# %%
import argparse
from tkinter import E


import tensorflow as tf
import numpy as np
from Bio import SeqIO
from tqdm import tqdm

import igrads

# %%
base2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

# %%
def sequence2int(sequence):
    return [base2int.get(base, 999) for base in sequence]

# %%
def sequence2onehot(sequence):
    return tf.one_hot(sequence2int(sequence), depth=4)

# %%
def count_fasta_seqs(fasta):
    n = 0
    with open(fasta) as f:
        for line in f:
            if line[0] == '>':
                n += 1
    return n

# %%
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta')
    parser.add_argument('-m', '--model')
    parser.add_argument('-o', '--output')
    parser.add_argument('--method', default='grad_x_input', help="'grad_x_input' or 'integrated_gradients'")
    args = parser.parse_args()

    model = tf.keras.models.load_model(args.model)
    n = count_fasta_seqs(args.fasta)

    if args.method == 'grad_x_input':
        afun = igrads.grad_x_input
    elif args.method == 'integrated_gradients':
        afun = igrads.integrated_gradients
    else:
        raise ValueError(f'Unkown method: {args.method}')

    with open(args.output, 'w') as f, tqdm(total=n) as pbar:
        for s in SeqIO.parse(args.fasta, 'fasta'):
            print('>' + s.description, file=f)
            print(s.seq.upper(), file=f)
            inputs = sequence2onehot(s.seq.upper())
            attribution = np.sum(afun(inputs, model).numpy(), axis=1) # np.max makes the attribution non-negative, use else np.sum
            print(' '.join(map(lambda x: f'{float(x):.5f}', attribution)), file=f)
            
            pbar.update(1)

# %%
if __name__ == '__main__':
    main()