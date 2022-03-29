# %%
import argparse

import tensorflow as tf
import numpy as np

import igrads

from rbpnet.attribution import attribution

# %%
base2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

# %%
def sequence2int(sequence):
    return [base2int.get(base, 999) for base in sequence]

# %%
def sequence2onehot(sequence):
    return tf.one_hot(sequence2int(sequence), depth=4)

# %%
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('afasta')
    parser.add_argument('-o', '--output', default='{name}.png')
    args = parser.parse_args()

    with open(args.afasta) as f:
        for line in f:
            if line[0] != '>':
                raise ValueError()
            name = line.strip()[1:].split(' ')[0]
            print(name)

            inputs = sequence2onehot(f.readline().strip().upper()).numpy()
            attribution = np.array(list(map(float, f.readline().strip().split(' '))))
            attribution = tf.expand_dims(attribution, axis=1).numpy()

            # print(inputs.shape)
            # print(attribution.shape)

            assert inputs.shape[0] == attribution.shape[0]

            logo = igrads.plot_sequence_attribution(inputs * attribution, width=40, height=4)
            logo.fig.savefig(args.output.format(name=name))

# %%
if __name__ == '__main__':
    main()