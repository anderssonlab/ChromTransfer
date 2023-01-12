# %%
import argparse
import tensorflow as tf
from traitlets import default

# %%
base2int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

def sequence2int(sequence, mapping=base2int):
    return [mapping.get(base, 999) for base in sequence]

def sequence2onehot(sequence, mapping=base2int):
    return tf.one_hot(sequence2int(sequence, mapping), depth=4)

def load_fasta(fasta, positive_name):
    with open(fasta) as f:
        for line in f:
            assert line[0] == '>'
            header, sequence = line.strip(), f.readline().strip()
            label = int(header[1:].split(':')[0] == positive_name)
            yield tf.cast(sequence2onehot(sequence), tf.float32), tf.cast(label, tf.float32)

def load_dataset(fasta, positive_name, batch=128, cache=False, shuffle=False):
    dataset = tf.data.Dataset.from_generator(lambda: load_fasta(fasta, positive_name), output_types=(tf.float32, tf.float32))
    if cache:
        dataset = dataset.cache()
    if shuffle:
        dataset = dataset.shuffle(shuffle)
    dataset = dataset.batch(batch)
    return dataset

# %%
def finetune(dataset, pretrained_model_filepath, epochs):
    # models/tensorflow1/tl_model_tf1/model_compiled.h5

    pretrained_model = tf.keras.models.load_model(pretrained_model_filepath)
    
    # frozen body model
    body = tf.keras.Model(inputs=pretrained_model.inputs, outputs=pretrained_model.layers[-8].output)
    body.trainable = False
    
    # trainable model head
    head = tf.keras.layers.Dense(1024, name='head_dense_1', activation='relu')(body.outputs[0])
    head = tf.keras.layers.Dropout(0.1, name='head_dropout_1')(head)
    head = tf.keras.layers.Dense(32, name='head_dense_2', activation='relu')(head)
    head = tf.keras.layers.Dropout(0.1, name='head_dropout_2')(head)
    head = tf.keras.layers.Dense(1, name='head_dense_3', activation='sigmoid')(head)
    
    model = tf.keras.Model(inputs=body.inputs, outputs=head)
    
    model.compile(loss='binary_crossentropy', optimizer=tf.keras.optimizers.Adam(learning_rate=0.001))
    model.fit(dataset, epochs=epochs)
    
    return model

# %%
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', metavar='<sequences.fasta>', help='Training sequences in FASTA format.')
    parser.add_argument('--positive-sequence-name', required=True, help='Name of positive-class samples. E.g. \'A549\' for headers of form \'>A549::chr2:1527286-1527886\'. Other names are treated as negatve instances.')
    parser.add_argument('--pretrained-model', metavar='<pretrained-model.h5>', required=True, help='Pretrained model.')
    parser.add_argument('--epochs', default=20, type=int, help='Number of finetuning epochs.')
    parser.add_argument('--shuffle', default=100_000, type=int, help='Dataset shuffle buffer size.')
    parser.add_argument('--cache', action='store_true', default=False, help='Cache dataset in memory.')
    parser.add_argument('--output', metavar='<model.h5>', required=True, help='Finetuned model destination.')
    args = parser.parse_args()
    
    # load_data
    dataset = load_dataset(args.fasta, positive_name=args.positive_sequence_name, cache=args.cache, shuffle=args.shuffle)
    
    model = finetune(dataset, args.pretrained_model, args.epochs)
    model.save(args.output)
    
# %%
if __name__ == '__main__':
    main()