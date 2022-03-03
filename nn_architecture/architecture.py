#!/usr/bin/env python
from tensorflow.keras import backend as K
import tensorflow as tf
from keras.utils import np_utils
from keras.layers import MaxPooling1D, GlobalAveragePooling1D, GlobalMaxPooling1D, Conv1D
from keras.models import Model
from keras.layers.merge import concatenate, add
from keras.layers import Dense, Input, Dropout, Activation
import keras
from keras.layers.normalization import BatchNormalization
from keras.regularizers import l2
from keras.models import load_model

inizialiser = "he_normal"
regularization_v = 1e-4


def conv1d_module(input_layer, window_size, filter_size, regularization_v):
	layer1d = Conv1D(filter_size, window_size, padding='same', kernel_initializer=inizialiser,
	               kernel_regularizer=l2(regularization_v), activation='relu')(input_layer)
	layer1d = BatchNormalization()(layer1d)
	layer1d = Dropout(0.1)(layer1d)

	return layer1d


def shortcut_module(input_layer, window_size, filter_size, regularization_v):
	block_1 = Conv1D(filter_size, window_size, padding='same',kernel_initializer=inizialiser,
					 kernel_regularizer=l2(regularization_v), activation='relu')(input_layer)
	block_1 = BatchNormalization()(block_1)
	block_1 = Dropout(0.1)(block_1)

	block_1 = Conv1D(filter_size, window_size, padding='same',kernel_initializer=inizialiser,
					 kernel_regularizer=l2(regularization_v), activation='relu')(block_1)
	block_1 = BatchNormalization()(block_1)
	block_1 = Dropout(0.1)(block_1)

	shortcut = Conv1D(filter_size, window_size, padding='same', kernel_initializer=inizialiser,
					  kernel_regularizer=l2(regularization_v), activation='relu')(input_layer)
	shortcut = BatchNormalization()(shortcut)
	shortcut = Dropout(0.1)(shortcut)
	connection = add([block_1, shortcut])

	return connection


def no_shortcut_module(input_layer, window_size, filter_size, regularization_v):
	block_1 = Conv1D(filter_size, window_size, padding='same',kernel_initializer=inizialiser,
					 kernel_regularizer=l2(regularization_v), activation='relu')(input_layer)
	block_1 = BatchNormalization()(block_1)
	block_1 = Dropout(0.1)(block_1)

	block_1 = Conv1D(filter_size, window_size, padding='same',kernel_initializer=inizialiser,
					 kernel_regularizer=l2(regularization_v), activation='relu')(block_1)
	block_1 = BatchNormalization()(block_1)
	block_1 = Dropout(0.1)(block_1)

	connection = add([block_1, input_layer])
	return connection


def build_model(regularization_v, inner_filters):
	inputs = Input(shape=(None, 4))

	layer_i = conv1d_module(inputs, 25, 64, regularization_v)
	layer_sm = shortcut_module(layer_i, 20, inner_filters // 2, regularization_v)

	layer1_nsm = no_shortcut_module(layer_sm, 15, inner_filters // 2, regularization_v)
	layer2_nsm = no_shortcut_module(layer1_nsm, 15, inner_filters // 2, regularization_v)
	layer3_nsm = no_shortcut_module(layer2_nsm, 15, inner_filters // 2, regularization_v)

	merged_layers = concatenate([layer1_nsm, layer2_nsm,layer3_nsm], axis=-1)
	merged_layers_o = shortcut_module(merged_layers, 10, inner_filters, regularization_v)
	merged_layers_o1 = no_shortcut_module(merged_layers_o, 5, inner_filters, regularization_v)
	merged_layers_o2 = no_shortcut_module(merged_layers_o1, 5, inner_filters, regularization_v)
	merged_layers_o3 = no_shortcut_module(merged_layers_o2, 5, inner_filters, regularization_v)

	last_1d = conv1d_module(merged_layers_o3, 10, 64, regularization_v)
    
	last_cnn = Conv1D(64, 5, padding='same', activation='relu',
	                  kernel_initializer=inizialiser, name='last_CNN_layer')(last_1d)

	layer = GlobalAveragePooling1D()(last_cnn)
    
	layer_d1 = Dense(512, activation='relu')(layer)
	layer_d1 = BatchNormalization()(layer_d1)
	layer_d1 = Dropout(0.1)(layer_d1)

	layer_d2 = Dense(128, activation='relu')(layer_d1)
	layer_d2 = BatchNormalization()(layer_d2)
	layer_d2 = Dropout(0.1)(layer_d2)
  
	layer_o = Dense(1, activation='sigmoid')(layer_d2)

	final_model = Model(inputs=inputs,outputs=layer_o)
	final_model.summary()
	return final_model


def transfer_learning_model():
    input_layer = Input(shape=(None, 4))
    model_file = 'models/tensorflow1/tl_model_tf1/model_compiled.h5'
    base_model = load_model(model_file, compile=False)
    base_model.trainable = True # for inference mode
    #base_model.trainable = False # completelly frozen
    deepmodel = Model(inputs=base_model.inputs, outputs=base_model._layers_by_depth[7][0].output)
    #for l in deepmodel.layers: # completelly frozen
    #    l.trainable = False # completelly frozen
    #deepmodel.trainable = False # completelly frozen

    output_ = deepmodel(input_layer)#, training=False)

    layer = Dense(1024, activation='relu')(output_)
    layer = BatchNormalization()(layer)
    layer = Dropout(0.1)(layer)

    layer = Dense(32, activation='relu')(layer)
    layer = BatchNormalization()(layer)
    layer = Dropout(0.1)(layer)
    
    
    layer = Dense(1, activation='sigmoid')(layer)
    tl_model = Model(inputs=visible, outputs=layer)
    tl_model.summary()
    return tl_model
