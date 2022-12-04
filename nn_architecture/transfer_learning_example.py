### CODE CREATED BY MARCO SALVATORE ###

### THE CODE IS INTENDED TO GIVE AN EXAMPLE ON HOW TO USE THE PRE-TRAINED MODEL AND FINE-TUNE WITH NEW DATA
### THE CODE NEED SEQUENCE DATA IN ONE-HOT ENCODING FORMAT

### STEPS TO FOLLOW:

  ### 1- LOAD PRE-TRAINED MODEL
  ### 2- DECIDE AT WHICH LEVEL YOU WANT TO FINE-TUNED THE ORIGINAL ARCHITECTURE (SEE COMMENTS IN transfer_learning_model()
  ### 3- LOAD THE FILE WITH SEQUENCES AND LABELS OF INTEREST 
  ### 4- RUN TRAINING/VALIDATION/TEST OR K-FOLD CROSS VALIDATION 

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
from sklearn.model_selection import train_test_split

inizialiser = "he_normal"
regularization_v = 1e-4

# Function to load the pre-trained model and add/change layers. 

def transfer_learning_model():
    input_layer = Input(shape=(None, 4))
    model_file = 'models/tensorflow1/tl_model_tf1/model_compiled.h5'
    base_model = load_model(model_file, compile=False)
    
    # FOR ADDITIONAL REFERENCE FOLLOW https://keras.io/guides/transfer_learning/
    # First freeze the model completely and do a first Fine-tuning of the model. Check stability and improvement of the model.
    #base_model.trainable = False # completelly frozen
    
    # Then run in "inference mode" ---> this is used to do a round of fine-tuning of the entire model and to try to improve performance even more
    base_model.trainable = True     
    
    # DECIDE AT WHICH LAYER OF THE NETWORK THE MODEL NEEDS TO BE FINE-TUNED (THE FINE-TUNING CAN BE DONE IN DIFFERENT WAYS)
    # FOR ADDITIONAL REFERENCE FOLLOW https://keras.io/guides/transfer_learning/
    # THIS IS THE MOST IMPORTANT LINE TO TWICK AND/OR TAKE INTO CONSIDERATION
    
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
    
    # IN THE ORIGINAL PUBLICATION WE FINE-TUNED FOR BINARY (CELL TYPE SPECIFIC CLASSIFICATION) 
    # YOU CAN CHANGE THIS LINE OF CODE IF YOU WISH TO DO MULTI-CLASS BY CHANGING NUMBER OF OUTPUTS NODES
    
    layer = Dense(1, activation='sigmoid')(layer)
    tl_model = Model(inputs=visible, outputs=layer)
    tl_model.summary()
    return tl_model
  
  
# Change file npz with your preferred files. The input files should one-hot encoded sequences.
# You can decide how to do train/val/test splits or CV (Note in the original publication, we have done 3-CV and test on Chr2/Chr3 in all the comparisons)

data = np.load('file.npz', allow_pickle=True, encoding='bytes')

# this is just an example: 

Y = data['label_'] # LABEL FOR THE SEQUENCES
X = data['sequence_'] # ONE-HOT ENCODING SEQUENCES 

X_train, X_val, y_train, y_val = train_test_split(X, Y, test_size=0.33, random_state=42)

learning_rate = 0.000005 

class_weights = class_weight.compute_class_weight('balanced',
                                                  np.unique(y_train),
                                                  y_train)

# Call the model for fine-tuning the pre-trained model 
model = deepmodel()

opt = Adam(lr=learning_rate)
model.save(specify_model_file)

# IN THE ORIGINAL PUBLICATION WE FINE-TUNED FOR BINARY (CELL TYPE SPECIFIC CLASSIFICATION) 
# YOU CAN CHANGE THIS LINE OF CODE IF YOU WISH TO DO MULTI-CLASS BY ADDING CATEGORICAL_CROSSENTROPY
# YOU CAN CHANGE THIS LINE OF CODE IF YOU WISH TO DO REGRESSION 

model.compile(optimizer=opt, loss='binary_crossentropy', metrics=METRICS)


# Tensorflow/Keras in built functions to monitor training of the model.

checkpointer = ModelCheckpoint(filepath=model_file3,
                               verbose=0,
                               save_best_only=True)

earlyStopping = keras.callbacks.EarlyStopping(monitor='val_loss',
                                              patience=10, verbose=1,
                                              mode='auto')
lr_reducer = keras.callbacks.ReduceLROnPlateau(monitor='val_loss',
                                               factor=0.1, patience=5,
                                               verbose=0, mode='auto',
                                               min_delta=0.0001,
                                               cooldown=0, min_lr=0)

# Save log file with info
logger_csv = keras.callbacks.CSVLogger(logger_file, separator=",", append=False)

# Use tensorboard to monitor the model training
tensorboard = TensorBoard(log_dir=log_file,
                          histogram_freq=0,
                          write_graph=True,
                          write_images=False)

callbacks = [checkpointer, earlyStopping, lr_reducer, tensorboard, logger_csv]


history = model.fit(X_train, y_train, validation_data=(X_train, y_train), epochs=100, verbose=1, batch_size=128,
                    class_weight=class_weights, 
                    callbacks=callbacks).history
