#!/usr/bin/env python

import os
import glob
import lycon
import cv2
import pandas as pd
import tarfile
import shutil
import datetime
import numpy as np
import math
import pickle
import yaml

import lib.read_settings as read_settings
import lib.datasets as datasets
import lib.models as models
import matplotlib.pyplot as plt
import tensorflow as tf


gpus = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_memory_growth(gpus[0], True)
tf.config.experimental.set_virtual_device_configuration(gpus[0],
    [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=30480)])

#################################### Settings ####################################
## Read settings
global_settings = read_settings.check_global('config.yaml')
cnn_settings = read_settings.check_cnn('config.yaml')
with open(r'config.yaml') as config_file:
    cfg = yaml.safe_load(config_file)

# Input data
dataset = cfg['dataset']
data_dir = cfg['base_dir']
frac = global_settings['input_data']['frac']

# Random state
random_state = global_settings['random_state']

# CNN settings
batch_size    = cnn_settings['data']['batch_size']
px_del        = cnn_settings['data']['px_del']
preserve_size = cnn_settings['data']['preserve_size']
augment       = cnn_settings['data']['augment']
use_weights   = cnn_settings['data']['use_weights']
weights       = cnn_settings['data']['weights']

fc_layers_nb          = cnn_settings['architecture']['fc_layers_nb']
fc_layers_size        = cnn_settings['architecture']['fc_layers_size']
fc_layers_dropout     = cnn_settings['architecture']['fc_layers_dropout']
classif_layer_dropout = cnn_settings['architecture']['classif_layer_dropout']
train_fe              = cnn_settings['architecture']['train_fe']

lr_method  = cnn_settings['compilation']['lr_method']
initial_lr = cnn_settings['compilation']['initial_lr']
decay_rate = cnn_settings['compilation']['decay_rate']
loss       = cnn_settings['compilation']['loss']

resume  = cnn_settings['training']['resume']
epochs  = cnn_settings['training']['epochs']
workers = cnn_settings['training']['workers']


## Output
os.makedirs('model_weights', exist_ok=True)
# Generate output directory pattern
if use_weights: # if using weigths 
    output_dir_patt = os.path.join('model_weights', '_'.join(['cnn', 'w', dataset]))
else: # if not using weigths 
    output_dir_patt = os.path.join('model_weights', '_'.join(['cnn', 'nw', dataset]))

# Look for previous outputs with same pattern
prev_output = glob.glob(output_dir_patt + '*')
prev_output.sort()

# Case of not resuming from previous CNN training
if not resume:        
    # Create a new output directory
    output_dir = '_'.join([output_dir_patt, datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')])
    print(f'Creating new output directory {output_dir}')
    os.mkdir(output_dir)

# Case of resuming from previous CNN training
else:
    # Check that settings are similar as previous ones
    read_settings.check_previous_cnn_settings(global_settings, cnn_settings, prev_output)
    
    # Look for saved model in most recent output
    saved_models = glob.glob(os.path.join(prev_output[-1], 'model.last.epoch.*.hdf5'))
    saved_models.sort()
    
    # Case of an existing previous output to resume from
    if len(saved_models) > 0:
        # Set output dir to the most recent previous output
        output_dir = prev_output[-1]
        # Choose most recent model
        saved_model = saved_models[-1]
        print(f'Resuming training from {output_dir}, found saved model {os.path.basename(saved_model)}')
        
    # Case of no existing previous output to resume from
    else:
        # Create a new output directory
        output_dir = '_'.join([output_dir_patt, datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')])
        print(f'No previous training to resume from, creating new output directory {output_dir}')
        os.mkdir(output_dir)
        
# Write settings to output directory
read_settings.write_cnn_settings(global_settings, cnn_settings, output_dir)


##################################################################################

## Read data for CNN
df_train, df_valid, _, df_classes, df_comp = datasets.read_data_cnn(
    path=os.path.join(data_dir, 'rhizaria_learning_set.csv'),
    frac=frac,
    random_state=random_state
)

# Write dataset composition to output_dir
df_comp.to_csv(os.path.join(output_dir, 'df_comp.csv'), index=True)

# Number of plankton classes to predict
nb_classes = len(df_classes)


# Generate class weights
class_weights = None
if use_weights:
    class_counts = df_train.groupby('classif_id').size()
    count_max = 0
    class_weights = {}
    for idx in class_counts.items():
        count_max = (idx[1], count_max) [idx[1] < count_max]
    for i,idx in enumerate(class_counts.items()):
        if weights == 'i_f': # Weights computed with inverse frequency
            class_weights.update({i : count_max / idx[1]})
        elif weights == 'sqrt_i_f': # # Weights computed with square root of inverse frequency
            class_weights.update({i : math.sqrt(count_max / idx[1])})

train_batches = datasets.EcoTaxaGenerator(
    images_paths = df_train['path_to_img'].values,
    input_shape = (224, 224, 3),
    labels = df_train['classif_id'].values,
    classes = df_classes.classif_id.tolist(),
    batch_size = batch_size,
    augment = True,
    shuffle = True,
    crop = [0,0,31,0]
    )

valid_batches = datasets.EcoTaxaGenerator(
    images_paths = df_valid['path_to_img'].values,
    input_shape = (224, 224, 3),
    labels = df_valid['classif_id'].values,
    classes = df_classes.classif_id.tolist(),
    batch_size = batch_size,
    augment = False,
    shuffle = False,
    crop = [0,0,31,0]
    )

for image_batch, label_batch in train_batches:
    print('Image batch shape: ', image_batch.shape)
    print('Label batch shape: ', label_batch.shape)
    break

    
## Case of not resuming from previous training
if not resume:
    ## Generate CNN
    my_cnn = models.create_cnn(
        fc_layers_nb,
        fc_layers_dropout, 
        fc_layers_size, 
        classif_layer_dropout, 
        classif_layer_size=nb_classes, 
        train_fe=train_fe, 
        glimpse=True
    )
    
    ## Compile CNN
    my_cnn = models.compile_cnn(
        my_cnn, 
        lr_method=lr_method, 
        initial_lr=initial_lr, 
        steps_per_epoch=len(train_batches), 
        decay_rate=decay_rate, 
        loss=loss
    )
    
    # Set first epoch to 0
    initial_epoch = 0
    # Declare no previous training history
    prev_history = None


## Case of resuming from previous training
else:
    my_cnn, initial_epoch, prev_history = models.load_cnn(saved_model, glimpse = True)


## Train CNN
history, best_epoch = models.train_cnn(
    model=my_cnn, 
    prev_history=prev_history,
    train_batches=train_batches, 
    valid_batches=valid_batches, 
    batch_size=batch_size, 
    initial_epoch=initial_epoch,
    epochs=epochs+initial_epoch, 
    class_weights=class_weights, 
    output_dir=output_dir,
    workers=workers,
)

    
## Predict test batches and evaluate CNN
models.predict_evaluate_cnn(
    model=my_cnn, 
    best_epoch=best_epoch,
    batches=valid_batches, 
    true_classes = np.array(df_valid.classif_id.tolist()),
    df_classes=df_classes, 
    output_dir=output_dir,
    workers=workers,
)
