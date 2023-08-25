#!/usr/bin/env python
# coding: utf-8


""" 
This is a script to predict a classification for vignettes exported from ecotaxa.
"""


import os
import glob
import cv2
import pandas as pd
import tarfile
import shutil
import datetime
import numpy as np
import math
import pickle
import time
import sys
from pyarrow.parquet import read_table
import yaml
from collections import Counter

import matplotlib.pyplot as plt
import tensorflow as tf
import tensorflow_hub as hub

import lib.read_settings as read_settings
import lib.datasets as datasets
import lib.models as models

# set a memory limit on tensorflow, to allow others to use the GPU too
# it may be slightly less efficient memory-wise
gpus = tf.config.experimental.list_physical_devices('GPU')
tf.config.experimental.set_virtual_device_configuration(gpus[0],
    [tf.config.experimental.VirtualDeviceConfiguration(memory_limit=20480)])

###################################### Settings ######################################
# read config
with open(r'config.yaml') as config_file:
    cfg = yaml.safe_load(config_file)

batch_size = 32  
workers = cfg['cnn']['training']['workers']
data_dir = 'data'
######################################################################################



# Read objects extracted from EcoTaxa to be predicted
df = read_table(os.path.join(data_dir, 'all_rhizaria.parquet')).to_pandas()
# Ignore validated and dubious objects
df = df[df['classif_qual'] == 'P'].reset_index(drop = True)
# Fix path to image
df['path_to_img'] = [os.path.join('data/images', p) for p in df['path_to_img']]


## Choose training output
# Find last training output
weights_dirs = glob.glob('model_weights/*')
weights_dirs.sort()
weights_dir = weights_dirs[-1]
weights_file = glob.glob(os.path.join(weights_dir, 'model.last.epoch*'))[0]

# List of classes
classes = pd.read_csv(os.path.join(weights_dir, 'df_comp.csv')).classif_id.tolist()

## Load CNN model
my_cnn = tf.keras.models.load_model(weights_file, custom_objects={'KerasLayer':hub.KerasLayer}, compile = False)
# get model input shape
input_shape = my_cnn.layers[0].input_shape
# remove the None element at the start (which is where the batch size goes)
input_shape = tuple(x for x in input_shape if x is not None)
my_cnn.summary()

# Prepare batches
batches = datasets.EcoTaxaGenerator(
    images_paths = df['path_to_img'].values,
    input_shape = input_shape,
    labels = None,
    classes = None,
    batch_size = batch_size,
    augment = False,
    shuffle = False,
    crop = [0,0,31,0]
    )

for image_batch, label_batch in batches:
    print('Image batch shape: ', image_batch.shape)
    #print('Label batch shape: ', label_batch.shape)
    break


## Run prediction
# Predict test batches and convert logits to plankton classes and scores
logits = my_cnn.predict(batches, max_queue_size=max(10, workers*2))
predicted_classes = np.array(classes)[np.argmax(logits, axis=1)]
scores = (np.max(tf.nn.softmax(logits), axis=1))

# Count predictions in each classes
Counter(predicted_classes)

## Inspect predictions in a given class
#for i in np.where(predicted_classes == 'juvenile<Euchaeta')[0]:
#    batch_nb = i // batch_size
#    plt.imshow(batches[batch_nb][0][i - batch_nb*batch_size][:,:,0], cmap = 'gray')
#    plt.title(predicted_classes[i])
#    plt.show()


# Store this in a dataframe
df_pred = pd.DataFrame()
df_pred['object_id'] = df['orig_id']
df_pred['object_annotation_category'] = predicted_classes
df_pred['object_annotation_status'] = 'predicted'
df_pred['object_cnn_score'] = scores
df_pred['object_cnn_precition'] = 'True'

# first_row as Dataframe row with appropriate headers
first_row = ['[t]', '[t]', '[t]', '[f]', '[t]']
first_row = pd.DataFrame(first_row).T
first_row.columns = df_pred.columns

# concat first_row and dataframe
df_pred = pd.concat([first_row, df_pred], ignore_index = True)

# save table
df_pred.to_csv('data/ecotaxa_rhizaria_predictions.tsv', index=False, sep='\t', header=True)

