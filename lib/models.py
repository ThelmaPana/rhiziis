import pandas as pd
import numpy as np
import itertools as it
import os
import pickle
import glob

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, recall_score, log_loss

import tensorflow as tf
import tensorflow_hub as hub
os.environ['TFHUB_CACHE_DIR'] = '.tf_models'
import tensorflow_addons as tfa
from tensorflow.keras import layers, optimizers, losses, callbacks 
import tensorflow_addons as tfa


def gridsearch_rf(df_train, df_valid, classes, eval_metric, max_features_try, min_samples_leaf_try, n_estimators_try, output_dir, n_jobs, class_weights=None, random_state=None):
    """
    Perform a grid search to find best hyperparameters for random forest model.
    
    Args:
        df_train (DataFrame): training data to use to fit grid search
        df_valid (DataFrame): validation data to use to evaluate grid search
        classes (list, array): name of classes
        eval_metric (str): metric to use for hyperparameters selection ('accuracy', 'balanced_accuracy' or 'log_loss')
        max_features_try (list): tries for number of variables per node; default sqrt(nb of vars)
        min_samples_leaf_try (list): tries for min number of objects in leaf; default for classif = 5
        n_estimators_try (list): tries for number of estimators (usually between 100 and 500)
        output_dir (str): directory where to save gridsearch results
        n_jobs (int): number of cores to use 
        class_weights (dict): weights for classes
        random_state (int or RandomState): controls both the randomness of the bootstrapping and features sampling; default=None

    
    Returns:
        results (DataFrame): results of grid search
        best_params (dict): best parameters based on evaluation metric
    """
    
    # Shuffle data
    df_train = df_train.sample(frac=1, random_state=random_state).reset_index(drop=True)
    df_valid = df_valid.sample(frac=1, random_state=random_state).reset_index(drop=True)
    
    # Split data and labels
    y_train = df_train['classif_id']
    X_train = df_train.drop('classif_id', axis=1)
    
    y_valid = df_valid['classif_id']
    X_valid = df_valid.drop('classif_id', axis=1)
    
    # Prepare one-hot encoding to compute log-loss for validation data
    mlb = MultiLabelBinarizer(classes=classes)
    y_true = mlb.fit_transform([[l] for l in y_valid])
    
    # Build grid of hyperparameters to explore
    grid = {
        'max_features': max_features_try, 
        'min_samples_leaf': min_samples_leaf_try,
    }
    # Make a list of all parameters combinations
    keys = list(grid.keys())
    grid_list = list(it.product(*(grid[key] for key in keys)))

    # Initiate empty dict for results
    results = {
        'n_estimators': [],
        'max_features': [],
        'min_samples_leaf': [],
        'valid_accuracy': [],
        'valid_balanced_accuracy':[],
        'valid_log_loss':[]
    }

    # First loop on parameters other than n_estimators
    for max_features, min_samples_leaf in grid_list:
        print(f"Trying parameters max_features = {max_features} and min_samples_leaf = {min_samples_leaf}.")
    
        # Initiate a RF model with warm start
        rf = RandomForestClassifier(
            criterion='entropy', 
            min_samples_split=2, 
            max_features=max_features, 
            min_samples_leaf=min_samples_leaf,
            warm_start=True,
            n_jobs=n_jobs,
            class_weight=class_weights,
            random_state=random_state
        )
        
        # Second loop on n_estimators
        for n_estimators in n_estimators_try:
            print(f"Number of estimators = {n_estimators}")
            
            # Set number of estimators in RF model
            rf.n_estimators = n_estimators
            
            # Fit on training data
            rf.fit(X=X_train, y=y_train)
            
            # Compute accuracy on validation data
            valid_accuracy = accuracy_score(y_valid, rf.predict(X_valid))
            valid_balanced_accuracy = balanced_accuracy_score(y_valid, rf.predict(X_valid))
            
            # Compute log loss on validation data  
            # log_loss only accepts weights as sample_weights and not as class_weights, compute sample_weights
            if class_weights is not None:
                sample_weight = [class_weights[c] for c in y_valid]
            else:
                sample_weight = None
            valid_log_loss = log_loss(y_true, rf.predict_proba(X_valid), sample_weight=sample_weight)

            # Store results in dict
            results['n_estimators'].append(n_estimators)
            results['max_features'].append(max_features)
            results['min_samples_leaf'].append(min_samples_leaf)
            results['valid_accuracy'].append(valid_accuracy) 
            results['valid_balanced_accuracy'].append(valid_balanced_accuracy) 
            results['valid_log_loss'].append(valid_log_loss) 

    # Write gridsearch results
    with open(os.path.join(output_dir, 'train_results.pickle'),'wb') as results_file:
        pickle.dump(results, results_file)
        
    # Convert to datfarame
    results = pd.DataFrame(results)
    
    # Extract best parameters based on evaluation metric value on validation data
    if eval_metric == 'log_loss':
        # if evaluation metric is log loss, look for the smallest value
        best_params = results.nsmallest(1, 'valid_'+ eval_metric).reset_index(drop=True).drop(['valid_accuracy', 'valid_balanced_accuracy', 'valid_log_loss'], axis=1)
    else:
        # in other cases, look for the largest value
        best_params = results.nlargest(1, 'valid_'+ eval_metric).reset_index(drop=True).drop(['valid_accuracy', 'valid_balanced_accuracy', 'valid_log_loss'], axis=1)
    best_params = best_params.iloc[0].to_dict()
        
    # Print GS results
    print(f"Hyperparameters selection using {eval_metric} value.")
    print(f"Selected hyperparameters are: n_estimators = {best_params['n_estimators']}, max_features = {best_params['max_features']}, min_samples_leaf = {best_params['min_samples_leaf']}")

    return results, best_params


def train_rf(df, n_estimators, max_features, min_samples_leaf, n_jobs, class_weights, random_state=None):
    """
    Fit a random forest model on data.
    
    Args:
        df (DataFrame): data to use for training
        tree_nb (int): number of trees for the RF model
        max_features (int): number of variables per node
        min_samples_leaf (int): min number of objects in leaf
        n_jobs (int): number of cores to use 
        class_weights(dict): weights for classes
        random_state (int or RandomState): controls both the randomness of the bootstrapping and features sampling; default=None
    
    Returns:
        rf (RandomForestClassifier): fitted random forest model
    """
    
    # Shuffle data
    df = df.sample(frac=1).reset_index(drop=True)
    
    # Split data and labels
    y_train = df['classif_id']
    X_train = df.drop('classif_id', axis=1)
    
    # Initiate RF model
    rf = RandomForestClassifier(
        n_estimators=n_estimators, 
        criterion='entropy', 
        min_samples_split=2, 
        min_samples_leaf=min_samples_leaf, 
        max_features=max_features,
        n_jobs=n_jobs,
        class_weight=class_weights,
        random_state=random_state
    )
    
    # Fit the RF model
    rf = rf.fit(X=X_train, y=y_train)
    
    return rf


def predict_evaluate_rf(rf_model, df, df_classes, output_dir):
    """
    Evaluate a random forest model.
    
    Args:
        rf_model (RandomForestClassifier): random forest model to evaluate
        df (DataFrame): data to use for model evaluation
        df_classes (DataFrame): dataframe of classes with living attribute
        output_dir (str): directory where to save prediction results
    
    Returns:
        nothing
    """
    
    # Split data and labels
    y = df['classif_id'].tolist()
    y = np.array(y)
    X = df.drop('classif_id', axis=1)

    # Make a list of classes
    classes = df_classes['classif_id'].tolist()
    classes.sort()
    classes = np.array(classes)
    
    # and of regrouped classes
    classes_g = df_classes['classif_id_2'].tolist()
    classes_g = list(set(classes_g))
    classes_g.sort()
    classes_g = np.array(classes_g)
    
    # Make a list of plankton classes
    plankton_classes = df_classes[df_classes['plankton']]['classif_id'].tolist()
    plankton_classes = np.array(plankton_classes)
    
    # Make a list of plankton classes for grouped classes
    plankton_classes_g = df_classes[df_classes['plankton']]['classif_id_2'].tolist()
    plankton_classes_g = np.array(plankton_classes_g)
    
    # Predict test data
    y_pred = rf_model.predict(X)
    
    # Compute accuracy between true labels and predicted labels
    accuracy = accuracy_score(y, y_pred)
    balanced_accuracy = balanced_accuracy_score(y, y_pred)
    plankton_precision = precision_score(y, y_pred, labels=plankton_classes, average='weighted', zero_division=0)
    plankton_recall = recall_score(y, y_pred, labels=plankton_classes, average='weighted', zero_division=0)
    
    # Display results
    print(f'Test accuracy = {accuracy}')
    print(f'Balanced test accuracy = {balanced_accuracy}')
    print(f'Weighted plankton precision = {plankton_precision}')
    print(f'Weighted plankton recall = {plankton_recall}')
     
    ## Now do the same after regrouping objects to larger classes
    # Generate taxonomy match between taxo used for classif and larger ecological classes 
    taxo_match = df_classes.set_index('classif_id').to_dict('index')
    
    # Convert true classes to larger ecological classes
    y_g = np.array([taxo_match[t]['classif_id_2'] for t in y])
    
    # Convert predicted classes to larger ecological classes
    y_pred_g = np.array([taxo_match[p]['classif_id_2'] for p in y_pred])
    
    # Compute accuracy, precision and recall for living classes and loss from true labels and predicted labels
    accuracy_g = accuracy_score(y_g, y_pred_g)
    balanced_accuracy_g = balanced_accuracy_score(y_g, y_pred_g)
    plankton_precision_g = precision_score(y_g, y_pred_g, labels=plankton_classes_g, average='weighted', zero_division=0)
    plankton_recall_g = recall_score(y_g, y_pred_g, labels=plankton_classes_g, average='weighted', zero_division=0)
    
    # Display results
    print(f'Grouped test accuracy = {accuracy_g}')
    print(f'Grouped balanced test accuracy = {balanced_accuracy_g}')
    print(f'Grouped weighted plankton precision = {plankton_precision_g}')
    print(f'Grouped weighted plankton recall = {plankton_recall_g}')

    # Write classes and test metrics into a test file
    with open(os.path.join(output_dir, 'test_results.pickle'),'wb') as test_file:
        pickle.dump({
            'classes': classes,
            'classes_g': classes_g,
            'plankton_classes': plankton_classes,
            'plankton_classes_g': plankton_classes_g,
            'true_classes': y,
            'predicted_classes': y_pred,
            'true_classes_g': y_g,
            'predicted_classes_g': y_pred_g,
            'accuracy': accuracy,
            'balanced_accuracy': balanced_accuracy,
            'plankton_precision': plankton_precision,
            'plankton_recall': plankton_recall,
            'accuracy_g': accuracy_g,
            'balanced_accuracy_g': balanced_accuracy_g,
            'plankton_precision_g': plankton_precision_g,
            'plankton_recall_g': plankton_recall_g,
        },
        test_file)

    pass
    

def create_cnn(fc_layers_nb, fc_layers_dropout, fc_layers_size, classif_layer_dropout, classif_layer_size,  train_fe = False, glimpse = True):

    """
    Generates a CNN model. 
    
    Args:
        fc_layers_nb (int): number of fully connected layers 
        fc_layers_dropout (float): dropout of fully connected layers
        fc_layers_size (int): size of fully connected layers 
        classif_layer_dropout (float): dropout of classification layer
        classif_layer_size (int): size of classification layer (i.e. number of classes to predict)
        train_fe (bool): whether to train the feature extractor (True) or only classification head (False)
        glimpse(bool): whether to show a model summary
    
    Returns:
        model (tensorflow.python.keras.engine.sequential.Sequential): CNN model
    """
    
    ## Initiate empty model
    model = tf.keras.Sequential()
    
    ## MobileNet V2 feature extractor
    fe_url = "https://tfhub.dev/google/imagenet/mobilenet_v2_140_224/feature_vector/4"
    fe_layer = hub.KerasLayer(fe_url, input_shape=(224, 224, 3))
    # Set feature extractor trainability
    fe_layer.trainable=train_fe
    model.add(fe_layer)
    
    ## Fully connected layers
    if fc_layers_nb:
        for i in range(fc_layers_nb):
            if fc_layers_dropout:
                model.add(layers.Dropout(fc_layers_dropout))
            model.add(layers.Dense(fc_layers_size, activation='relu'))
    
    ## Classification layer
    if classif_layer_dropout:
        model.add(layers.Dropout(classif_layer_dropout))
    model.add(layers.Dense(classif_layer_size))
    
    if glimpse:
        model.summary()

    return model


def compile_cnn(model, lr_method, initial_lr, steps_per_epoch, decay_rate=None, loss='cce'):
    """
    Compiles a CNN model. 
    
    Args:
        model (tensorflow.python.keras.engine.sequential.Sequential): CNN model to compile
        lr_method (str): method for learning rate. 'constant' for a constant learning rate, 'decay' for a decay
        initial_lr (float): initial learning rate. If lr_method is 'constant', set learning rate to this value
        steps_per_epochs (int): number of training steps at each epoch. Usually number_of_samples // batch_size or len(train_batches)
        decay_rate (float): rate for learning rate decay
        loss (str): method to compute loss. 'cce' for CategoricalCrossentropy (see https://www.tensorflow.org/api_docs/python/tf/keras/losses/CategoricalCrossentropy), 'sfce' for SigmoidFocalCrossEntropy (see https://www.tensorflow.org/addons/api_docs/python/tfa/losses/SigmoidFocalCrossEntropy), usefull for unbalanced classes
    
    Returns:
        model (tensorflow.python.keras.engine.sequential.Sequential): compiled CNN model
    """
    # TODO if lr_method='decay', decay_rate in mandatory

    ## Learning rate
    if lr_method == 'decay':
        lr_schedule = optimizers.schedules.InverseTimeDecay(
            initial_learning_rate=initial_lr, 
            decay_steps=steps_per_epoch, 
            decay_rate=decay_rate, 
            staircase=False
        )
    else: # Keep constant learning rate
        lr_schedule = initial_lr

    ## Optimizer: use Adam
    optimizer = optimizers.Adam(learning_rate=lr_schedule)
    
    ## Loss
    if loss == 'cce':
        loss = losses.CategoricalCrossentropy(from_logits=True,reduction=losses.Reduction.SUM_OVER_BATCH_SIZE)
    elif loss == 'sfce':
        loss = tfa.losses.SigmoidFocalCrossEntropy(from_logits=True,reduction=losses.Reduction.SUM_OVER_BATCH_SIZE)
    
    
    model.compile(
      optimizer=optimizer,
      loss=loss,
      metrics=['accuracy']
    )
    
    return model


def load_cnn(saved_model, glimpse=True):
    """
    Load a CNN model from a previous model, return index of last epoch and training history
    
    Args:
        saved_model (str): path to model to load
        glimpse(bool): whether to show a model summary
    
    Returns:
        model (tensorflow.python.keras.engine.sequential.Sequential): CNN model
        initial_epoch (int): last training epoch
        history (dict): training history
    """
    
    # Extract number of last training epoch
    initial_epoch = int(saved_model.split('.')[-2])
    
    # Recreate the exact same model, including its weights and the optimizer
    model = tf.keras.models.load_model(saved_model, custom_objects={'KerasLayer':hub.KerasLayer})
    
    if glimpse:
        model.summary()
        
    # Load previous history
    with open(os.path.join(os.path.dirname(saved_model), 'train_results.pickle'),'rb') as file:
        history = pickle.load(file)
    
    return model, initial_epoch, history


def train_cnn(model, prev_history, train_batches, valid_batches, batch_size, initial_epoch, epochs, class_weights, output_dir, workers):
    """
    Trains a CNN model. 
    
    Args:
        model (tensorflow.python.keras.engine.sequential.Sequential): CNN model to train
        prev_history (dict or None): previous training history
        train_batches: batches of training data
        valid_batches: batches of validation data
        batch_size (int): size of batches
        initial_epoch (int): epoch to start from
        epochs (int): number of epochs to train for
        class_weight(dict): weights for classes
        output_dir (str): directory where to save model weights
        workers (int): number of parallel threads for data generators
    
    Returns:
        history (tensorflow.python.keras.callbacks.History): training history
        best_epoch (int): index of best epoch
    """
    ## Set callbacks
    # Checkpoint callback
    cp_callback = callbacks.ModelCheckpoint(
        filepath=os.path.join(output_dir, 'weights.{epoch:02d}.hdf5'),
        monitor='val_loss',
        save_best_only=False,
        mode='min',
        save_weights_only=False,
        save_freq='epoch',
        verbose=1)
    
    # Learning rate callback
    lr_history = [] # initiate empty list to store lr values
    class LearningRateLoggingCallback(callbacks.Callback):
        # Get learning rate value at the beginning of each epoch
        def on_epoch_begin(self, epoch, logs=lr_history):
            optimizer = self.model.optimizer # get optimizer
            lr_value = optimizer.lr(optimizer.iterations).numpy() # get lr value
            lr_history.append(lr_value) # append lr value to list of lr values
            print(f'LR is {lr_value}') # print lr value
    
    # Fit the model
    history = model.fit(
        train_batches, 
        epochs=epochs,
        initial_epoch=initial_epoch,
        validation_data=valid_batches,
        callbacks=[cp_callback, LearningRateLoggingCallback()],
        class_weight=class_weights,
        max_queue_size=max(10, workers*2),
        workers=workers
    )
    training_history = history.history
    
    # Add learning rate logging to training history
    training_history['lr'] = lr_history
    
    # Look for previous history to merge with
    if prev_history:
        for k in prev_history.keys():
            prev_history[k].extend(training_history[k])
        training_history = prev_history
    
    # Write training history 
    with open(os.path.join(output_dir, "train_results.pickle"),"wb") as results_file:
        pickle.dump(training_history, results_file)
    
    # Save the last epoch to resume training if needed
    model.save(os.path.join(output_dir, 'model.last.epoch.' + str(epochs).zfill(3) + '.hdf5'))
    
    # Compute index of best epoch based on validation accuracy
    # NB: add 1 because python indexing starts from 0 and epochs start from 1
    best_epoch = np.argmin(training_history['val_loss']) + 1
    
    return history, best_epoch


def predict_evaluate_cnn(model, best_epoch, batches, true_classes, df_classes, output_dir, workers):
    """
    Predict batches and evaluate a CNN model.
    Predict images from test set and compute accuracy, balanced_accuracy, precision and recall on relevant classes.
    Regroup objects into larger ecological classes and re-evaluate model.
    
    Args:
        model (tensorflow.python.keras.engine.sequential.Sequential): CNN model to evaluate
        best_epoch(int): index of the best epoch, used to test the model
        batches (datasets.DataGenerator): batches of test data to predict
        true_classes (list): true classes of test images
        df_classes (DataFrame): dataframe of classes with living attribute
        output_dir (str): directory where to save prediction results
        workers(int): number of parallel threads for data generators

    Returns:
        Nothing
    """

    # Make a list of classes
    classes = df_classes['classif_id'].tolist()
    classes.sort()
    classes = np.array(classes)
    
    # Make a list of plankton classes
    plankton_classes = df_classes[df_classes['plankton']]['classif_id'].tolist()
    plankton_classes = np.array(plankton_classes)
    
    # Load best epoch weights to CNN model
    model.load_weights(os.path.join(output_dir, (f'weights.{best_epoch:02d}.hdf5')))
    print(f'Loading weights of epoch number {best_epoch} to evaluate CNN model.')
    
    # Predict test batches and convert predictions to plankton classes
    logits = model.predict(batches, max_queue_size=max(10, workers*2), workers=workers)
    predicted_classes = classes[np.argmax(logits, axis=1)]

    # Compute accuracy, precision and recall for living classes and loss from true labels and predicted labels
    accuracy = accuracy_score(true_classes, predicted_classes)
    balanced_accuracy = balanced_accuracy_score(true_classes, predicted_classes)
    plankton_precision = precision_score(true_classes, predicted_classes, labels=plankton_classes, average='weighted', zero_division=0)
    plankton_recall = recall_score(true_classes, predicted_classes, labels=plankton_classes, average='weighted', zero_division=0)
    
    # Display results
    print(f'Test accuracy = {accuracy}')
    print(f'Balanced test accuracy = {balanced_accuracy}')
    print(f'Weighted plankton precision = {plankton_precision}')
    print(f'Weighted plankton recall = {plankton_recall}')
        
    # Write classes and test metrics into a test file
    with open(os.path.join(output_dir, 'test_results.pickle'),'wb') as test_file:
        pickle.dump({
            'classes': classes,
            'plankton_classes': plankton_classes,
            'true_classes': true_classes,
            'predicted_classes': predicted_classes,
            'accuracy': accuracy,
            'balanced_accuracy': balanced_accuracy,
            'plankton_precision': plankton_precision,
            'plankton_recall': plankton_recall,
        },
        test_file)

    pass
