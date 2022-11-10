#keras -> stuff
from keras.layers import Dense, Dropout, Activation, Flatten
from keras.layers import Conv1D, Conv2D, MaxPooling1D, MaxPooling2D
from keras.utils import np_utils
from keras.models import load_model
from keras.callbacks import ModelCheckpoint
from keras.activations import relu
from keras.utils import plot_model
from keras.optimizers import SGD, Adam, Adamax
from keras.models import Sequential
#keras -< stuff
#keras -> scipy
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
#keras -< scipy
#keras -> others
from matplotlib import pyplot as plt
import astropy.io.fits as fits
import sys,os,time,math
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
np.random.seed(123)
#keras <- others

# import keras.callbacks as cb
# from keras.datasets import mnist
# from keras.layers.core import Activation, Dense, Dropout
# from keras.regularizers import l1, l2
# from keras.utils import np_utils


# INPUT[32x32x3] will hold the raw pixel values of the image, in this case an image 
# of width 32, height 32, and with three color channels R, G, B.
# CONV layer will compute the output of neurons that are connected to local regions 
# in the input, each computing a dot product between their weights and a small region
#  they are connected to in the input volume. This may result in volume such as [32x32x12] 
#  if we decided to use 12 filters.
# RELU layer will apply an elementwise activation function, such as the max(0, x)
# thresholding at zero. This leaves the size of the volume unchanged([32x32x12]).
# POOL layer will perform a downsampling operation along the spatial dimensions(width,
#  height), resulting in volume such as [16x16x12].
# FC(i.e. fully-connected) layer will compute the class scores, resulting in volume of 
# size[1x1x10], where each of the 10 numbers correspond to a class score, such as among 
# the 10 categories of CIFAR-10. As with ordinary Neural Networks and as the name implies,
#  each neuron in this layer will be connected to all the numbers in the previous volume.

# 1D CNN!!!! dimen = 2 -> 2D
def Gross_neuron(model_in_dim, model_out_dim, dimen=1, batch_size=128):

    if dimen == 1:
        #model_in_dim = x-axis*channels
        #model_out_dim = ¿? (rows*cols*4x4) in 2D
        activation_func1 = 'linear'#linear
        activation_func2 = 'relu'  # relu
    #    metrics = 'Accuracy'
        metrics = 'mse'
        loss_function = 'mean_squared_error'
    #    loss_function = 'mean_absolute_error'
        loss_function = 'mean_squared_logarithmic_error'
        optimizer = 'Adamax'
        dropout_rate = 0.1  # defecto 0.4 Lo pongo menor y es mejor 48%
        learning_rate = 1e-2  # 0.005
        weight_regularizer = None
        number_of_filters = batch_size#16*4  # (reduce the dimensionalyty to 36) Mayor -> mejor
        size_of_filters= 3 #(reduce the dimensionalyty to 36)
        padding = (size_of_filters - 1)/2  # same as adding padding='same'
        dense_neuron = model_in_dim#4*16  # model_out_dim
    ## Initialize model.
        model = Sequential()
        # model.add(Conv1D(batch_size, (7,), activation=activation_func)) #conv filter size is 7
        # conv filter size is 7
        # Dimensions (batch, rows, cols, channels) or (x-axis*channels)
        model.add(Conv1D(filters=number_of_filters, kernel_size=size_of_filters,
                        activation=activation_func2, padding='same', data_format='channels_last'))
        # reduces dimensionality a factor 2
        model.add(MaxPooling1D(pool_size=2, padding='same')) #4
        #model.add(Dropout(dropout_rate))
        model.add(Conv1D(filters=int(number_of_filters), kernel_size=int(size_of_filters*2),
                        activation=activation_func2, padding='same', data_format='channels_last'))  # conv filter size is 5
        # reduces dimensionality a factor 2
        model.add(MaxPooling1D(pool_size=3, padding='same')) #3
        #model.add(Dropout(dropout_rate))
        model.add(Conv1D(filters=int(number_of_filters), kernel_size=int(size_of_filters*3),
                        activation=activation_func2, padding='same', data_format='channels_last'))  # conv filter size is 5
        # reduces dimensionality a factor 2
        model.add(MaxPooling1D(pool_size=10, padding='same')) #10
        model.add(Dropout(dropout_rate))
        model.add(Flatten()) #flatten rows*channels
        model.add(Dense(dense_neuron, activation=activation_func1,
                        W_regularizer=weight_regularizer))  # inutdim es la dimension de entrada comprimida 360x4 = 
        # outdim es la dimension de matriz 4x4 en 1d = 16
        model.add(Dense(model_out_dim, activation=activation_func2))
        ## Define optimizer. we select Adamax
        ## opt = Adamax(lr=learning_rate, beta=0.9, epsilon=1e-08, decay=0.0)
        model.compile(loss=loss_function, optimizer=optimizer,metrics=[metrics], learning_rate=1e-2)
    elif dimen == 2:
        print('2d')
        activation_func1 = 'linear'  # linear
        activation_func2 = 'relu'  # relu
        metrics = 'mse'
        loss_function = 'mean_squared_logarithmic_error'
        optimizer = 'Adamax'
        dropout_rate = 0.1  # defecto 0.4 Lo pongo menor y es mejor 48%
        learning_rate = 1e-2  # 0.005
        weight_regularizer = None
        number_of_filters = batch_size
        size_of_filters = (3,3)  # (reduce the dimensionalyty to 36)
        #padding = (size_of_filters - 1)/2  # same as adding padding='same'
        dense_neuron = 16*16#model_in_dim  # 4*16  # model_out_dim
        ## Initialize model.
        model = Sequential()
        # Dimensions (batch, rows, cols, channels) 
        # this applies 16 convolution filters of size 3x3 each.
        model.add(Conv2D(filters=number_of_filters, kernel_size=size_of_filters,
                         activation=activation_func2, padding='same', data_format='channels_last'))
        # reduces dimensionality a factor 2
        model.add(MaxPooling2D(pool_size=(2,2), padding='same'))  # 4
        model.add(Dropout(dropout_rate))
        model.add(Conv2D(filters=int(number_of_filters), kernel_size=size_of_filters,
                         activation=activation_func2, padding='same', data_format='channels_last'))  # conv filter size is 5
        # reduces dimensionality a factor 2
        model.add(MaxPooling2D(pool_size=(2,2), padding='same'))  # 3
        model.add(Dropout(dropout_rate))
        model.add(Conv2D(filters=int(number_of_filters), kernel_size=size_of_filters,
                         activation=activation_func2, padding='same', data_format='channels_last'))  # conv filter size is 5
        # reduces dimensionality a factor 2
        model.add(MaxPooling2D(pool_size=(2,2), padding='same'))  # 10
        model.add(Dropout(dropout_rate))
        model.add(Flatten())  # flatten rows*channels
        model.add(Dense(dense_neuron, activation=activation_func1,
                        W_regularizer=weight_regularizer))  # inutdim es la dimension de entrada comprimida 360x4 =
        # outdim es la dimension de matriz 4x4 en 1d = 16
        model.add(Dense(model_out_dim, activation=activation_func2))
        ## Define optimizer. we select Adamax
        ## opt = Adamax(lr=learning_rate, beta=0.9, epsilon=1e-08, decay=0.0)
        model.compile(loss=loss_function, optimizer=optimizer,
                      metrics=[metrics], learning_rate=1e-2)

    else:
        pass
    return model


def TrainModel(train_in, train_out, data=None, epochs=20, batch_size=128,dimen=1):
    # #Input
    # train_in = np.zeros((n_samples, x_dimension, 4)) #OJO  x_dimension* 4 = out dimension
    # #output
    # train_out = np.zeros((n_samples, out_dimension))
    if dimen == 1:
        if train_in is None:
            print("Must provide data.")
            return
        start_time = time.time()
        model_in_dim = train_in.shape[1:]
        model_in_dim = model_in_dim[0]*model_in_dim[1]
        model_out_dim = 16
        
        model = Gross_neuron(model_in_dim, model_out_dim,
                            dimen=1, batch_size=batch_size)
        print(model_in_dim, model_out_dim)
        print('Start training.')
        history = model.fit(train_in, train_out, batch_size=batch_size,
                            epochs=epochs, verbose=1, validation_split=0.2)
        model.summary()
    elif dimen == 2:
        print('2d')
        if train_in is None:
            print("Must provide data.")
            return
        start_time = time.time()
    # Dimensions (batch, rows, cols, channels)
        model_in_dim = train_in.shape[1:]
        model_in_dim = model_in_dim[0]*model_in_dim[1]*model_in_dim[2]
        model_out_dim = train_out.shape[1:]
        model_out_dim = model_out_dim[0]

        model = Gross_neuron(model_in_dim, model_out_dim,
                             dimen=2, batch_size=batch_size)
        print(model_in_dim, model_out_dim)
        print('Start training.')
        history = model.fit(train_in, train_out, batch_size=batch_size,
                            epochs=epochs, verbose=1, validation_split=0.2)
        model.summary()
    else:
        pass


    #model.save(filepath)

    print("Training took {0} seconds.".format(time.time() - start_time))
    return model, history

def model_synthesize(model,data):
    #model = load_model(filepath)
    return model.predict(data,verbose=1)


def PlotHistory(train_value, test_value, value_is_loss_or_acc):
    f, ax = plt.subplots()
    ax.plot([None] + train_value, 'o-')
    ax.plot([None] + test_value, 'x-')
    ## Plot legend and use the best location automatically: loc = 0.
    ax.legend(['Train ' + value_is_loss_or_acc, 'Validation ' + value_is_loss_or_acc], loc = 0) 
    ax.set_title('Training/Validation ' + value_is_loss_or_acc + ' per Epoch')
    ax.set_xlabel('Epoch')
    ax.set_ylabel(value_is_loss_or_acc)

#PlotHistory(training_history.history['loss'],
#            training_history.history['val_loss'], 'Loss')
#PlotHistory(training_history.history['acc'],
#            training_history.history['val_acc'], 'Accuracy')

def drawWeightHistogram(model):
    x = model.layers[0].get_weights()[0].flatten()
    ## the histogram of the data
    fig = plt.subplots()
    n, bins, patches = plt.hist(x, 50)
    plt.xlim(-0.5, 0.5)
    plt.xlabel('Weight')
    plt.ylabel('Count')
    zero_counts = (x == 0.0).sum()
    plt.title("Weight Histogram. Num of '0's: %d" % zero_counts)

#w1 = trained_model.layers[0].get_weights()[0].flatten()
#drawWeightHistogram(w1)


def TestModel(model=None, data_in=None,data_out=None):
    if model is None:
        print("Must provide a trained model.")
        return
    if data_in is None:
        print("Must provide data.")
        return
    #x_test, y_test = data
    scores = model.evaluate(data_in, data_out)
    print(scores)
    print('-----------------------------')
    print("Test loss {:.4f}, accuracy {:.2f}%".format(scores[0], scores[1] * 100))
    print('-----------------------------')

    return scores

