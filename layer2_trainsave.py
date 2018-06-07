 #
# deep_learning_cross_validation_Reg_dropout.py
#
import numpy, argparse
from pandas import read_csv
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras.wrappers.scikit_learn import KerasClassifier
from keras.constraints import maxnorm
from keras.optimizers import SGD
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from itertools import cycle

from scipy import interp
from keras import regularizers
from sklearn.metrics import roc_curve, auc

import pandas as pd


import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

prog_name = "layer2_trainsave.py"

parser = argparse.ArgumentParser(description='Train and save', prog = prog_name)
parser.add_argument('-f1', '--training', required = True, metavar = 'The training dataset', type = str, help ='Upload the training dataset')
parser.add_argument('-epochs', '--nb_epoch',  default=100, metavar = 'Number of epoch', type = int, help ='Please select the number of epoch')
parser.add_argument('-batch', '--batch_size', default=200, metavar = 'Size of the batch', type = int, help ='The size of the batch to be train by model')
parser.add_argument('-d', '--data_dimension', default=5, metavar = 'Dimension of input data', type = int, help ='The dimension of input data.')
parser.add_argument('-pf', '--file_prefix', required = True, metavar = 'The file prefix of output', type = str, help ='The file prefix of output.')

args = parser.parse_args()

training = os.path.abspath(args.training)
nb_epoch = args.nb_epoch
batch_size = args.batch_size

seed = 7
numpy.random.seed(seed)

# load pima indians dataset

#dataset = numpy.loadtxt("Final_SZ_genes_dataset")
dataset = numpy.loadtxt(training)

X = dataset[:,0:args.data_dimension]
Y = dataset[:,args.data_dimension]

# define 10-fold cross validation test harness

model = Sequential()
model.add(Dropout(0.5, input_shape=(args.data_dimension,)))
model.add(Dense(12, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
model.add(Dense(10, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
model.add(Dense(1, kernel_initializer='normal', activation='sigmoid'))

sgd = SGD(lr=0.1, momentum=0.9, decay=0.0, nesterov=False)
#model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
print('epochs=%d batch_size=%d' % (nb_epoch, batch_size))
model.fit(X, Y, epochs=nb_epoch, batch_size=batch_size) #, verbose=0)
model.evaluate(X, Y)

model.save(args.file_prefix+"_model_f"+str(args.data_dimension)+"_layer2.h5")

print("NOTICE: Training is done! \n -------------------------------------------------------")









