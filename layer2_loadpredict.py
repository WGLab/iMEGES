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
from keras.models import load_model


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

prog_name = "layer2_loadpredict.py"

parser = argparse.ArgumentParser(description='Load and predict', prog = prog_name)
parser.add_argument('-f1', '--training', required = True, metavar = 'The testing dataset', type = str, help ='Upload the testing dataset')
parser.add_argument('-d', '--data_dimension', default=5, metavar = 'Dimension of input data', type = int, help ='The dimension of input data.')
parser.add_argument('-pf', '--file_prefix', required = True, metavar = 'The file prefix of output', type = str, help ='The file prefix of output.')

args = parser.parse_args()

training = os.path.abspath(args.training)

seed = 7
numpy.random.seed(seed)


#dataset = numpy.loadtxt("Final_SZ_genes_dataset")
dataset = numpy.loadtxt(training)

X = dataset[:,0:args.data_dimension]
Y = dataset[:,args.data_dimension]

# define 10-fold cross validation test harness
if True:
    model = load_model(args.file_prefix+"_model_f"+str(args.data_dimension)+"_layer2.h5")

    print('\nNOTICE: Evaluating the model! \n -------------------------------------------------------')
    scores = model.evaluate(X, Y)
    print("\n %s: %.2f%%" % (model.metrics_names[1], scores[1]*100))

    print("NOTICE: Predicting on testing data! \n -------------------------------------------------------")
    Y_score = model.predict(X)
    numpy.savetxt(args.file_prefix+'_ncGene_score.txt',Y_score)
    rounded = [round(x[0]) for x in Y_score]

    fprl, tprl, thresholds = roc_curve(Y, Y_score)
    roc_auc = auc(fprl, tprl)
    print('AUC: %f \n -------------------------------------------------------' % roc_auc)

    print("NOTICE: Prediction is done! \n -------------------------------------------------------")










