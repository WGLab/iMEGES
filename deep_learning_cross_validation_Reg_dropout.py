 #MLP for Pima Indians Dataset with 10-fold cross validation
import numpy
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

seed = 7
numpy.random.seed(seed)

# load pima indians dataset


dataset = numpy.loadtxt("ASD_dataset_frequency_match.txt")
#dataset = numpy.loadtxt("DeepSea_GRASP_eQTLs_dataset.txt")

#dataset = numpy.loadtxt("DeepSea_training_GWAS_Catalog_dataset.txt")

#dataset = numpy.loadtxt("gwava_dataset_new_marks_v02_without_gnomad.txt")

#dataset = numpy.loadtxt("gwava_dataset_new_marks_v02.txt")
#dataset = numpy.loadtxt("cuarted_SNPs_dataset.txt")

#dataset = numpy.loadtxt("eQTLs_SNPs_dataset.txt")

#dataset = numpy.loadtxt("ASD_dataset.txt")
#dataset = numpy.loadtxt("ASD_dataset_withoutgenomad.txt")

#dataset = numpy.loadtxt("deltaSVM_dataset_new_marks_v02.txt")
#dataset = numpy.loadtxt("deltaSVM_dataset_new_marks.txt")
#dataset = numpy.loadtxt("deltaSVM_dataset.txt")
#dataset = numpy.loadtxt("Final_SCHI_dataset.txt")
#dataset = numpy.loadtxt("Final_SCHI_dataset_frequency_match.txt")
#dataset = numpy.loadtxt("PRVCS_refined_training_dataset.txt")
#PRVCS_refined_final_dataset_v01.txt
# split into input (X) and output (Y) variables

#X = dataset[:,0:214]
#Y = dataset[:,214]

X = dataset[:,0:11]
Y = dataset[:,11]

# define 10-fold cross validation test harness

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
model = Sequential()
model.add(Dropout(0.2, input_shape=(11,)))
model.add(Dense(50, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
model.add(Dense(40, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
model.add(Dense(30, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
#model.add(Dense(30, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
#model.add(Dense(20, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
model.add(Dense(10, kernel_initializer='normal', activation='relu', kernel_constraint=maxnorm(3)))
model.add(Dense(1, kernel_initializer='normal', activation='sigmoid'))
#model.add(Dense(50, input_dim=214, init='uniform', activation='relu',  kernel_regularizer=regularizers.l2(0.01)))
#model.add(Dense(40, init='uniform', activation='relu', kernel_regularizer=regularizers.l2(0.01)))
#model.add(Dense(30, init='uniform', activation='relu', kernel_regularizer=regularizers.l2(0.01)))
#model.add(Dense(20, init='uniform', activation='relu', kernel_regularizer=regularizers.l2(0.01)))
#model.add(Dense(10, init='uniform', activation='relu', kernel_regularizer=regularizers.l2(0.01)))
#model.add(Dense(1, init='uniform', activation='sigmoid'))

sgd = SGD(lr=0.1, momentum=0.9, decay=0.0, nesterov=False)
model.compile(loss='binary_crossentropy', optimizer=sgd, metrics=['accuracy'])
#model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])

mean_tpr = 0.0
mean_fpr = numpy.linspace(0, 1, 100)
all_tpr = []

colors = cycle(['cyan', 'indigo', 'seagreen', 'yellow', 'blue', 'darkorange', 'seagreen', 'yellow', 'blue', 'darkorange'])

lw = 2

cvscores = []

fig= plt.figure()
i=0
for (train, test), color in zip(cv.split(X, Y),colors):
        model.fit(X[train], Y[train], epochs=100, batch_size=200, verbose=0)
	scores = model.evaluate(X[test], Y[test])
	print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))
	cvscores.append(scores[1] * 100)
        Y_score = model.predict(X[test])
        numpy.savetxt('test.txt',Y_score)
        rounded = [round(x[0]) for x in Y_score]
        print(rounded)
        fpr, tpr, _ = roc_curve(numpy.round(Y[test]), Y_score)      
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = auc(fpr, tpr)
        

	plt.plot(fpr, tpr, lw=lw, color=color, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))
        i += 1
	print('AUC: %f \n -------------------------------------------------------' % roc_auc)
plt.plot([0, 1], [0, 1], linestyle='--', lw=lw, color='k', label='Luck')
mean_tpr /= cv.get_n_splits(X, Y)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
plt.plot(mean_fpr, mean_tpr, color='g', linestyle='--', label='Mean ROC (area = %0.2f)' % mean_auc, lw=lw)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic (ROC) curve')
plt.legend(loc="lower right")
plt.show()
fig.savefig('cross_validation_score_eQTLs.png')


print("%.2f%% (+/- %.2f%%)" % (numpy.mean(cvscores), numpy.std(cvscores)))












