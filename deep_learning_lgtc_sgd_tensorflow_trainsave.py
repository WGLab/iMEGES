


import sys, argparse, os
seed = 7

prog_name = 'deep_learning_lgtc_sgd_tensorflow_pred.py'

def main():
    parser = argparse.ArgumentParser(description='Train and predict', prog = prog_name)
    parser.add_argument('-f1', '--training', required = True, metavar = 'The training dataset', type = str, help ='Upload the training dataset')
    parser.add_argument('-epochs', '--nb_epoch', required = True, metavar = 'Number of epoch', type = int, help ='Please select the number of epoch')
    parser.add_argument('-batch', '--batch_size', required = True, metavar = 'Size of the batch', type = int, help ='The size of the batch to be train by model')
    parser.add_argument('-d', '--data_dimension', required = True, metavar = 'Dimension of input data', type = int, help ='The dimension of input data.')
    parser.add_argument('-pf', '--file_prefix', required = True, metavar = 'The file prefix of output', type = str, help ='The file prefix of output.')

    args = parser.parse_args()
    
    training = os.path.abspath(args.training)
    nb_epoch = args.nb_epoch
    batch_size = args.batch_size

    print("\nNOTICE: Keras packages are importing! \n -------------------------------------------------------")
    
    from keras.layers import Dense, Dropout 
    from keras.models import Sequential
    from keras.layers import Dense
    
    from argparse import ArgumentParser
    from sklearn.utils import shuffle
    from sklearn.linear_model import SGDClassifier
    from sklearn.preprocessing import StandardScaler
    import scipy
    from sklearn.metrics import f1_score
    import numpy
    from sklearn.metrics import roc_curve, auc
    import csv
    import sys
    import matplotlib as mpl
    if os.environ.get('DISPLAY','') == '':
        print('NOTICE: no display found. Using non-interactive Agg backend print \n -------------------------------------------------------')
        mpl.use('Agg')
    import matplotlib.pyplot as plt

    print(" ------------------------------------------------------- \nNOTICE: The data is loading! \n -------------------------------------------------------")

    dataset_tr = numpy.loadtxt(training)
    input_dim=args.data_dimension
    X_tr = dataset_tr[:,0:input_dim]
    Y_tr = dataset_tr[:,input_dim]

    print('NOTICE: Creating the model! \n -------------------------------------------------------')
    model = Sequential()
    model.add(Dense(12, input_dim=input_dim, kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(8, kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(15, kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(12, kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(10, kernel_initializer='uniform', activation='relu'))
    model.add(Dropout(0.2))
    model.add(Dense(1, kernel_initializer='uniform', activation='sigmoid'))

    print('NOTICE: Compiling the model! \n -------------------------------------------------------')
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
    
    print('NOTICE: Training the model! \n -------------------------------------------------------')
    model.fit(X_tr, Y_tr, epochs=nb_epoch, batch_size=batch_size)

    model.save(args.file_prefix+"_model_f"+str(input_dim)+".h5")
     
    print("NOTICE: Training is done! \n -------------------------------------------------------")

if __name__ == '__main__':
        main()

