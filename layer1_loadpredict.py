


import sys, argparse, os
seed = 7

prog_name = 'layer1_loadpredict.py'

def main():
    parser = argparse.ArgumentParser(description='Train and predict', prog = prog_name)
    parser.add_argument('-f2', '--testing', required = True, metavar = 'The testing dataset', type = str, help ='Upload the testing dataset')
    parser.add_argument('-d', '--data_dimension', required = True, metavar = 'Dimension of input data', type = int, help ='The dimension of input data.')
    parser.add_argument('-pf', '--file_prefix', required = True, metavar = 'The file prefix of output', type = str, help ='The file prefix of output.')

    args = parser.parse_args()
    
    testing = os.path.abspath(args.testing)

    print("\nNOTICE: Keras packages are importing! \n -------------------------------------------------------")
    
    from keras.layers import Dense, Dropout 
    from keras.models import Sequential
    from keras.layers import Dense
    from keras.models import load_model
    
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

    dataset_te= numpy.loadtxt(testing)
    
    input_dim=args.data_dimension
    X_te = dataset_te[:,0:input_dim]
    Y_te = dataset_te[:,input_dim]

    model = load_model(args.file_prefix+"_model_f"+str(input_dim)+"_layer1.h5")    

    print('\nNOTICE: Evaluating the model! \n -------------------------------------------------------')
    scores = model.evaluate(X_te, Y_te)
    print("\n %s: %.2f%%" % (model.metrics_names[1], scores[1]*100))

    print("NOTICE: Predicting on testing data! \n -------------------------------------------------------")
    Y_score = model.predict(X_te)
    numpy.savetxt(args.file_prefix+'_ncDeepBrain_score.txt',Y_score)
    rounded = [round(x[0]) for x in Y_score]

    fprl, tprl, thresholds = roc_curve(Y_te, Y_score)
    roc_auc = auc(fprl, tprl)
    print('AUC: %f \n -------------------------------------------------------' % roc_auc)
    
    print("NOTICE: Prediction is done! \n -------------------------------------------------------")

if __name__ == '__main__':
        main()

