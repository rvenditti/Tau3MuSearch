#!/bin/env python2
import sys
import os
import argparse
from argparse import RawTextHelpFormatter

import pickle
import root_pandas
import numpy            as np
import pandas           as pd
import matplotlib.cm    as cm

from sklearn            import preprocessing
from sklearn.externals  import joblib
from sklearn.model_selection    import train_test_split

sys.path.append('/lustrehome/fsimone/MVA_2018/Tau3Mu_BDT_XGBoost')
from prepare_dataset  import sig, bkg, train_features, feat_labels

from trainer    import start_XGBoost_trainer, add_bdt_score
from plotter    import plot_overtraining, plot_ROC, plot_correlation_matrix, plot_features

#argument parser
parser = argparse.ArgumentParser('XGBoost training script')
parser.add_argument('--categ'    , default = 0           , help = 'event category')
parser.add_argument('--load'     , default = None        , help = 'load an existing training from the specified file')
parser.add_argument('--save_tree', action  = 'store_true', help = 'save enriched ntuples after training or loading'  )
args = parser.parse_args()

category    = args.categ
#tag to be used for output plots and files
tag = os.popen('date +%s').readlines()[0].strip('\n')
tag = tag+"_"+category

## assign an id to the test and train sets seprately to avoid mismatch when folding
sig.insert(len(sig.columns), 'id', np.arange(len(sig)))
bkg.insert(len(bkg.columns), 'id', np.arange(len(bkg)))

#prepare train and test dataframes
print("category "+category)
if category is not 0:
        sig = sig.loc[(sig['category'] == category)]
        bkg  = bkg.loc[(bkg['category'] == category)]
print("massrange for signal is [1.75;1.80], for background use SB")
sig_train  = sig.loc[(sig['massrange'] == "tau")]
bkg_train  = bkg.loc[(bkg['massrange'] == "SB")]

data = pd.concat([sig_train, bkg_train], ignore_index = True)
#data['id'] = np.arange(len(data))

train, test = train_test_split(data, test_size=0.4, random_state=1986)
## assign an id to the test and train sets seprately to avoid mismatch when folding
#train.insert(len(train.columns), 'id', np.arange(len(train)))
#test .insert(len(test .columns), 'id', np.arange(len(test )))

## train a new classifier or load an existing one from a .pck file
if args.load is None:
    classifiers = start_XGBoost_trainer(
        train       = train     ,
        test        = test      ,
        tag         = tag       ,
        features    = train_features  ,
        category    = args.categ,
    )

    if not os.path.exists('./pck'): os.mkdir('./pck')

    classifier_file = open('./pck/classifier_%s.pck' %tag, 'w+')
    pickle.dump(classifiers, classifier_file)
    classifier_file.close()
else: 
    classifiers = joblib.load(args.load)

## add the BDT score to the training and test dataset
for jj, clf in enumerate(classifiers):
    add_bdt_score(classifier = clf, sample = train, features = train_features, scale = len(classifiers), index = jj)
    add_bdt_score(classifier = clf, sample = test , features = train_features, scale = len(classifiers), index = jj)

## plot ROCs, correlation matrices and overtraining tests
if not os.path.exists('./pdf'):         os.mkdir('./pdf')
if not os.path.exists('./pdf/%s' %tag): os.mkdir('./pdf/%s' %tag)

plot_ROC(y = test ['target'], score = test ['bdt'], title = 'ROC curve (test)' , filename = './pdf/%s/roc_test.pdf'  %tag, color = 'r', label = 'BDT score', xlab = 'FPR', ylab = 'TPR')
plot_ROC(y = train['target'], score = train['bdt'], title = 'ROC curve (train)', filename = './pdf/%s/roc_train.pdf' %tag, color = 'b', label = 'BDT score', xlab = 'FPR', ylab = 'TPR', save_file = True)

plot_overtraining(  train  = train, 
                    test   = test , 
                    target = 1, score = 'bdt' ,
                    title  = 'bkg proba distribution' , filename = './pdf/%s/overtraining_bkg_proba.pdf' %tag)

plot_correlation_matrix( sample = pd.concat([train[train.target == 0], test[test.target == 0]]), features = train_features + ['bdt', 'tripletMass'], labels = feat_labels + ['bdt', 'tripletMass'], label = '', filename = './pdf/%s/corr_mat_bkg.pdf' %tag)
plot_correlation_matrix( sample = pd.concat([train[train.target == 1], test[test.target == 1]]), features = train_features + ['bdt', 'tripletMass'], labels = feat_labels + ['bdt', 'tripletMass'], label = '', filename = './pdf/%s/corr_mat_signal.pdf' %tag)
plot_correlation_matrix( sample = pd.concat([train, test]), features = train_features + ['bdt', 'tripletMass'], labels = feat_labels + ['bdt', 'tripletMass'], label = '', filename = './pdf/%s/corr_mat_bkg.pdf' %tag)

#plot_features(classifiers = classifiers, labels = feat_labels, filename = './pdf/%s/f_score.pdf' %tag)


### ON FULL MASS RANGE DATASETS

## add the BDT score to the original dataset
for jj, clf in enumerate(classifiers):
    add_bdt_score(classifier = clf, sample = sig, features = train_features, scale = len(classifiers), index = jj)
    add_bdt_score(classifier = clf, sample = bkg , features = train_features, scale = len(classifiers), index = jj)

## save the enriched ntuples
if args.save_tree:

    if not os.path.exists('./ntuples'): os.mkdir('./ntuples')

    #os.remove('./ntuples/signal_{TAG}.root'    .format(TAG = tag))
    #os.remove('./ntuples/background_{TAG}.root'.format(TAG = tag))

    sig.to_root('./ntuples/signal_{TAG}.root'    .format(TAG = tag), key = 'tree')
    bkg.to_root('./ntuples/background_{TAG}.root'.format(TAG = tag), key = 'tree')

print 'all done'
