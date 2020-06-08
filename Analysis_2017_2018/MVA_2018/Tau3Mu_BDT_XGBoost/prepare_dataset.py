#!/bin/env python2
import sys
import os
import argparse

import pickle
import root_pandas
import numpy            as np
import pandas           as pd
import matplotlib.cm    as cm
import root_numpy

from sklearn.model_selection    import train_test_split
from collections                import OrderedDict

def common_sel(df):
    df_cut = df.loc[(df['fv_nC'] > 0) & (df['fv_nC'] < 100) & (df['TreeMu1.mu_trackerLayersWithMeasurement'] > 8) & (df['TreeMu2.mu_trackerLayersWithMeasurement'] > 8) & (df['TreeMu3.mu_trackerLayersWithMeasurement'] > 8) & (df['TreeMu1.mu_innerTrack_validFraction'] > 0.5) & (df['TreeMu2.mu_innerTrack_validFraction'] > 0.5) & (df['TreeMu3.mu_innerTrack_validFraction'] > 0.5) & (df['TreeMu1.mu_innerTrack_normalizedChi2'] < 40) & (df['TreeMu2.mu_innerTrack_normalizedChi2'] < 40) & (df['TreeMu3.mu_innerTrack_normalizedChi2'] < 40) & (df['TreeMu1.mu_combinedQuality_trkKink'] < 900) & (df['TreeMu2.mu_combinedQuality_trkKink'] < 900) & (df['TreeMu3.mu_combinedQuality_trkKink'] < 900)]
    return df_cut

def signal_sel(df):
    df_cut = df.loc[(df['tripletMass'] > 1.75) & (df['tripletMass'] < 1.80)]
    return df_cut
def bkg_sel(df):
    df_cut = df.loc[(df['tripletMass'] > 1.60) & (df['tripletMass'] < 2.02)]
    return df_cut

@np.vectorize
def cat(tripletMassReso):
    if   (tripletMassReso < 0.007 )                               : return "A"
    elif (tripletMassReso >= 0.007 and tripletMassReso  < 0.0105) : return "B"
    elif (tripletMassReso  >= 0.0105)                             : return "C"
    else                                                          : return NaN

@np.vectorize
def mass(tripletMass):
    if   (tripletMass > 1.75 and tripletMass < 1.80)             : return "tau"
    else                                                         : return "SB"

#labels = OrderedDict()
#
#labels['cLP'                            ] = '$\\tau$ $p_{T}$'
#labels['tKink'                             ] = '$m_{T}(\\tau, MET)$'

def gini(actual, pred, cmpcol = 0, sortcol = 1):
    assert( len(actual) == len(pred) )
    all = np.asarray(np.c_[ actual, pred, np.arange(len(actual)) ], dtype=np.float)
    all = all[ np.lexsort((all[:,2], -1*all[:,1])) ]
    totalLosses = all[:,0].sum()
    giniSum = all[:,0].cumsum().sum() / totalLosses
    
    giniSum -= (len(actual) + 1) / 2.
    return giniSum / len(actual)
 
def gini_normalized(a, p):
    return gini(a, p) / gini(a, a)

def gini_xgb(preds, dtrain):
    #labels = dtrain.get_label()
    gini_score = gini_normalized(labels, preds)
    return 'gini', gini_score

#features for training
features = []
feat_labels   = []
with open('/lustrehome/fsimone/MVA_2018/BDTinputVar_dataset_2018_6may_optimised_v2_A.txt', 'rb') as inputvar_file:
    for line in inputvar_file:
        fields = line.split('\t')
        features.append(fields[0])    #definition
        feat_labels.append(fields[0]) #names

train_features = [ff for ff in features] ## import purpose


signalDs_2018 = [
    '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1715/AnalysedTree_MC_2018Ds_tau3mu_6may.root',   ## MiniTree DsTau3Mu
]

signalB0_2018 = [
    '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1716/AnalysedTree_MC_2018B0_tau3mu_6may.root',
]

signalBp_2018 = [
    '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1716/AnalysedTree_MC_2018Bp_tau3mu_6may.root',
]

background_2018 = [
           '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200512_1532/AnalysedTree_data_2018A_tau3mu_12may.root',
           '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1713/AnalysedTree_data_2018B_tau3mu_6may.root',
           '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1713/AnalysedTree_data_2018C_tau3mu_6may.root',
           '/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1713/AnalysedTree_data_2018D_tau3mu_6may.root'
]

#TTree enriched with muonID value
signalDs_2018_muID = []
for treename in signalDs_2018:
   signalDs_2018_muID.append(treename.replace('.root', '_MuonID.root'))
signalB0_2018_muID = []
for treename in signalB0_2018:
   signalB0_2018_muID.append(treename.replace('.root', '_MuonID.root'))
signalBp_2018_muID = []
for treename in signalBp_2018:
   signalBp_2018_muID.append(treename.replace('.root', '_MuonID.root'))
background_2018_muID = []
for treename in background_2018:
   background_2018_muID.append(treename.replace('.root', '_MuonID.root'))


#branches to be taken from "main" minitree
main_branches = [
    'run',
    'evt',
    'lumi',
    'puFactor',
    'tripletMass',
    'tripletMassReso',
    'P_tripl',
    'Pt_tripl',
    'Eta_tripl',
    'Pmu3',
    'cLP',
    'tKink',
    'segmComp',
    'fv_nC',
    'fv_dphi3D',
    'fv_d3Dsig',
    'fv_d3D',
    'bs_sv_d3Dsig',
    'bs_sv_d3D',
    'd0sig',
    'd0',
    'mindca_iso',
    'trkRel',
    'RefVx1',
    'RefVy1',
    'RefVz1',
    'RefVx2',
    'RefVy2',
    'RefVz2',
    'RefVx3',
    'RefVy3',
    'RefVz3',
    'SVx',
    'SVy',
    'SVz',
    'dxy1',
    'dxy2',
    'dxy3',
    'dxyErr1',
    'dxyErr2',
    'dxyErr3',
]

#main ntuples
sigDs_2018_main = pd.DataFrame( root_numpy.root2array(signalDs_2018  , 'FinalTree', branches = main_branches) )
sigB0_2018_main = pd.DataFrame( root_numpy.root2array(signalB0_2018  , 'FinalTree', branches = main_branches) )
sigBp_2018_main = pd.DataFrame( root_numpy.root2array(signalBp_2018  , 'FinalTree', branches = main_branches) )
bkg_2018_main  =  pd.DataFrame( root_numpy.root2array(background_2018, 'FinalTree', branches = main_branches) )

#single muon info
sigDs_2018_mu1 = pd.DataFrame( root_numpy.root2array(signalDs_2018   , 'TreeMu1') )
sigDs_2018_mu2 = pd.DataFrame( root_numpy.root2array(signalDs_2018   , 'TreeMu2') )
sigDs_2018_mu3 = pd.DataFrame( root_numpy.root2array(signalDs_2018   , 'TreeMu3') )
sigDs_2018_mu1.columns = ['TreeMu1.'+x for x in sigDs_2018_mu1.columns]
sigDs_2018_mu2.columns = ['TreeMu2.'+x for x in sigDs_2018_mu2.columns]
sigDs_2018_mu3.columns = ['TreeMu3.'+x for x in sigDs_2018_mu3.columns]

sigB0_2018_mu1 = pd.DataFrame( root_numpy.root2array(signalB0_2018   , 'TreeMu1') )
sigB0_2018_mu2 = pd.DataFrame( root_numpy.root2array(signalB0_2018   , 'TreeMu2') )
sigB0_2018_mu3 = pd.DataFrame( root_numpy.root2array(signalB0_2018   , 'TreeMu3') )
sigB0_2018_mu1.columns = ['TreeMu1.'+x for x in sigB0_2018_mu1.columns]
sigB0_2018_mu2.columns = ['TreeMu2.'+x for x in sigB0_2018_mu2.columns]
sigB0_2018_mu3.columns = ['TreeMu3.'+x for x in sigB0_2018_mu3.columns]

sigBp_2018_mu1 = pd.DataFrame( root_numpy.root2array(signalBp_2018   , 'TreeMu1') )
sigBp_2018_mu2 = pd.DataFrame( root_numpy.root2array(signalBp_2018   , 'TreeMu2') )
sigBp_2018_mu3 = pd.DataFrame( root_numpy.root2array(signalBp_2018   , 'TreeMu3') )
sigBp_2018_mu1.columns = ['TreeMu1.'+x for x in sigBp_2018_mu1.columns]
sigBp_2018_mu2.columns = ['TreeMu2.'+x for x in sigBp_2018_mu2.columns]
sigBp_2018_mu3.columns = ['TreeMu3.'+x for x in sigBp_2018_mu3.columns]

bkg_2018_mu1  =  pd.DataFrame( root_numpy.root2array(background_2018,  'TreeMu1') )
bkg_2018_mu2  =  pd.DataFrame( root_numpy.root2array(background_2018,  'TreeMu2') )
bkg_2018_mu3  =  pd.DataFrame( root_numpy.root2array(background_2018,  'TreeMu3') )
bkg_2018_mu1.columns = ['TreeMu1.'+x for x in bkg_2018_mu1.columns]
bkg_2018_mu2.columns = ['TreeMu2.'+x for x in bkg_2018_mu2.columns]
bkg_2018_mu3.columns = ['TreeMu3.'+x for x in bkg_2018_mu3.columns]

#MuonID evaluation
sigDs_2018_muonid1 = pd.DataFrame( root_numpy.root2array(signalDs_2018_muID  , 'MuonIDeval_Mu1') )
sigDs_2018_muonid2 = pd.DataFrame( root_numpy.root2array(signalDs_2018_muID  , 'MuonIDeval_Mu2') )
sigDs_2018_muonid3 = pd.DataFrame( root_numpy.root2array(signalDs_2018_muID  , 'MuonIDeval_Mu3') )
sigDs_2018_muonid1.rename(columns={'MuonID':'MuonIDeval_Mu1.MuonID'}, inplace=True)
sigDs_2018_muonid2.rename(columns={'MuonID':'MuonIDeval_Mu2.MuonID'}, inplace=True)
sigDs_2018_muonid3.rename(columns={'MuonID':'MuonIDeval_Mu3.MuonID'}, inplace=True)

sigB0_2018_muonid1 = pd.DataFrame( root_numpy.root2array(signalB0_2018_muID  , 'MuonIDeval_Mu1') )
sigB0_2018_muonid2 = pd.DataFrame( root_numpy.root2array(signalB0_2018_muID  , 'MuonIDeval_Mu2') )
sigB0_2018_muonid3 = pd.DataFrame( root_numpy.root2array(signalB0_2018_muID  , 'MuonIDeval_Mu3') )
sigB0_2018_muonid1.rename(columns={'MuonID':'MuonIDeval_Mu1.MuonID'}, inplace=True)
sigB0_2018_muonid2.rename(columns={'MuonID':'MuonIDeval_Mu2.MuonID'}, inplace=True)
sigB0_2018_muonid3.rename(columns={'MuonID':'MuonIDeval_Mu3.MuonID'}, inplace=True)

sigBp_2018_muonid1 = pd.DataFrame( root_numpy.root2array(signalBp_2018_muID  , 'MuonIDeval_Mu1') )
sigBp_2018_muonid2 = pd.DataFrame( root_numpy.root2array(signalBp_2018_muID  , 'MuonIDeval_Mu2') )
sigBp_2018_muonid3 = pd.DataFrame( root_numpy.root2array(signalBp_2018_muID  , 'MuonIDeval_Mu3') )
sigBp_2018_muonid1.rename(columns={'MuonID':'MuonIDeval_Mu1.MuonID'}, inplace=True)
sigBp_2018_muonid2.rename(columns={'MuonID':'MuonIDeval_Mu2.MuonID'}, inplace=True)
sigBp_2018_muonid3.rename(columns={'MuonID':'MuonIDeval_Mu3.MuonID'}, inplace=True)

bkg_2018_muonid1 = pd.DataFrame( root_numpy.root2array(background_2018_muID  , 'MuonIDeval_Mu1') )
bkg_2018_muonid2 = pd.DataFrame( root_numpy.root2array(background_2018_muID  , 'MuonIDeval_Mu2') )
bkg_2018_muonid3 = pd.DataFrame( root_numpy.root2array(background_2018_muID  , 'MuonIDeval_Mu3') )
bkg_2018_muonid1.rename(columns={'MuonID':'MuonIDeval_Mu1.MuonID'}, inplace=True)
bkg_2018_muonid2.rename(columns={'MuonID':'MuonIDeval_Mu2.MuonID'}, inplace=True)
bkg_2018_muonid3.rename(columns={'MuonID':'MuonIDeval_Mu3.MuonID'}, inplace=True)

#All infos in dataframe
sigDs_2018 = pd.concat([sigDs_2018_main, sigDs_2018_mu1, sigDs_2018_mu2, sigDs_2018_mu3, sigDs_2018_muonid1, sigDs_2018_muonid2, sigDs_2018_muonid3], axis=1, sort=False)
sigB0_2018 = pd.concat([sigB0_2018_main, sigB0_2018_mu1, sigB0_2018_mu2, sigB0_2018_mu3, sigB0_2018_muonid1, sigB0_2018_muonid2, sigB0_2018_muonid3], axis=1, sort=False)
sigBp_2018 = pd.concat([sigBp_2018_main, sigBp_2018_mu1, sigBp_2018_mu2, sigBp_2018_mu3, sigBp_2018_muonid1, sigBp_2018_muonid2, sigBp_2018_muonid3], axis=1, sort=False)
bkg_2018   = pd.concat([bkg_2018_main, bkg_2018_mu1, bkg_2018_mu2, bkg_2018_mu3, bkg_2018_muonid1, bkg_2018_muonid2, bkg_2018_muonid3], axis=1, sort=False)

#print(sigDs_2018_main.head())
#print(sigDs_2018_mu1.head() )
#print(sigDs_2018_mu2.head() )
#print(sigDs_2018_mu3.head() )
#print(sigDs_2018_muonid1.head() )
#print(sigDs_2018_muonid2.head() )
#print(sigDs_2018_muonid3.head() )

#clean memory from partial dataframes
del sigDs_2018_main
del sigB0_2018_main
del sigBp_2018_main
del bkg_2018_main 
del sigDs_2018_mu1 
del sigDs_2018_mu2 
del sigDs_2018_mu3 
del sigB0_2018_mu1 
del sigB0_2018_mu2 
del sigB0_2018_mu3 
del sigBp_2018_mu1 
del sigBp_2018_mu2 
del sigBp_2018_mu3 
del bkg_2018_mu1
del bkg_2018_mu2
del bkg_2018_mu3
del sigDs_2018_muonid1
del sigDs_2018_muonid2
del sigDs_2018_muonid3
del sigB0_2018_muonid1
del sigB0_2018_muonid2
del sigB0_2018_muonid3
del sigBp_2018_muonid1
del sigBp_2018_muonid2
del sigBp_2018_muonid3
del bkg_2018_muonid1
del bkg_2018_muonid2
del bkg_2018_muonid3

#useful for future combination
features.append('year')
sigDs_2018['year'] = np.full(sigDs_2018.shape[0], 2018)
sigB0_2018['year'] = np.full(sigB0_2018.shape[0], 2018)
sigBp_2018['year'] = np.full(sigBp_2018.shape[0], 2018)
bkg_2018  ['year'] = np.full(bkg_2018.shape [0], 2018)

#scale factors for signal normalization
features.append('norm')
sigDs_2018['norm'] = np.full(sigDs_2018.shape[0], 1.323E-03)
sigB0_2018['norm'] = np.full(sigB0_2018.shape[0], 4.784E-04)
sigBp_2018['norm'] = np.full(sigBp_2018.shape[0], 1.437E-03)
bkg_2018  ['norm'] = np.full(bkg_2018.shape [0],  1.0)

#weights to be used during training will reflect signals normalisation
features.append('weight')
sigDs_2018['weight'] = np.full(sigDs_2018.shape[0], 1.323E-03/1.323E-03)
sigB0_2018['weight'] = np.full(sigB0_2018.shape[0], 4.784E-04/1.323E-03)
sigBp_2018['weight'] = np.full(sigBp_2018.shape[0], 1.437E-03/1.323E-03)
bkg_2018  ['weight'] = np.full(bkg_2018.shape [0],  1.0)

#concatenate signals and background(s)
sig = pd.concat([sigDs_2018, sigB0_2018, sigBp_2018], ignore_index = True)
bkg = bkg_2018 #bkg  = pd.concat([bkg_2017, bkg_2018 ], ignore_index = True)

#clean memory
del sigDs_2018
del sigB0_2018
del sigBp_2018
del bkg_2018  

## correctly set the target shape as (#events, #classes)
sig['target'] =  np.ones (sig.shape[0]).astype(np.int)
bkg['target'] =  np.zeros(bkg.shape[0]).astype(np.int)

#append category
for dd in [bkg, sig]:
    dd['category'] = cat(dd['tripletMassReso'])

#append mass range
for dd in [bkg, sig]:
    dd['massrange'] = mass(dd['tripletMass'])

#append derived variables
for dd in [bkg, sig]:
    dd['abs(dxy3/dxyErr3)'] = abs(dd['dxy3']/dd['dxyErr3'])
    dd['PS_SV_dxy'        ] = np.sqrt( (dd['RefVx1']-dd['SVx']) * (dd['RefVx1']-dd['SVx']) + (dd['RefVy1']-dd['SVy']) * (dd['RefVy1']-dd['SVy']) )
    dd['PS_SV_dz'         ] = abs(dd['RefVz1']-dd['SVz'])


#apply selections
print("common preselections")
print(sig.shape)
sig = common_sel(sig)
print(bkg.shape)
bkg = common_sel(bkg)
