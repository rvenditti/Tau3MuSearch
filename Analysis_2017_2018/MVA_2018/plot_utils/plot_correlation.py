#system
import sys, os
import re

#ROOT
import ROOT
from root_pandas import read_root

#data visualization
import seaborn as sn
import matplotlib.pyplot as plt

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("inputFile", type=str, help="Specify path of TMVA rootfile")
parser.add_argument("category", type=str, help="Specify category")
args = parser.parse_args()

datasetName = args.inputFile.replace('.root','')
datasetName = datasetName.split('TMVA_',1)[1]
print(datasetName)
df = read_root(args.inputFile, '/'+datasetName+'/'+"TestTree")

df_signal = df[df.classID==0]
df_bkg = df[df.classID==1]
print(df_signal)
print(df_bkg)
drop_columns = [
                'classID',
                'className', 
#                'Pmu3_45_45:Pmu3_Pmu3',
                'cLP_30_30:cLP_cLP',
                'tKink_80_80:tKink_tKink',
                'fv_nC_25_25:fv_nC_fv_nC',
                'fv_dphi3D_0.15_0.15:fv_dphi3D_fv_dphi3D',
                'fv_d3Dsig_100_100:fv_d3Dsig_fv_d3Dsig',
                'segmComp_0.2_0.2:segmComp_segmComp',
                'd0sig_15_15:d0sig_d0sig',
                'mindca_iso_0.5_0.5:mindca_iso_mindca_iso',
                'trkRel_10_10:trkRel_trkRel', 
                'weight', 
                'puFactor', 
                'evt', 
                'tripletMassReso', 
                'BDT'
                ]
df_signal = df_signal.drop(drop_columns, axis=1)
df_bkg = df_bkg.drop(drop_columns, axis=1)
print(df_signal)
print(df_bkg)

#print(df_bkg)
rename_columns = [
 #                 "Pmu3",
                  "cLP",
                  "tKink", 
                  "fv_nC",
                  "fv_dphi3D",
                  "d3Dsig",
                  "d0sig",
                  "mindca_iso",
                  "trkRel", 
                  "segmComp",
                  "MuID1",
                  "MuID2",
                  "MuID3",
                  "Pt_tripl",
                  "m3m",
                  ]
df_signal.columns = rename_columns
df_bkg.columns = rename_columns
#print(df_bkg)
corrMatrix_signal = df_signal.corr()
corrMatrix_bkg = df_bkg.corr()

#print(corrMatrix_signal)

g_signal = sn.heatmap(
        corrMatrix_signal,
        square=True,
        cbar_kws={'shrink':.9 },
        annot=True,
        linewidths=0.1,vmax=1.0, linecolor='white',
        annot_kws={'fontsize':8 },
        fmt=".1%"
    )
g_signal.set_yticklabels(g_signal.get_yticklabels(), rotation = 0, fontsize = 8)
g_signal.set_xticklabels(g_signal.get_xticklabels(), rotation = 45, fontsize = 8)

plt.title('Pearson Correlation of Features (Signal) - cat. '+args.category, y=1.05, size=20)
plt.savefig(datasetName+'/correlationSignal.png', verbose=True)
plt.close()

g_bkg = sn.heatmap(
        corrMatrix_bkg,
        square=True,
        cbar_kws={'shrink':.9 },
        annot=True,
        linewidths=0.1,vmax=1.0, linecolor='white',
        annot_kws={'fontsize':8 },
        fmt=".1%"
    )
g_bkg.set_yticklabels(g_bkg.get_yticklabels(), rotation = 0, fontsize = 8)
g_bkg.set_xticklabels(g_bkg.get_xticklabels(), rotation = 45, fontsize = 8)

plt.title('Pearson Correlation of Features (Background) - cat. '+args.category, y=1.05, size=20)
plt.savefig(datasetName+'/correlationBkg.png', verbose=True)
plt.close()
