import sys
import os
import csv
import string

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("dataset", type=str, choices=['RECO', 'DIGI', 'SIM'], help="Specify which list")
args = parser.parse_args()


if args.dataset == 'SIM':
      listName = 'CRAB3_DsPhiPi_13TeV_SIM_190212_140348'
      path = '/lustre/cms/store/user/fsimone/DsPhiPi/CRAB3_DsPhiPi_13TeV_SIM/190212_140348/'
if args.dataset == 'DIGI':
      listName = 'crab_crab_DsPhiPi_13TeV_DIGI_190214_100253'
      path = '/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_DIGI/190214_100253/'
if args.dataset == 'RECO':
      listName = 'crab_crab_DsPhiPi_13TeV_RECO_190225_141742'
      path = '/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/'

#generating the list of all .root files in given directory and subdirectories
files_in_folder = []
for r, d, f in os.walk(path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            files_in_folder.append(os.path.join(r, file))
print('Total number of files in folder is {0:2d}'.format(len(files_in_folder)))
#prepare final script
good_fileList = open("good_fileList_"+args.dataset+".txt", "w")

for line in files_in_folder:
    found = 0
    with open(listName) as search:
        if line in search.read():
           found = 1
    if found == 0:
        good_fileList.write(line+'\n')

good_fileList.close()
