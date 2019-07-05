# Tau3MuSearch
```bash
mkdir Tau3Mu
cd Tau3Mu
cmsrel CMSSW_9_4_4
cd CMSSW_9_4_4/src
cmsenv
git cms-init
git clone https://github.com/rvenditti/Tau3MuSearch
```
# Content

- MiniAna2017Tree module that makes trees storing variables related to the single reco muons and the triplet candidates

- Tau3MuSkim_cff.py is an example of configuration file selecting the triplet candidates (e.g. 3 loose muons with same vertex)

- run_Tau3MuSkim_cfg.py is the file that you have to "cmsRun" to generate the tree

- "Analysis" folder contains some tools useful for the analysis

In order to perform analysis:
* Compile the macro **Analysis.cpp**
* Execute it giving as input the 2 proper arguments according to the analysis you want to perform. In particular :
* **1st input arg** :  _Type_ of analysis you want to perform; the possibilities are :
* "MC" : Analysis DsTau3Mu on MC datasets 
* "MC_sgn" : Analysis DsTau3Mu on MC datasets selecting only the signal region
* "MC_control" :  Analysis DsPhiPi only on MC DsPhiPi dataset
* "data" : Analysis DsTau3Mu on data 
* "data_bkg" : Analysis DsTau3Mu on data selecting only the sideband (background) region
* "data_control" :   Analysis DsPhiPi on data 
* **2nd input arg** :  _Dataset_ on which the analysis has to be performed; the possibilities are :
* "Ds", "B0", "Bp", "Minibias" : for the "MC" or "MC_sgn" type of analysis
* "2017B", "2017C", ... : for the "data", "data_bkg", "data_control" type of analysis


Example:
In order to perform the _DsPhiPi_ analysis on _2017B_ dataset :
* Compile the Analysis macro : 
```
g++ -I $ROOTSYS/include -I $ROOTSYS/include Analysis.cpp `root-config --glibs` `root-config --libs` `root-config --cflags`  -L $ROOTSYS/lib -o AnalysisCode
```
* Execute :
```
./AnalysisCode "data_control" "2017B"
```

In order to perform the same analysis on RECAS (computer farm in Bari) using the batch system, it's necessary to create 2 .job files, and then, to submit the job:
```
condor_submit my_HTCondor_Control_data_2017B.job -name ettore
```

    
