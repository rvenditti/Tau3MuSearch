# MVA training and evaluation

Set of macros for BDT training, BDT score evaluation on data and MC, categorisation and mass distribution production


### Installing
```
cd $WORKDIR/CMSSW_10_6_8/src
cmsenv
```

## Example

BDT training
- set input files in `T3M_common.h` together with BDT input variables and spectators
- set BDT folder in `T3M_common.h` > 'TMVA_outputpath='
- `source run_MVA.sh` to launch BDT training. The weights will be produced in folders 'TMVA_outputpath'+category 

BDT evaluation
- set BDT folder in `T3M_common.h` > 'TMVA_inputpath=' for the desired weights. List of variables should match the one used for training
- 'root -l evaluation_fillbdt.C' evaluates the BDT running on input files and produces BDT score distributions for data and signal
- 'root -l -b BDT_optimal_cut.C' returns optimal cuts for categorisation running on BDT score distributions previously produced. Useful plots are produces in 'TMVA_inputpath'+category folders
- 'root -l evaluation_fillmass.C' fill mass distribution separately for the 6 categories. Output file: 'MVA_2018/datacardT3Mu_'+TMVA_inputpath+'_A.root'
