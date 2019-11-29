**Train and Test MVA methods: example of call**
cd <work_directory>/CMSSW_9_4_4/src
cmsenv
cd <work_directory>/CMSSW_9_4_4/src/Analysis/MVAeval_datacards
root -l MVA_Code.cpp\(\”cat\”\)
where cat = “A”, “B”, “C”, “B+C” (different combinations are easy to implement)

Input files for MVA_Code.cpp are the “AnalysedTree_*.root” produced analysing the ntuples. The training is done using the small ntuples stored in the root files containing variables for BDT, pileup-reweighting factor and triplet invariant mass.

**Evaluation of BDT for event classification and final invariant mass distributions:**

Macro RunT3Mu.C:
- read the BDT decision distribution on test samples and set the thresholds based on combined significance using function defined in Get_BDT_cut.C

- loops on analysed ntuples, read the weights for the given category (dataset_"+cat+"/weights/TMVA_new_BDT.weights.xml), evaluate the BDT decision and categorise the event depending on the cut set before
-  fill the triplet invariant mass in output histogram "datacardT3Mu_"+cat+".root”, which is normalised in case of MC signal.

Code RunT3Mu_launcher.cpp:
- Book the category name, prepare the input ntuples as ROOT TChain
- launch “Loop” method on each group of events (data, signal Ds, B0, Bp)
- Create inclusive distribution of the 3 signal samples and close output files

Usage example:
compile RunT3Mu.C

g++ -I $ROOTSYS/include -I $ROOTSYS/include RunT3Mu_launcher.cpp `root-config --glibs` `root-config --libs` `root-config --cflags` -lTMVA -L $ROOTSYS/lib -o RunT3Mu_exe

execute ./RunT3Mu_exe “cat”

Macro BDT_optimal_cut.C:
- Set cut values on BDT scuse (same as Get_BDT_cut.C)
- Plot BDT score on test sample normalized to 1 with the cut values
- Plot combined significance vs cut values used for optimization
