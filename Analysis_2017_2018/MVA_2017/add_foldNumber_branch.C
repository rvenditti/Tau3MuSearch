#include "TH1F.h"
#include <cmath>
#include <string>
#include "T3M_common.h"

void add_eventID_branch() {
//input arguments: numer of eventID for assignement in TMVA
//this will run on signal and MC files defined in T3M_common.h
    TString cat[] = {"A", "B", "C"};
    TString inputpath_datarun[] = {
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1436/AnalysedTree_data_2018A_tau3mu_12march.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1035/AnalysedTree_data_2018B_tau3mu_12march.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1036/AnalysedTree_data_2018C_tau3mu_12march.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1037/AnalysedTree_data_2018D_tau3mu_12march.root",
           };

    //Has to be the merging of the previous 4
    //    TString inputpath_data = "AnalysedTree_data_2018_merged_12march.root";
    //
    //        TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1038/AnalysedTree_MC_2018Ds_tau3mu_12march.root";
    //            TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1038/AnalysedTree_MC_2018B0_tau3mu_12march.root";
    //                TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1039/AnalysedTree_MC_2018Bp_tau3mu_12march.root";

    for(int j=0;j<4;j++){
        TFile *f = new TFile(inputpath_datarun[j],"update"); 
        for(int z=0;z<3;z++){
            TTree *T = (TTree*)f->Get("FinalTree"+cat[z]+"_sgn"); 
            int eventID; 
            TBranch *eventID = T->Branch("eventID",&eventID,"eventID/I"); 
            T->SetBranchAddress("px",&px); 
            T->SetBranchAddress("py",&py); 
            Long64_t nentries = T->GetEntries(); 
            for (Long64_t i=0;i<nentries;i++) { 
                T->GetEntry(i); 
                eventID = generate_random(nFolds); 
                bpt->Fill(); 
            } 
            T->Print(); 
            T->Write();
        } 
        delete f;
    }
 }



