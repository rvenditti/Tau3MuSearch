#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>


void cut_eff() 
{
    //Print to terminal selection efficiencies
    TString run_name[]={"2018A",
                        "2018B",
                        "2018C",
                        "2018D",
                        "2018_MC_Ds",
                        "2018_MC_B0", 
                        "2018_MC_Bp"};

    TString file_name_data[]={

    //21april
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_data_2018A_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018B_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018C_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018D_tau3mu_21april.root",

           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018Ds_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018B0_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_MC_2018Bp_tau3mu_21april.root"
                              };
    for(int i=0; i<7; i++){

        TFile *f = new TFile(file_name_data[i],"READ");
        TH1F *h_eff_evt_data;
        h_eff_evt_data = (TH1F*)f->Get("CutEff_NEvents");
        cout<<run_name[i]<<endl;
        cout<<file_name_data[i]<<endl;
        for(int j = 0; j<h_eff_evt_data->GetNbinsX(); j++) {
            cout<<(Double_t)h_eff_evt_data->GetBinContent(j)<<endl;
        }
    }    
}
