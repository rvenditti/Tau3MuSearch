#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>


void cut_eff() 
{
    //Print to terminal selection efficiencies
    TString run_name[]={"2017B",
                        "2017C",
                        "2017D",
                        "2017E",
                        "2017F",
                        "2017_MC_Ds",
                        "2017_MC_B0", 
                        "2017_MC_Bp"};

    TString file_name_data[]={

    //22april
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1255/AnalysedTree_data_2017B_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1255/AnalysedTree_data_2017C_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1256/AnalysedTree_data_2017D_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1258/AnalysedTree_data_2017E_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1258/AnalysedTree_data_2017F_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1454/AnalysedTree_MC_2017Ds_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1455/AnalysedTree_MC_2017B0_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1455/AnalysedTree_MC_2017Bp_tau3mu_22april.root"
                              };
    for(int i=0; i<8; i++){

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
