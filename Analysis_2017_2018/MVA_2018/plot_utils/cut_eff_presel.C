#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>


void cut_eff_presel() 
{
    //Print to terminal selection efficiencies
    TString file_name[]={
                         //"/lustrehome/fsimone/Run2017B_AOD_v8_merged.root",
                         "/lustrehome/fsimone/2017C_AOD.root",
                         //"/lustrehome/fsimone/Run2017D_AOD_v8_merged.root",
                         //"/lustrehome/fsimone/Run2017E_AOD_v8_merged.root",
                         //"/lustrehome/fsimone/Run2017F_AOD_v8_merged.root",
                         "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200424_1854/MiniAOD_check.root",
                         //"/lustrehome/fsimone/2017MCDs_v8_merged.root",
                         "/lustrehome/fsimone/2017MCB0_v8_merged.root",
                         "/lustrehome/fsimone/2017MCBp_v8_merged.root"

                        };

    TString plots[] = {
                       "InitialPlots/hEvtCount",
                       "PlotsAfterTrigger/hEvtCount",
                       "PlotsAfterLooseMuon/hEvtCount",
                       "PlotsAfter3Muons/hEvtCount",
                       "PlotsAfterTauCand/hEvtCount",
                       //"Tree3Mu/hEventsAfterGoodCand"
                       "TreeMakerBkg/hEventsAfterGoodCand"
                       };
    //number of input files
    size_t n_file = sizeof(file_name)/sizeof(file_name[0]);
    //number of histograms
    size_t n_histo = sizeof(plots)/sizeof(plots[0]);
    for(auto i = 0; i<n_file; i++){
        cout<<file_name[i]<<endl;
        for(auto j = 0; j<n_histo; j++){
            TFile *f = new TFile(file_name[i],"READ");
            TH1F *h_eff;
            h_eff = (TH1F*)f->Get(plots[j]);
            cout<<(Double_t)h_eff->GetEntries()<<endl;
        }
    }
}
