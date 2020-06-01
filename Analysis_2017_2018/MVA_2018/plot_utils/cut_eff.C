#include "TH1F.h"
#include <cmath>
#include <iomanip>
#include <sstream>


void cut_eff() 
{
    //open root files where to store all plots!!
    //TFile *fout = new TFile("dsphipi_all_runs_output.root", "RECREATE");
    //fout->cd();
    //TString run_name[]={"2018A","2018B","2018C","2018D", "2018_MC_Ds","2018_MC_B0", "2018_MC_Bp"};
    //TString run_name[]={"2016B",
    //                    "2016C",
    //                    "2016D",
    //                    "2016E",
    //                    "2016F",
    //                    "2016G",
    //                    "2016H",
    //                    "2016_MC_Ds",
    //                    "2016_MC_B0", 
    //                    "2016_MC_Bp"};

    TString run_name[]={"2018A",
                        "2018B",
                        "2018C",
                        "2018D",
                        "2018_MC_Ds",
                        "2018_MC_B0", 
                        "2018_MC_Bp"};

    TString file_name_data[]={
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1158/AnalysedTree_data_2016B_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1218/AnalysedTree_data_2016C_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1221/AnalysedTree_data_2016D_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1223/AnalysedTree_data_2016E_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1224/AnalysedTree_data_2016F_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1233/AnalysedTree_data_2016G_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1236/AnalysedTree_data_2016H_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1240/AnalysedTree_MC_2016Ds_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1243/AnalysedTree_MC_2016B0_tau3mu_noChi_isPF_9feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200209_1246/AnalysedTree_MC_2016Bp_tau3mu_noChi_isPF_9feb.root"

//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200206_1108/AnalysedTree_data_2016B_tau3mu_noChi_6feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200206_1110/AnalysedTree_data_2016C_tau3mu_noChi_6feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200206_1111/AnalysedTree_data_2016D_tau3mu_noChi_6feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200206_1117/AnalysedTree_data_2016E_tau3mu_noChi_6feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200206_1119/AnalysedTree_data_2016F_tau3mu_noChi_6feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200206_1123/AnalysedTree_data_2016G_tau3mu_noChi_6feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200206_1127/AnalysedTree_data_2016H_tau3mu_noChi_6feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200207_2314/AnalysedTree_MC_2016Ds_tau3mu_noChi_7feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200207_2315/AnalysedTree_MC_2016B0_tau3mu_noChi_7feb.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200207_2316/AnalysedTree_MC_2016Bp_tau3mu_noChi_7feb.root"

     //      "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1154/AnalysedTree_data_2018A_tau3mu_17march.root",
     //      "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1155/AnalysedTree_data_2018B_tau3mu_17march.root",
     //      "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1156/AnalysedTree_data_2018C_tau3mu_17march.root",
     //      "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1201/AnalysedTree_data_2018D_tau3mu_17march.root",
    
     //      "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1211/AnalysedTree_MC_2018Ds_tau3mu_17march.root",
     //      "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1212/AnalysedTree_MC_2018B0_tau3mu_17march.root",
     //      "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1215/AnalysedTree_MC_2018Bp_tau3mu_17march.root",

    ////22 march
    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1957/AnalysedTree_data_2018A_tau3mu_22march.root",
    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_2003/AnalysedTree_data_2018B_tau3mu_22march.root",
    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_2005/AnalysedTree_data_2018C_tau3mu_22march.root",
    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1948/AnalysedTree_data_2018D_tau3mu_22march.root",
    //
    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200322_0935/AnalysedTree_MC_2018Ds_tau3mu_22march.root",
    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200322_1517/AnalysedTree_MC_2018B0_tau3mu_22march.root",
    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200322_1517/AnalysedTree_MC_2018Bp_tau3mu_22march.root"

    //12 may
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200512_1532/AnalysedTree_data_2018A_tau3mu_12may.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1713/AnalysedTree_data_2018B_tau3mu_6may.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1713/AnalysedTree_data_2018C_tau3mu_6may.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1713/AnalysedTree_data_2018D_tau3mu_6may.root",

           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1715/AnalysedTree_MC_2018Ds_tau3mu_6may.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1716/AnalysedTree_MC_2018B0_tau3mu_6may.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200506_1716/AnalysedTree_MC_2018Bp_tau3mu_6may.root"
    //21april
         //  "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_data_2018A_tau3mu_21april.root",
         //  "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018B_tau3mu_21april.root",
         //  "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018C_tau3mu_21april.root",
         //  "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018D_tau3mu_21april.root",

         //  "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018Ds_tau3mu_21april.root",
         //  "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018B0_tau3mu_21april.root",
         //  "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_MC_2018Bp_tau3mu_21april.root"

//    //12 march
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1436/AnalysedTree_data_2018A_tau3mu_12march.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1035/AnalysedTree_data_2018B_tau3mu_12march.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1036/AnalysedTree_data_2018C_tau3mu_12march.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1037/AnalysedTree_data_2018D_tau3mu_12march.root",
//
//      //     "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200414_1613/AnalysedTree_MC_2018Ds_tau3mu_14april_l1filter.root",
//      //     "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200414_1613/AnalysedTree_MC_2018B0_tau3mu_14april_l1filter.root",
//      //     "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200414_1613/AnalysedTree_MC_2018Bp_tau3mu_14april_l1filter.root"    
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200405_2009/AnalysedTree_MC_2018Ds_tau3mu_5april.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200405_2018/AnalysedTree_MC_2018B0_tau3mu_5april.root",
//    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1038/AnalysedTree_MC_2018Ds_tau3mu_12march.root",
//    //       "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1038/AnalysedTree_MC_2018B0_tau3mu_12march.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1039/AnalysedTree_MC_2018Bp_tau3mu_12march.root"
//
          // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200226_1228/AnalysedTree_data_2018A_tau3mu_26feb.root",
          // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200226_1002/AnalysedTree_data_2018B_tau3mu_26feb.root",
          // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200225_1749/AnalysedTree_data_2018C_tau3mu_26feb.root",
          // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200226_1225/AnalysedTree_data_2018D_tau3mu_26feb.root",

          // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200226_1120/AnalysedTree_MC_2018Ds_tau3mu_26feb.root",
          // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200226_1122/AnalysedTree_MC_2018B0_tau3mu_26feb.root",
          // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200226_1124/AnalysedTree_MC_2018Bp_tau3mu_26feb.root",

           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1825/AnalysedTree_data_2016B_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1826/AnalysedTree_data_2016C_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1828/AnalysedTree_data_2016D_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1829/AnalysedTree_data_2016E_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1849/AnalysedTree_data_2016F_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200215_1204/AnalysedTree_data_2016G_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1855/AnalysedTree_data_2016H_tau3mu_14feb_noOmega.root",

           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1856/AnalysedTree_MC_2016Ds_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1913/AnalysedTree_MC_2016B0_tau3mu_14feb_noOmega.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200214_1913/AnalysedTree_MC_2016Bp_tau3mu_14feb_noOmega.root"


//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200120_1032/AnalysedTree_data_2016B_tau3mu_20jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200120_1036/AnalysedTree_data_2016C_tau3mu_20jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200120_1039/AnalysedTree_data_2016D_tau3mu_20jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200120_1042/AnalysedTree_data_2016E_tau3mu_20jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200120_2023/AnalysedTree_data_2016F_tau3mu_20jan.root",
//                          //   "",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200121_0908/AnalysedTree_data_2016Hv2_tau3mu_20jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200120_2045/AnalysedTree_data_2016Hv3_tau3mu_20jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200128_1719/AnalysedTree_MC_2016Ds_tau3mu_28jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200128_1731/AnalysedTree_MC_2016B0_tau3mu_28jan.root",
//                             "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200128_1732/AnalysedTree_MC_2016Bp_tau3mu_28jan.root",
                             // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200124_1032/AnalysedTree_data_2018A_tau3mu_noChi_24jan.root",
                             // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200124_1036/AnalysedTree_data_2018B_tau3mu_noChi_24jan.root",
                             // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200124_1037/AnalysedTree_data_2018C_tau3mu_noChi_24jan.root",
                             // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200124_1038/AnalysedTree_data_2018D_tau3mu_noChi_24jan.root",
                             // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200124_0940/AnalysedTree_MC_2018Ds_tau3mu_noChi_24jan.root",
                             // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200124_0941/AnalysedTree_MC_2018B0_tau3mu_noChi_24jan.root",
                             // "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200124_0941/AnalysedTree_MC_2018Bp_tau3mu_noChi_24jan.root",
                             // "../ControlPlots/Control_v3/AnalysedTree_data_2017B_control_19dic.root",
                             // "../ControlPlots/Control_v3/AnalysedTree_data_2017C_control_21nov.root",
                             // "../ControlPlots/Control_v3/AnalysedTree_data_2017D_control_21nov.root",
                             // "../ControlPlots/Control_v3/AnalysedTree_data_2017E_control_21nov.root",
                             // "../ControlPlots/Control_v3/AnalysedTree_data_2017F_control_19dic.root"
                              };
    TString file_name_MC = "../ControlPlots/Control_v3/AnalysedTree_MC_DsPhiPi_control_21nov.root";

    for(int i=0; i<7; i++){

     //   TFile *f = new TFile(file_name_MC,"READ");
        TFile *f = new TFile(file_name_data[i],"READ");
        TH1F *h_eff_evt_data;
        h_eff_evt_data = (TH1F*)f->Get("CutEff_NEvents");
 //       TH1F *h_eff_trip_data;
 //       h_eff_trip_data = (TH1F*)f->Get("CutEff_Ntriplets");
 //       TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
 //       h_eff_evt_data->Draw();
 //       h_eff_trip_data->Draw("same");
        cout<<run_name[i]<<endl;
        cout<<file_name_data[i]<<endl;
        for(int j = 0; j<h_eff_evt_data->GetNbinsX(); j++) {
            cout<<(Double_t)h_eff_evt_data->GetBinContent(j)<<endl;
        }
    }    
    //fout->Write();
    //fout->Close();
}
