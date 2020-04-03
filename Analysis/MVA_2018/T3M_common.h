#include <iostream>

//TMVA Training options
    TString TMVA_outputpath = "dataset_12march_2018_CV_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = true;
   // TString method = "BDT";
   // TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   TString method = "BDTG";
   TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath = "dataset_12march_2018_CV_";  //name to load TMVA results for evaluation

//data rootfiles
   // //17 march
   // TString inputpath_datarun[] = {
   //        "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1154/AnalysedTree_data_2018A_tau3mu_17march.root",
   //        "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1155/AnalysedTree_data_2018B_tau3mu_17march.root",
   //        "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1156/AnalysedTree_data_2018C_tau3mu_17march.root",
   //        "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1201/AnalysedTree_data_2018D_tau3mu_17march.root",
   //        };

   // //Has to be the merging of the previous 4
   // TString inputpath_data = "AnalysedTree_data_2018_merged_17march.root";
   // 
   // TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1211/AnalysedTree_MC_2018Ds_tau3mu_17march.root";
   // TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1212/AnalysedTree_MC_2018B0_tau3mu_17march.root";
   // TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200316_1215/AnalysedTree_MC_2018Bp_tau3mu_17march.root";

   // //Has to be the merging of the previous 3 //just for plotting!!!
   // TString inputpath_MC = "AnalysedTree_MC_2018_merged_tau3mu_17march.root";

    //12 march
    TString inputpath_datarun[] = {
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1436/AnalysedTree_data_2018A_tau3mu_12march.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1035/AnalysedTree_data_2018B_tau3mu_12march.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1036/AnalysedTree_data_2018C_tau3mu_12march.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1037/AnalysedTree_data_2018D_tau3mu_12march.root",
           };

    //Has to be the merging of the previous 4
    TString inputpath_data = "AnalysedTree_data_2018_merged_12march.root";
    
    TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1038/AnalysedTree_MC_2018Ds_tau3mu_12march.root";
    TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1038/AnalysedTree_MC_2018B0_tau3mu_12march.root";
    TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200312_1039/AnalysedTree_MC_2018Bp_tau3mu_12march.root";

    ////Has to be the merging of the previous 3 //just for plotting!!!
    TString inputpath_MC = "AnalysedTree_MC_2018_merged_tau3mu_12march.root";

//    //22 march
//    TString inputpath_datarun[] = {
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1957/AnalysedTree_data_2018A_tau3mu_22march.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_2003/AnalysedTree_data_2018B_tau3mu_22march.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_2005/AnalysedTree_data_2018C_tau3mu_22march.root",
//           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200321_1948/AnalysedTree_data_2018D_tau3mu_22march.root",
//           };
//
//    //Has to be the merging of the previous 4
//    TString inputpath_data = "AnalysedTree_data_2018_merged_tau3mu_22march.root";
//    
//    TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200322_0935/AnalysedTree_MC_2018Ds_tau3mu_22march.root";
//    TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200322_1517/AnalysedTree_MC_2018B0_tau3mu_22march.root";
//    TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200322_1517/AnalysedTree_MC_2018Bp_tau3mu_22march.root";
//
//    ////Has to be the merging of the previous 3 //just for plotting!!!
//    TString inputpath_MC = "AnalysedTree_MC_2018_merged_tau3mu_22march.root";

//Coefficients for signal normalisation
    Double_t Dplus_correction = 1.05; // to be applied to D signal  
    Double_t Bs_correction = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction = 0.93; //2018 updated on 22 march

    Double_t wNormDs = 7.15293E-03; //MC D3mu Jian
    Double_t wNormB0 = 3.0252E-03;
    Double_t wNormBp = 1.4374E-03;
    Double_t sig_norm = 0;
//(wNormDs+wNormB0+wNormBp) * Dplus_correction / 3.; //average normalization factor for the three signal samples


//TMVA settings
//
// Variables declaration
    TString var_Pmu3 =       "Pmu3>45?45:Pmu3";
    TString var_cLP =        "cLP>30?30:cLP";
    TString var_tKink =      "tKink>80?80:tKink";
    TString var_segmComp =   "segmComp<0.2?0.2:segmComp";
    TString var_fv_nC =      "fv_nC>25?25:fv_nC";
    TString var_fv_dphi3D =  "fv_dphi3D>0.15?0.15:fv_dphi3D"; 
    TString var_fv_d3Dsig =  "fv_d3Dsig>100?100:fv_d3Dsig";
    TString var_d0sig =      "d0sig>15?15:d0sig";
    TString var_mindca_iso = "mindca_iso>0.5?0.5:mindca_iso";
    TString var_trkRel =     "trkRel>10?10:trkRel";
