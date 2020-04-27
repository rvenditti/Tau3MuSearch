#include <iostream>

//TMVA Training options
    TString TMVA_outputpath = "dataset_2017_22april_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = false;
    TString method = "BDT";
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath = "dataset_2017_22april_";  //name to load TMVA results for evaluation

//data rootfiles

    //22 april
    TString inputpath_datarun[] = {
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1255/AnalysedTree_data_2017B_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1255/AnalysedTree_data_2017C_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1256/AnalysedTree_data_2017D_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1258/AnalysedTree_data_2017E_tau3mu_22april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1258/AnalysedTree_data_2017F_tau3mu_22april.root",
           };

    TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1454/AnalysedTree_MC_2017Ds_tau3mu_22april.root";
    TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1455/AnalysedTree_MC_2017B0_tau3mu_22april.root";
    TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200422_1455/AnalysedTree_MC_2017Bp_tau3mu_22april.root";


//Coefficients for signal normalisation
    Double_t Dplus_correction = 1.05; // to be applied to D signal  
    Double_t Bs_correction = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction = 0.775; //2017 Ds correction factor updated on 28 November

    Double_t wNormDs = 1.242E-03; // Ds 3.67M initial events
    Double_t wNormB0 = 4.160E-04; // B0 3.00M initial events
    Double_t wNormBp = 6.203E-04; // Bp 2.01M initial events
    Double_t sig_norm = 0;


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



//TMVA settings
// Variables
    std::vector<TString> var_BDT = {
                                      "Pmu3>45?45:Pmu3",
                                      "cLP>30?30:cLP",
                                      "tKink>80?80:tKink",
                                      "segmComp<0.2?0.2:segmComp",
                                      "fv_nC>25?25:fv_nC",
                                      "fv_dphi3D>0.15?0.15:fv_dphi3D",
                                      "fv_d3Dsig>100?100:fv_d3Dsig",
                                      "d0sig>15?15:d0sig",
                                      "mindca_iso>0.5?0.5:mindca_iso",
                                      "trkRel>10?10:trkRel",
                                      "nMatchesMu3"
                                      };
    std::vector<TString> var_names = {
                                      "Pmu3",
                                      "cLP",
                                      "tKink",
                                      "segmComp",
                                      "fv_nC",
                                      "fv_dphi3D",
                                      "fv_d3Dsig",
                                      "d0sig",
                                      "mindca_iso",
                                      "trkRel",
                                      "nMatchesMu3"
                                      };


    std::vector<TString> var_spec_names = {
                                      "tripletMass",
                                      "puFactor",
                                     // "evt",
                                     // "evt % 8192",
                                      };

