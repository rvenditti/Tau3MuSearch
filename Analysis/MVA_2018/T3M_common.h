#include <iostream>

//TMVA Training options
    TString TMVA_outputpath = "dataset_2018_22april_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = false;
    TString method = "BDT";
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath = "dataset_2018_22april_";  //name to load TMVA results for evaluation

//data rootfiles

    //21 april: muonID in miniTree
    TString inputpath_datarun[] = {
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_data_2018A_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018B_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018C_tau3mu_21april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018D_tau3mu_21april.root",
           };

    TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018Ds_tau3mu_21april.root";
    TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018B0_tau3mu_21april.root";
    TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_MC_2018Bp_tau3mu_21april.root";


//Coefficients for signal normalisation
    Double_t Dplus_correction = 1.05; // to be applied to D signal  
    Double_t Bs_correction = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction = 0.93; //2018 updated on 22 march

    Double_t wNormDs = 1.323E-03; //new 4.6E+06 MC DsTau
  //        Double_t wNormDs = 7.15E-03; //old 2.13E+06 MC DsTau
    Double_t wNormB0 = 4.784E-04;  //new 3.5E+06 MC B0Tau
  //        Double_t wNormB0 = 3.03E-03;  //old 0.655E+06 MC B0Tau
    Double_t wNormBp = 1.437E-03;  //usual 1.33E+06 MC BpTau
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

