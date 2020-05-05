#include <iostream>
#include <fstream>

using namespace std;

//TMVA Training options

    TString chi2cut = "";
    TString TMVA_outputpath = "dataset_2018_4may_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = false;
    TString method = "BDT";
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
 
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath = "dataset_2018_4may_";  //name to load TMVA results for evaluation

//data rootfiles

    //21 april: muonID in miniTree
    TString inputpath_datarun[] = {
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200430_1121/AnalysedTree_data_2018B_tau3mu_30april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200430_1121/AnalysedTree_data_2018C_tau3mu_30april.root",
           "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200430_1121/AnalysedTree_data_2018D_tau3mu_30april.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_data_2018A_tau3mu_21april.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018B_tau3mu_21april.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018C_tau3mu_21april.root",
           //"/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0826/AnalysedTree_data_2018D_tau3mu_21april.root",
           };

    //TString inputpath_data = "AnalysedTree_data_2018_merged_21april.root";
    TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200430_1123/AnalysedTree_MC_2018Ds_tau3mu_30april.root";
    TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200430_1123/AnalysedTree_MC_2018B0_tau3mu_30april.root";
    TString inputpath_Bp = "";
    //TString inputpath_Ds = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018Ds_tau3mu_21april.root";
    //TString inputpath_B0 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0824/AnalysedTree_MC_2018B0_tau3mu_21april.root";
    //TString inputpath_Bp = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20200421_0825/AnalysedTree_MC_2018Bp_tau3mu_21april.root";


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
    TString BDTinVar_A = "BDTinputVar_"+TMVA_outputpath+"A.txt"; 
    TString BDTinVar_B = "BDTinputVar_"+TMVA_outputpath+"B.txt"; 
    TString BDTinVar_C = "BDTinputVar_"+TMVA_outputpath+"C.txt"; 
    TString BDTspecVar = "BDTspecVar_"+TMVA_outputpath+".txt";

//Utility to read variables from text file
void readVarName(std::vector<TString> &var_name, std::vector<TString> &var_def,  TString filename ){
     TString var[2];
     ifstream inputVar;
     inputVar.open(filename);
     while (!inputVar.fail() && !inputVar.eof()){
         inputVar >> var[0] >> var[1];
         var_name.push_back(var[0]);
         var_def.push_back(var[1]);
     }
     inputVar.close();
}
