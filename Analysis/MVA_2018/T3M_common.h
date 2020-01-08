//TMVA Training options
    TString TMVA_outputpath = "dataset_02dic"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = true;
    TString method = "BDT";

//TMVA Evaluating options
    TString TMVA_inputpath = "dataset_02dic";  //name to load TMVA results for evaluation
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDT.weights.xml"; //name given training BDT with crossvalidation

//data rootfiles
    TString inputpath_datarun[] = {
           "AnalysedTree_data_2018A_tau3mu_02dic.root",
           "AnalysedTree_data_2018B_tau3mu_02dic.root",
           "AnalysedTree_data_2018C_tau3mu_02dic.root",
           "AnalysedTree_data_2018D_tau3mu_02dic.root"
           };
    //Has to be the merging of the previous 5
    TString inputpath_data = "AnalysedTree_data_2018_merged_02dic.root";
    
    TString inputpath_Ds = "AnalysedTree_MC_2018Ds_tau3mu_04dic.root";
    TString inputpath_B0 = "";
    TString inputpath_Bp = "";

//Coefficients for signal normalisation
    Double_t Ds_correction = 1;
    Double_t Dplus_correction = 1; 
    Double_t wNormDs = 0.0142;
    Double_t wNormB0 = 0;
    Double_t wNormBp = 0;
    Double_t sig_norm = 0.0142*Dplus_correction; //average normalization factor for the three signal samples

