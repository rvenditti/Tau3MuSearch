#include <iostream>
#include <fstream>

using namespace std;

//TMVA Training options
    TString TMVA_outputpath = "dataset_20apr"; //name to give to TMVA output files
    TString BDT_inVar = "BDT_inputVar.txt"; // input variables file name
    //change it to perform 5-fold Cross Validation
    bool doCV = false;
    TString method = "BDT";
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way

   // if(doCV)
//   TString method = "BDTG";
//   TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation

//TMVA Evaluating options
    TString TMVA_inputpath = "dataset_20apr";  //name to load TMVA results for evaluation

//data rootfiles (merging of the previous 4)

    TString inputpath_datarun[] = { "/Users/caterinaaruta/Desktop/LastTreeParking/Cate4/AnalysedTree_data_ParkingBkg_All_tau3mu_Cate4.root" };

    //Has to be the merging of the previous 4
    TString inputpath_data = "/Users/caterinaaruta/Desktop/LastTreeParking/Cate4/AnalysedTree_data_ParkingBkg_All_tau3mu_Cate4.root";

// MC rootfiles    
    TString inputpath_Ds = "/Users/caterinaaruta/Desktop/LastTreeParking/Cate4/AnalysedTree_data_ParkingDs_tau3mu_Cate4_new.root";
    TString inputpath_B0 = "/Users/caterinaaruta/Desktop/LastTreeParking/Cate4/AnalysedTree_data_ParkingB0_tau3mu_Cate4_new.root";
    TString inputpath_Bp = "/Users/caterinaaruta/Desktop/LastTreeParking/Cate4/AnalysedTree_data_ParkingBp_tau3mu_Cate4.root";

    ////Has to be the merging of the previous 3 //just for plotting!!!
    TString inputpath_MC = "/Users/caterinaaruta/Desktop/LastTreeParking/Cate4/AnalysedTree_data_ParkingSgn_All_tau3mu_Cate4_new.root";

    Double_t Dplus_correction = 1.05; // to be applied to D signal  
    Double_t Bs_correction = 1.12; // to be applied to B0 and Bp signal  
    Double_t f_correction = 1.; // to be applied to B0 and Bp signal  
    Double_t Ds_correction = 1; //

    Double_t wNormDs = 6.49E-04; //MC D3mu Jian
    Double_t wNormB0 = 2.41E-04;
    Double_t wNormBp = 7.20E-04;
    Double_t sig_norm = 0;
//(wNormDs+wNormB0+wNormBp) * Dplus_correction / 3.; //average normalization factor for the three signal samples

    vector<TString> VarName, VarSel;
    vector<char> VarType;

    void ImportVar(vector<TString> &VarName, vector<TString> &VarSel, vector<char> &VarType){
        TString varArg[2]; char vartype;
        ifstream inputVar;
        inputVar.open(BDT_inVar);
        while (!inputVar.fail() && !inputVar.eof()){
            inputVar >> varArg[0] >> varArg[1] >> vartype;
            VarName.push_back(varArg[0]);
            VarSel.push_back(varArg[1]);
            VarType.push_back(vartype);
        }
        inputVar.close();
    }
