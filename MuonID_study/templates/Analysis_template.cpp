#include "ntupleClass_muonid.C"
#include <TROOT.h>
#include <stdio.h>
#include <iostream>
#include <string.h>

// ###### INSTRUCTIONS (input options)
// (after having compiled the code [.L Analysis.cpp])
//
// * type = "MC" -> for the 'standard analysis' (triplet mass quite large)
// * type = "MC_sgn" -> for the analysis for MVA  (triplet mass narrower)
// *** Ds->Tau->3Mu : Analysis("MC", "Ds")
// *** B0->Tau->3Mu : Analysis("MC", "B0")
// *** Bp->Tau->3Mu : Analysis("MC", "Bp")
// *** Ds->Phi->Pi  : Analysis("MC", "Ds_Phi")
// *** Minimum bias : Analysis("MC", "MiniBias")
//
// * type = "data" -> for the 'standard analysis' (triplet mass quite large)
// * type = "data_bkg" -> for the analysis for MVA (triplet mass complementary w.r.t. tau mass)
// *** 2017 B : Analysis("data", "2017B")
// *** 2017 C : Analysis("data", "2017C")
// *** 2017 D : Analysis("data", "2017D")
// *** 2017 F : Analysis("data", "2017F")

using namespace std;

int main(int narg, char** arg){
    TTree *tree;
    char type[10];
    strcpy(type, arg[1]);
    cout << "type : " << type << endl;
    char datasetName[15];
    strcpy(datasetName, arg[2]);
    cout << "datasetName : " << datasetName << endl << endl;
    TString fileout = "";
    
    // Check input arguments
    if(strcmp(type, "MC") != 0 && strcmp(type, "MC_sgn") != 0 && strcmp(type, "MC_control") != 0 && strcmp(type, "data") != 0 && strcmp(type, "data_bkg") != 0 && strcmp(type, "data_control_sgn") != 0 && strcmp(type, "data_control_bkg") != 0 && strcmp(type, "data_control") != 0){
        cout << "The first argument is wrong! Please choose between 'MC', 'MC_sgn', 'MC_control', 'data','data_bkg', 'data_control', 'data_control_sgn' and 'data_control_bkg'" << endl;
        return -1;
    }
    if((strcmp(type, "MC") == 0 || strcmp(type, "MC_sgn") == 0 ) && 
        strcmp(datasetName, "2018BdToKK") != 0 && strcmp(datasetName, "2018BdToPiPi") != 0 &&
        strcmp(datasetName, "2018Ds") != 0 && strcmp(datasetName, "2018B0") != 0 && strcmp(datasetName, "2018Bp") != 0 && 
        strcmp(datasetName, "MiniBias") != 0){
        cout << "The second argument is wrong!" << endl;
        return -1;
    }

    // ################ MC
    if (strcmp(type, "MC") == 0){
        cout << "This is a MC" << endl;

        // Ds -> Tau -> 3Mu
        if (strcmp(datasetName, "2018Ds") == 0){
            cout << "MC Dataset : Ds -> Tau -> 3Mu" << endl << endl;
            TChain* chain = new TChain("TreeMakerBkg/ntuple");
        //AddFile_MC2018Ds_muonid
        //OutFile_MC2018Ds_muonid

        ntupleClass_muonid class_MC(chain, fileout);
        class_MC.Loop(type, datasetName);
        }
        // B0 -> Tau -> 3Mu
        if (strcmp(datasetName, "2018B0") == 0){
            cout << "MC Dataset : B0 -> Tau -> 3Mu" << endl << endl;
            TChain* chain = new TChain("TreeMakerBkg/ntuple");
        //AddFile_MC2018B0_muonid
        //OutFile_MC2018B0_muonid

        ntupleClass_muonid class_MC(chain, fileout);
        class_MC.Loop(type, datasetName);
        }
        // Bp -> Tau -> 3Mu
        if (strcmp(datasetName, "2018Bp") == 0){
            cout << "MC Dataset : Bp -> Tau -> 3Mu" << endl << endl;
            TChain* chain = new TChain("TreeMakerBkg/ntuple");
        //AddFile_MC2018Bp_muonid
        //OutFile_MC2018Bp_muonid

        ntupleClass_muonid class_MC(chain, fileout);
        class_MC.Loop(type, datasetName);
        }
        // BdKK or BdPiPi
        if (strcmp(datasetName, "2018BdToPiPi") == 0 || strcmp(datasetName, "2018BdToKK") == 0){
            cout << "MC Dataset : Bd to KK or Bd to PiPi for muon ID studies" << endl << endl;
            TChain* chain = new TChain("TreeMakerBkg/ntuple");
        //AddFile_MC2018Bd_muonid
        //OutFile_MC2018Bd_muonid

        ntupleClass_muonid class_MC(chain, fileout);
        class_MC.Loop(type, datasetName);
        }
    }
    
    // ##################### Data
    if (strcmp(type, "data") == 0 || strcmp(type, "data_bkg") == 0){
        cout << "This is data" << endl;
        cout << "Data " << datasetName << endl << endl;
        TChain* chain = new TChain("TreeMakerBkg/ntuple");
        //AddFile_data_muonid
        //OutFile_data_muonid

        ntupleClass_muonid class_Data(chain, fileout);
        class_Data.Loop(type, datasetName);
    }
    
    return 0;
}
