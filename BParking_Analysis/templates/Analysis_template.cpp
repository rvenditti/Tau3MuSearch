#include "ntupleClass_Tau3mu.C"
#include "ntupleClass_Control.C"
#include <TROOT.h>
#include <stdio.h>
#include <iostream>
#include <string.h>

using namespace std;

int main(int narg, char** arg){
    TTree *tree;
    char type[10];
    strcpy(type, arg[1]);
    cout << "type : " << type << endl;
    char datasetName[15];
    strcpy(datasetName, arg[2]);
    cout << "datasetName : " << datasetName << endl << endl;
    
    // Check input arguments
    if(strcmp(type, "data") != 0 && strcmp(type, "data_control") != 0){ cout << "The first argument is wrong! Please choose 'data'" << endl; return -1; }
    
    // ##################### Data
    if (strcmp(type, "data") == 0){
        cout << "This is data" << endl;
        cout << "Data " << datasetName << endl << endl;
//        TChain* chain = new TChain("TreeMakerBkg/ntuple");
        TChain* chain = new TChain("ntuple");
        //TChain* chain = new TChain("Tree3Mu/ntuple");
        //AddFile_data_tau3mu
            TString fileout = "AddOutput_data_tau3mu.root";

        ntupleClass_Tau3mu class_Data(chain, fileout);
        class_Data.LoopTau3mu(type, datasetName);
    }
    
    // ##################### Data Control
    
    if (strcmp(type, "data_control") == 0){
        cout << "Control channel analysis on data" << endl;
        cout << "Data " << datasetName << endl << endl;
        //TChain* chain = new TChain("Tree3Mu/ntuple");
        TChain* chain = new TChain("ntuple");
        //AddFile_data_control
            TString fileout = "AddOutput_data_control.root";

        ntupleClass_Control class_Data(chain, fileout);
        class_Data.LoopControl(type, datasetName);
    }
    
    return 0;
}
