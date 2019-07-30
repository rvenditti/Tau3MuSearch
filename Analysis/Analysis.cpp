#include "ntupleClass_MC.C"
#include "ntupleClass_data.C"
#include "ntupleClass_Control.C"
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
    
    // Check input arguments
    if(strcmp(type, "MC") != 0 && strcmp(type, "MC_sgn") != 0 && strcmp(type, "MC_control") != 0 && strcmp(type, "data") != 0 && strcmp(type, "data_bkg") != 0 && strcmp(type, "data_control") != 0){
        cout << "The first argument is wrong! Please choose between 'MC', 'MC_sgn', 'MC_control', 'data', 'data_bkg'and 'data_control'" << endl;
        return -1;
    }
    if((strcmp(type, "MC") == 0 || strcmp(type, "MC_sgn") == 0 ) && strcmp(datasetName, "Ds") != 0 && strcmp(datasetName, "B0") != 0 && strcmp(datasetName, "Bp") != 0 && strcmp(datasetName, "MiniBias") != 0){
        cout << "The second argument is wrong! Please choose between 'Ds', 'B0', 'Bp' and 'MiniBias' " << endl;
        return -1;
    }
    /*
    if((strcmp(type, "data") == 0 || strcmp(type, "data_bkg") == 0 || strcmp(type, "data_control") == 0) && strcmp(datasetName, "2017B") != 0 && strcmp(datasetName, "2017B_0a") != 0 && strcmp(datasetName, "2017B_0b") != 0 && strcmp(datasetName, "2017B_0c") != 0 && strcmp(datasetName, "2017B_1") != 0 && strcmp(datasetName, "2017C") != 0 && strcmp(datasetName, "2017D") != 0 && strcmp(datasetName, "2017D_0") != 0 && strcmp(datasetName, "2017D_1") != 0  && strcmp(datasetName, "2017F") != 0 && strcmp(datasetName, "2017F_0a") != 0 && strcmp(datasetName, "2017F_0b") != 0 && strcmp(datasetName, "2017F_0c") != 0 && strcmp(datasetName, "2017F_0d") != 0 && strcmp(datasetName, "2017F_0e") != 0 && strcmp(datasetName, "2017F_1a") != 0 && strcmp(datasetName, "2017F_1b") != 0 && strcmp(datasetName, "2017F_1c") != 0 && strcmp(datasetName, "2017F_1d") != 0 && strcmp(datasetName, "2017F_1e") != 0 && strcmp(datasetName, "2017F_2") != 0 ){
        cout << "The second argument is wrong! Please choose between '2017B', '2017B_0a', '2017B_0b', '2017B_0c', '2017B_1', '2017C', '2017D', '2017F_0a', '2017F_0b', '2017F_0c', '2017F_0d', '2017F_0e', '2017F_1a', '2017F_1b', '2017F_1c', '2017F_1d', '2017F_1e', '2017F_2' and '2017F'" << endl;
        return -1;
    }
    */
    // ################ MC
    if ((strcmp(type, "MC") == 0 || strcmp(type, "MC_sgn") == 0) && (strcmp(datasetName, "Ds") == 0 || strcmp(datasetName, "B0") == 0 || strcmp(datasetName, "Bp") == 0 || strcmp(datasetName, "MiniBias") == 0)){
        cout << "This is a MC" << endl;
        // Ds -> Tau -> 3Mu
        if (strcmp(datasetName, "Ds") == 0){
            cout << "MC Dataset : Ds -> Tau -> 3Mu" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Tree/Tree_DsTau3Mu.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("../Tree/Tree_DsTau3Mu.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("../Tree/Tree_DsTau3Mu.root:/Tree3Mu");
            dir->GetObject("ntuple",tree);
        }
        // B0 -> Tau -> 3Mu
        if (strcmp(datasetName, "B0") == 0){
            cout << "MC Dataset : B0 -> Tau -> 3Mu" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Tree/Tree_B0Tau3Mu.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("../Tree/Tree_B0Tau3Mu.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("../Tree/Tree_B0Tau3Mu.root:/Tree3Mu");
            dir->GetObject("ntuple",tree);
        }
        // Bp -> Tau -> 3Mu
        if (strcmp(datasetName, "Bp") == 0){
            cout << "MC Dataset : Bp -> Tau -> 3Mu" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Tree/Tree_BpTau3Mu.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("../Tree/Tree_BpTau3Mu.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("../Tree/Tree_BpTau3Mu.root:/Tree3Mu");
            dir->GetObject("ntuple",tree);
        }
        // Minimum bias
        if(strcmp(datasetName, "MiniBias") == 0){
            cout << "Minimum Bias Dataset " << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Tree/TreeMinBias.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("../Tree/TreeMinBias.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("../Tree/TreeMinBias.root:/TreeMakerBkg");
            dir->GetObject("ntuple",tree);
        }
        ntupleClass_MC class_MC(tree);
        class_MC.LoopMC_New(type, datasetName);
        
    }
    // Ds -> Phi -> Pi
    if(strcmp(type, "MC_control") == 0){
        cout << "MC Dataset : Ds -> Phi -> Pi" << endl << endl;
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../Tree/TreeDsPhiPi.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("../Tree/TreeDsPhiPi.root");
        }
        TDirectory * dir = (TDirectory*)f->Get("../Tree/TreeDsPhiPi.root:/Tree3Mu");
        dir->GetObject("ntuple",tree);
        
        ntupleClass_Control class_Control(tree);
        class_Control.LoopControl();
    }
    
    // ##################### Data
    if (strcmp(type, "data") == 0 || strcmp(type, "data_bkg") == 0){
        cout << "This is data" << endl;
        cout << "Data " << datasetName << endl << endl;
        TString treeName = "/lustre/cms/store/user/rosma/DoubleMuonLowMass/TreeDsTau3Mu/TreeData_Run";  treeName += datasetName; treeName += ".root";
        TString treeMakerName = "/lustre/cms/store/user/rosma/DoubleMuonLowMass/TreeDsTau3Mu/TreeData_Run";  treeMakerName += datasetName;
        treeMakerName += ".root:/TreeMakerBkg";
        
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(treeName);
        if (!f || !f->IsOpen()) {
            f = new TFile(treeName);
        }
        TDirectory * dir = (TDirectory*)f->Get(treeMakerName);
        dir->GetObject("ntuple",tree);

        ntupleClass_MC class_Data(tree);
        class_Data.LoopData_New(type, datasetName);
    }
    
    // ##################### Data Control
    
    if (strcmp(type, "data_control") == 0){
        cout << "Control channel analysis on data" << endl;
        cout << "Data " << datasetName << endl << endl;
        TString treeName = "/lustre/cms/store/user/rosma/DoubleMuonLowMass/TreeDsPhiPiData/Tree_DsPhiPi_Data_Run"; treeName += datasetName; treeName += "_new.root";
        TString treeMakerName = "/lustre/cms/store/user/rosma/DoubleMuonLowMass/TreeDsPhiPiData/Tree_DsPhiPi_Data_Run"; treeMakerName += datasetName; treeMakerName += "_new.root:/Tree3Mu";
        
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(treeName);
        if (!f || !f->IsOpen()) {
            f = new TFile(treeName);
        }
        TDirectory * dir = (TDirectory*)f->Get(treeMakerName);
        dir->GetObject("ntuple",tree);
    
        ntupleClass_Control class_Data(tree);
        class_Data.LoopControl_Data(type, datasetName);
    }
    
    return 0;
}



