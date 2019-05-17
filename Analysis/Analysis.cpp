#include "ntupleClass_MC.C"
#include "ntupleClass_data.C"
#include "ntupleClass_Control.C"
#include <stdio.h>
#include <string.h>

// ###### INSTRUCTIONS (input options)
// (after having compiled the code [.L Analysis.cpp])
//
// * MC
// *** Ds->Tau->3Mu : Analysis("MC", "Ds")
// *** B0->Tau->3Mu : Analysis("MC", "B0")
// *** Bp->Tau->3Mu : Analysis("MC", "Bp")
// *** Ds->Phi->Pi  : Analysis("MC", "Ds_Phi")
// *** Minimum bias : Analysis("MC", "MiniBias")
// * Data
// *** 2017 B : Analysis("data", "2017B")
// *** 2017 F : Analysis("data", "2017F")


void Analysis(TString type, TString datasetName){
    TTree *tree;
    
    // Check input arguments
    if(strcmp(type, "MC") != 0 && strcmp(type, "data") != 0){
        cout << "The first argument is wrong! Please choose between 'MC' and 'data'" << endl;
        return;
    }
    if(strcmp(type, "MC") == 0 && strcmp(datasetName, "Ds") != 0 && strcmp(datasetName, "B0") != 0 && strcmp(datasetName, "Bp") != 0 && strcmp(datasetName, "Ds_Phi") != 0 && strcmp(datasetName, "MiniBias") != 0){
        cout << "The second argument is wrong! Please choose between 'Ds', 'B0', 'Bp', 'Ds_Phi' and 'MiniBias' " << endl;
        return;
    }
    if(strcmp(type, "data") == 0 && strcmp(datasetName, "2017B") != 0 && strcmp(datasetName, "2017F") != 0 ){
        cout << "The second argument is wrong! Please choose between '2017B' and '2017F'" << endl;
        return;
    }
    
    // ################ MC
    if (strcmp(type, "MC") == 0 && (strcmp(datasetName, "Ds") == 0 || strcmp(datasetName, "B0") == 0 || strcmp(datasetName, "Bp") == 0 || strcmp(datasetName, "MiniBias") == 0)){
        cout << "This is a MC" << endl;
        // Ds -> Tau -> 3Mu
        if (strcmp(datasetName, "Ds") == 0){
            cout << "MC Dataset : Ds -> Tau -> 3Mu" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree/Tree_DsTau3Mu.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("Tree/Tree_DsTau3Mu.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("Tree/Tree_DsTau3Mu.root:/Tree3Mu");
            dir->GetObject("ntuple",tree);
        }
        // B0 -> Tau -> 3Mu
        if (strcmp(datasetName, "B0") == 0){
            cout << "MC Dataset : B0 -> Tau -> 3Mu" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree/Tree_B0Tau3Mu.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("Tree/Tree_B0Tau3Mu.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("Tree/Tree_B0Tau3Mu.root:/Tree3Mu");
            dir->GetObject("ntuple",tree);
        }
        // Bp -> Tau -> 3Mu
        if (strcmp(datasetName, "Bp") == 0){
            cout << "MC Dataset : Bp -> Tau -> 3Mu" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree/Tree_BpTau3Mu.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("Tree/Tree_BpTau3Mu.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("Tree/Tree_BpTau3Mu.root:/Tree3Mu");
            dir->GetObject("ntuple",tree);
        }
        // Minimum bias
        if(strcmp(datasetName, "MiniBias") == 0){
            cout << "Minimum Bias Dataset " << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree/TreeMinBias.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("Tree/TreeMinBias.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("Tree/TreeMinBias.root:/TreeMakerBkg");
            dir->GetObject("ntuple",tree);
        }
//        ntupleClass_MC class_MC(tree);
        ntupleClass_MC class_MC(tree);
        class_MC.LoopMC_New(datasetName);
    }
    // Ds -> Phi -> Pi
    if(strcmp(type, "MC") == 0 && strcmp(datasetName, "Ds_Phi") == 0){
        cout << "MC Dataset : Ds -> Phi -> Pi" << endl << endl;
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree/TreeDsPhiPi.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("Tree/TreeDsPhiPi.root");
        }
        TDirectory * dir = (TDirectory*)f->Get("Tree/TreeDsPhiPi.root:/Tree3Mu");
        dir->GetObject("ntuple",tree);
        
        ntupleClass_Control class_Control(tree);
        class_Control.LoopControl();
    }
    
    // ##################### Data
    if (strcmp(type, "data") == 0){
        cout << "This is data" << endl;
        // 2017B
        if (strcmp(datasetName, "2017B") == 0){
            cout << "Data 2017 B" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree/Tree_DoubleMuonLowMass_Run2017B.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("Tree/Tree_DoubleMuonLowMass_Run2017B.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("Tree/Tree_DoubleMuonLowMass_Run2017B.root:/TreeMakerBkg");
            dir->GetObject("ntuple",tree);
        }
        // 2017F
        if (strcmp(datasetName, "2017F") == 0){
            cout << "Data 2017 F" << endl << endl;
            TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Tree/TreeData_Run2017F.root");
            if (!f || !f->IsOpen()) {
                f = new TFile("Tree/TreeData_Run2017F.root");
            }
            TDirectory * dir = (TDirectory*)f->Get("Tree/TreeData_Run2017F.root:/TreeMakerBkg");
            dir->GetObject("ntuple",tree);
        }
        
        ntupleClass_MC class_Data(tree);
        class_Data.LoopData_New(datasetName);
    }
    
}
