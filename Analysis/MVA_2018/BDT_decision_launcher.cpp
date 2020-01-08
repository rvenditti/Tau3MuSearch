#include "BDT_decision.C"
#include <TROOT.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include "BDT_optimal_cut.h"

using namespace std;

int main(int narg, char** cat){
    TTree *tree;
    int isMC;
    TString category = cat[1];
    std::cout << "category: " << category << std::endl;
    
    // Check input arguments
    if (!(category == "A" || category == "B" || category == "C")) {
        std::cout<<"Wrong category name, please chooose between: A B C"<<std::endl;
        return 0;
    }
    //prepare output file
    TString fout_path = TMVA_inputpath+category+"/BDTdecision_"+category+".root";
    TFile *fout = new TFile(fout_path, "recreate");
    fout->cd();
    fout->Close();

    //data
    isMC = 0;
    TFile *f0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_data);
    if (!f0 || !f0->IsOpen()) f0 = new TFile(inputpath_data);
    std::cout<<"Opened input file: "<<inputpath_data<<std::endl;
    TChain* chain_data = new TChain("FinalTree"+category+"_Bkg");
    chain_data->Add(inputpath_data);
    BDT_decision class_data(chain_data);
    class_data.Loop(isMC, category);

    //MC Ds
    isMC = 1;
    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f1 || !f1->IsOpen()) f1 = new TFile(inputpath_Ds);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    TChain* chain_ds = new TChain("FinalTree"+category+"_sgn");
    chain_ds->Add(inputpath_Ds);
    BDT_decision class_Ds(chain_ds);
    class_Ds.Loop(isMC, category);

//    //MC B0
//    isMC = 2;
//    TFile *f2 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
//    if (!f2 || !f2->IsOpen()) f2 = new TFile(inputpath_B0);
//    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
//    TChain* chain_b0 = new TChain("FinalTree"+category+"_sgn");
//    chain_b0->Add(inputpath_B0);
//    BDT_decision class_B0(chain_b0);
//    class_B0.Loop(isMC, category);
//
//    //MC Bp
//    isMC = 3;
//    TFile *f3 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
//    if (!f3 || !f3->IsOpen()) f3 = new TFile(inputpath_Bp);
//    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;
//    TChain* chain_bp = new TChain("FinalTree"+category+"_sgn");
//    chain_bp->Add(inputpath_Bp);
//    BDT_decision class_Bp(chain_bp);
//    class_Bp.Loop(isMC, category);

    //Add signals to inclusive distribution
    TFile *fout_update = new TFile(fout_path, "update");
    fout_update->cd();
    TH1F *hBDTdecisionDs, *hBDTdecisionB0, *hBDTdecisionBp;
    hBDTdecisionDs = (TH1F*)fout_update->Get("BDTdecision_signalDs"+category);
//    hBDTdecisionB0 = (TH1F*)fout_update->Get("BDTdecision_signalB0"+category);
//    hBDTdecisionBp = (TH1F*)fout_update->Get("BDTdecision_signalBp"+category);

    //inclusive distribution
    TH1F * hBDTdecisionSignal = new TH1F ("BDTdecision_signal"+category, "BDT decision "+category, 240, -0.6, 0.6);
    hBDTdecisionSignal->Add(hBDTdecisionDs);
//    hBDTdecisionSignal->Add(hBDTdecisionB0);
//    hBDTdecisionSignal->Add(hBDTdecisionBp);
    hBDTdecisionSignal->Write();
    fout_update->Close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    return 0;
}
