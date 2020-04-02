#include "RunT3Mu.C"
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
    if (!(category == "A" || category == "B" || category == "C"|| category == "B+C")) {
        std::cout<<"Wrong category name, please chooose between: A B C"<<std::endl;
        return 0;
    }
    
    //prepare output file
    TString fout_path = "datacardT3Mu_"+TMVA_inputpath+category+".root";
    TFile *fout = new TFile(fout_path, "recreate");
    fout->cd();
    fout->Close();

    //data
    isMC = 0;
    TFile *f0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_data);
    if (!f0 || !f0->IsOpen()) f0 = new TFile(inputpath_data);
    std::cout<<"Opened input file: "<<inputpath_data<<std::endl;
    if(category == "B+C") {
        TChain* chain_data = new TChain("FinalTreeB_sgn");
        chain_data->Add(inputpath_data);
        chain_data->Add(inputpath_data+"/FinalTreeC_sgn");
        RunT3Mu class_data(chain_data);
        class_data.Loop(isMC, category);
    }
    else{
        TChain* chain_data = new TChain("FinalTree"+category+"_sgn");
        chain_data->Add(inputpath_data);
        RunT3Mu class_data(chain_data);
        class_data.Loop(isMC, category);
    }

    //MC Ds
    isMC = 1;
    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f1 || !f1->IsOpen()) f1 = new TFile(inputpath_Ds);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    if(category == "B+C") {
        TChain* chain_ds = new TChain("FinalTreeB_sgn");
        chain_ds->Add(inputpath_Ds);
        chain_ds->Add(inputpath_Ds+"/FinalTreeC_sgn");
        RunT3Mu class_Ds(chain_ds);
        class_Ds.Loop(isMC, category);
    }
    else{
        TChain* chain_ds = new TChain("FinalTree"+category+"_sgn");
        chain_ds->Add(inputpath_Ds);
        RunT3Mu class_Ds(chain_ds);
        class_Ds.Loop(isMC, category);
    }
    //MC B0
    isMC = 2;
    TFile *f2 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
    if (!f2 || !f2->IsOpen()) f2 = new TFile(inputpath_B0);
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    if(category == "B+C") {
        TChain* chain_b0 = new TChain("FinalTreeB_sgn");
        chain_b0->Add(inputpath_B0);
        chain_b0->Add(inputpath_B0+"/FinalTreeC_sgn");
        RunT3Mu class_B0(chain_b0);
        class_B0.Loop(isMC, category);
    }
    else{
        TChain* chain_b0 = new TChain("FinalTree"+category+"_sgn");
        chain_b0->Add(inputpath_B0);
        RunT3Mu class_B0(chain_b0);
        class_B0.Loop(isMC, category);
    }
    //MC Bp
    isMC = 3;
    TFile *f3 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    if (!f3 || !f3->IsOpen()) f3 = new TFile(inputpath_Bp);
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;
    if(category == "B+C") {
        TChain* chain_bp = new TChain("FinalTreeB_sgn");
        chain_bp->Add(inputpath_Bp);
        chain_bp->Add(inputpath_Bp+"/FinalTreeC_sgn");
        RunT3Mu class_Bp(chain_bp);
        class_Bp.Loop(isMC, category);
    }
    else{
        TChain* chain_bp = new TChain("FinalTree"+category+"_sgn");
        chain_bp->Add(inputpath_Bp);
        RunT3Mu class_Bp(chain_bp);
        class_Bp.Loop(isMC, category);
    }

    //Add signals to inclusive distribution
    TFile *fout_update = new TFile(fout_path, "update");
    fout_update->cd();
    TH1F *hTripletMassDs1, *hTripletMassB01, *hTripletMassBp1;
    TH1F *hTripletMassDs2, *hTripletMassB02, *hTripletMassBp2;
    hTripletMassDs1 = (TH1F*)fout_update->Get("signalDs"+category+"1");
    hTripletMassDs2 = (TH1F*)fout_update->Get("signalDs"+category+"2");
    hTripletMassB01 = (TH1F*)fout_update->Get("signalB0"+category+"1");
    hTripletMassB02 = (TH1F*)fout_update->Get("signalB0"+category+"2");
    hTripletMassBp1 = (TH1F*)fout_update->Get("signalBp"+category+"1");
    hTripletMassBp2 = (TH1F*)fout_update->Get("signalBp"+category+"2");
    //inclusive distribution
    TH1F * hSignal1 = new TH1F ("signal"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
    TH1F * hSignal2 = new TH1F ("signal"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);
    hSignal1->Add(hTripletMassDs1);
    hSignal1->Add(hTripletMassB01);
    hSignal1->Add(hTripletMassBp1);
    hSignal2->Add(hTripletMassDs2);
    hSignal2->Add(hTripletMassB02);
    hSignal2->Add(hTripletMassBp2);
    hSignal1->Write();
    hSignal2->Write();
    fout_update->Close();
    cout<<"++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    return 0;
}
