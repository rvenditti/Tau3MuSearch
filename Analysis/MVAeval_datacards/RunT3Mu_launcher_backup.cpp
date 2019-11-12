#include "RunT3Mu.C"
#include <TROOT.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include <cmath>
#include <fstream>

using namespace std;

int main(int narg, char** cat){
    TTree *tree;
    int isMC;
    TString category = cat[1];
    std::cout << "category: " << category << std::endl;
    
    // Check input arguments
    if (!(category == "A" || category == "B" || category == "C" || category == "B+C")) {
        std::cout<<"Wrong category name, please chooose between: A B C B+C"<<std::endl;
        return 0;
    }

    TString inputpath_data = "AnalysedTree_data_2017_merged_tau3mu_29oct.root";
    TString inputpath_Ds = "AnalysedTree_MC_3globalpt2_Ds_tau3mu_24oct.root";
    TString inputpath_B0 = "AnalysedTree_MC_3globalpt2_B0_tau3mu_24oct.root";
    TString inputpath_Bp = "AnalysedTree_MC_3globalpt2_Bp_tau3mu_24oct.root";

    //prepare output file
    TFile *fout = new TFile("datacardT3Mu_"+category+".root","recreate");
    fout->cd();
    fout->Close();

    //data
    isMC = 0;
    TFile *f0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_data);
    if (!f0 || !f0->IsOpen()) f0 = new TFile(inputpath_data);
    std::cout<<"Opened input file: "<<inputpath_data<<std::endl;
    if(category=="B+C"){
        TTree *treeB;
        TDirectory * dir0B = (TDirectory*)f0->Get(inputpath_data+":/FinalTreeB_Bkg");
        dir0B->GetObject("FinalTreeB_Bkg",treeB);
        RunT3Mu class_data_B(treeB);
        class_data_B.Loop(isMC, category);

        TTree *treeC;
        TDirectory * dir0C = (TDirectory*)f0->Get(inputpath_data+":/FinalTreeC_Bkg");
        dir0C->GetObject("FinalTreeC_Bkg",treeC);
        RunT3Mu class_data_C(treeC);
        class_data_C.Loop(isMC, category);
    }
    else{
        TDirectory * dir0 = (TDirectory*)f0->Get(inputpath_data+":/FinalTree"+category+"_Bkg");
        dir0->GetObject("FinalTree"+category+"_Bkg",tree);
        RunT3Mu class_data(tree);
        class_data.Loop(isMC, category);
    }

    //MC Ds
    isMC = 1;
    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f1 || !f1->IsOpen()) f1 = new TFile(inputpath_Ds);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    if(category=="B+C"){
        TTree *treeB;
        TDirectory * dir1B = (TDirectory*)f1->Get(inputpath_Ds+":/FinalTreeB_sgn");
        dir1B->GetObject("FinalTreeB_sgn",treeB);
        RunT3Mu class_Ds_B(treeB);
        class_Ds_B.Loop(isMC, category);

        TTree *treeC;
        TDirectory * dir1C = (TDirectory*)f1->Get(inputpath_Ds+":/FinalTreeC_sgn");
        dir1C->GetObject("FinalTreeC_sgn",treeC);
        RunT3Mu class_Ds_C(treeC);
        class_Ds_C.Loop(isMC, category);
    }
    else{
        TDirectory * dir1 = (TDirectory*)f1->Get(inputpath_Ds+":/FinalTree"+category+"_sgn");
        dir1->GetObject("FinalTree"+category+"_sgn",tree);
        RunT3Mu class_Ds(tree);
        class_Ds.Loop(isMC, category);
    }
    //MC B0
    isMC = 2;
    TFile *f2 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
    if (!f2 || !f2->IsOpen()) f2 = new TFile(inputpath_B0);
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    if(category=="B+C"){
        TTree *treeB;
        TDirectory * dir2B = (TDirectory*)f2->Get(inputpath_B0+":/FinalTreeB_sgn");
        dir2B->GetObject("FinalTreeB_sgn",treeB);
        RunT3Mu class_B0_B(treeB);
        class_B0_B.Loop(isMC, category);

        TTree *treeC;
        TDirectory * dir2C = (TDirectory*)f2->Get(inputpath_B0+":/FinalTreeC_sgn");
        dir2C->GetObject("FinalTreeC_sgn",treeC);
        RunT3Mu class_B0_C(treeC);
        class_B0_C.Loop(isMC, category);
    }
    else{
        TDirectory * dir2 = (TDirectory*)f2->Get(inputpath_B0+":/FinalTree"+category+"_sgn");
        dir2->GetObject("FinalTree"+category+"_sgn",tree);
        RunT3Mu class_B0(tree);
        class_B0.Loop(isMC, category);
    }
    //MC Bp
    isMC = 3;
    TFile *f3 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    if (!f3 || !f3->IsOpen()) f3 = new TFile(inputpath_Bp);
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;
    if(category=="B+C"){
        TTree *treeB;
        TDirectory * dir3B = (TDirectory*)f3->Get(inputpath_Bp+":/FinalTreeB_sgn");
        dir3B->GetObject("FinalTreeB_sgn",treeB);
        RunT3Mu class_Bp_B(treeB);
        class_Bp_B.Loop(isMC, category);

        TTree *treeC;
        TDirectory * dir3C = (TDirectory*)f3->Get(inputpath_Bp+":/FinalTreeC_sgn");
        dir3C->GetObject("FinalTreeC_sgn",treeC);
        RunT3Mu class_Bp_C(treeC);
        class_Bp_C.Loop(isMC, category);
    }
    else{
        TDirectory * dir3 = (TDirectory*)f3->Get(inputpath_Bp+":/FinalTree"+category+"_sgn");
        dir3->GetObject("FinalTree"+category+"_sgn",tree);
        RunT3Mu class_Bp(tree);
        class_Bp.Loop(isMC, category);
    }

    //Add signals to inclusive distribution
    TFile *fout_update = new TFile("datacardT3Mu_"+category+".root","update");
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
    return 0;
}
