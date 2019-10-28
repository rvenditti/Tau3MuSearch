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
    if (!(category == "A" || category == "B" || category == "C")) {
        std::cout<<"Wrong category name, please chooose between: A B C"<<std::endl;
        return 0;
    }

    TString inputpath_data = "AnalysedTree_data_3globalpt2_2017_merged_tau3mu_24oct.root";
    TString inputpath_Ds = "AnalysedTree_MC_3globalpt2_Ds_tau3mu_24oct.root";
    TString inputpath_B0 = "AnalysedTree_MC_3globalpt2_B0_tau3mu_24oct.root";
    TString inputpath_Bp = "AnalysedTree_MC_3globalpt2_Bp_tau3mu_24oct.root";

    //data
    isMC = 0;
    TFile *f0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_data);
    if (!f0 || !f0->IsOpen()) f0 = new TFile(inputpath_data);
    std::cout<<"Opened input file: "<<inputpath_data<<std::endl;
    TDirectory * dir0 = (TDirectory*)f0->Get(inputpath_data+":/FinalTree"+category+"_Bkg");
    dir0->GetObject("FinalTree"+category+"_Bkg",tree);
    RunT3Mu class_data(tree);
    class_data.Loop(isMC, category);

    //MC Ds
    isMC = 1;
    TFile *f1 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f1 || !f1->IsOpen()) f1 = new TFile(inputpath_Ds);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    TDirectory * dir1 = (TDirectory*)f1->Get(inputpath_Ds+":/FinalTree"+category+"_sgn");
    dir1->GetObject("FinalTree"+category+"_sgn",tree);
    RunT3Mu class_Ds(tree);
    class_Ds.Loop(isMC, category);
    //MC B0
    isMC = 2;
    TFile *f2 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
    if (!f2 || !f2->IsOpen()) f2 = new TFile(inputpath_B0);
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    TDirectory * dir2 = (TDirectory*)f2->Get(inputpath_B0+":/FinalTree"+category+"_sgn");
    dir2->GetObject("FinalTree"+category+"_sgn",tree);
    RunT3Mu class_B0(tree);
    class_B0.Loop(isMC, category);
    //MC Bp
    isMC = 3;
    TFile *f3 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    if (!f3 || !f3->IsOpen()) f3 = new TFile(inputpath_Bp);
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;
    TDirectory * dir3 = (TDirectory*)f3->Get(inputpath_Bp+":/FinalTree"+category+"_sgn");
    dir3->GetObject("FinalTree"+category+"_sgn",tree);
    RunT3Mu class_Bp(tree);
    class_Bp.Loop(isMC, category);

    return 0;
}
