//root -l MVA_code_2018.cpp\(\"A\"\)

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

//#include "TMVA/CrossValidation.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

#include "T3M_common.h"

using namespace TMVA;


void MVA_muonid(TString categ){
    //Check on input argument
    //if(!categ.Contains("A") && !categ.Contains("B") && !categ.Contains("C")){
    //    cout << "Please choose between any combination of 'A', 'B' and 'C'" << endl;
    //    return;
    //}

    // Output file
    TFile *fout = new TFile("TMVA_"+TMVA_outputpath+categ+".root", "RECREATE");

    std::vector<TTree*> sigTree, bkgTree;

    // Get the signal and background trees from TFile source(s);
    // SIGNAL tree
    //signal
    TString treeName_sig = "FinalTree";
    //Ds
    TFile *f_sig_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f_sig_ds || !f_sig_ds->IsOpen()) {
         f_sig_ds = new TFile(inputpath_Ds);
    }
    sigTree.push_back( (TTree*)f_sig_ds->Get(treeName_sig));
    ////B0
    //TFile *f_sig_b0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
    //if (!f_sig_b0 || !f_sig_b0->IsOpen()) {
    //    f_sig_b0 = new TFile(inputpath_B0);
    //}
    //sigTree.push_back( (TTree*)f_sig_b0->Get(treeName_sig));
    ////Bp
    //TFile *f_sig_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    //if (!f_sig_bp || !f_sig_bp->IsOpen()) {
    //   f_sig_bp = new TFile(inputpath_Bp);
    //}
    //sigTree.push_back( (TTree*)f_sig_bp->Get(treeName_sig));

    //background
    TString treeName_bkg = "FinalTree";
    //BdToPiPi
    TFile *f_bkg1 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_BdToPiPi);
    if (!f_bkg1 || !f_bkg1->IsOpen()) {
        f_bkg1 = new TFile(inputpath_BdToPiPi);
    }
    bkgTree.push_back((TTree*)f_bkg1->Get(treeName_bkg));
    //BdToKK
    TFile *f_bkg2 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_BdToKK);
    if (!f_bkg2 || !f_bkg2->IsOpen()) {
        f_bkg2 = new TFile(inputpath_BdToKK);
    }
    bkgTree.push_back((TTree*)f_bkg2->Get(treeName_bkg));

    // Set the event weights per tree
    Double_t sigWeight  = 1.0;
    Double_t bkgWeight = 1.0;
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
   
    for(int i = 0; i<sigTree.size(); i++) {
        dataloader->AddSignalTree(sigTree.at(i), sigWeight);
    }
    for(int j = 0; j<bkgTree.size(); j++) {
        dataloader->AddBackgroundTree(bkgTree.at(j), bkgWeight);
    }

    // Spectators declaration
    for(int k = 0; k<var_spec_names.size(); k++){
        dataloader->AddSpectator(var_spec_names.at(k), var_spec_names.at(k), "", 'F');
    }

    // Variables declaration
    for(int k = 0; k<var_names.size(); k++){
        dataloader->AddVariable(var_names.at(k), var_names.at(k), "", 'F');
    }

    TCut cutS = "abs(mu_simPdgId)==13 && abs(mu_simMotherPdgId)==15"; 
    TCut cutB = "";
    TCut preselCut = "mu_isGlobal == 1 && mu_pt > 2 && abs(mu_eta)<2.4";

    //TCut preselCut = "";
    TString prepareTrainTestOptions = ":SplitMode=Random"
                                       ":NormMode=NumEvents"
                                       ":!V";
    dataloader->PrepareTrainingAndTestTree(preselCut&&cutS, preselCut&&cutB, prepareTrainTestOptions);
    
    // Booking of MVA methods : BDT
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", 
         "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50" );
         //"!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=50" );

    // Booking of MVA methods : MLP
    //factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
 
    // Training the MVA methods
    factory->TrainAllMethods();
    
    // Testing the MVA methods
    factory->TestAllMethods();
    
    // Evaluating the MVA methods
    factory->EvaluateAllMethods();

    // Save the output
    fout->Close();
    
    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()){
        TMVAGui("TMVA_"+TMVA_outputpath+categ+".root");
    }
}
