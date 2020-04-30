//root -l MVA_muonid.cpp\(\"category\"\)

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
    if(!categ.Contains("barrel") && !categ.Contains("endcap")){
        cout << "Please choose between any combination of 'barrel' and 'endcap'" << endl;
        return;
    }

    // Output file
    TFile *fout = new TFile("TMVA_"+TMVA_outputpath+categ+".root", "RECREATE");

    std::vector<TTree*> sigTree, bkgTree;

    // Get the signal and background trees from TFile source(s);
    //signal
    TString treeName = "FinalTree";
    //Ds
    TFile *f_sig_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
    if (!f_sig_ds || !f_sig_ds->IsOpen()) {
         f_sig_ds = new TFile(inputpath_Ds);
    }
    sigTree.push_back( (TTree*)f_sig_ds->Get(treeName));
    //B0
    TFile *f_sig_b0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
    if (!f_sig_b0 || !f_sig_b0->IsOpen()) {
        f_sig_b0 = new TFile(inputpath_B0);
    }
    sigTree.push_back( (TTree*)f_sig_b0->Get(treeName));
    ////Bp
    //TFile *f_sig_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    //if (!f_sig_bp || !f_sig_bp->IsOpen()) {
    //   f_sig_bp = new TFile(inputpath_Bp);
    //}
    //sigTree.push_back( (TTree*)f_sig_bp->Get(treeName_sig));

    //background
    size_t n_bkg = sizeof(inputpath_bkg)/sizeof(inputpath_bkg[0]);
    //Loop on samples
    for(int j = 0; j<n_bkg; j++){
        TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_bkg[j]);
        if (!f_bkg || !f_bkg->IsOpen()) {
            f_bkg = new TFile(inputpath_bkg[j]);
        }
        bkgTree.push_back((TTree*)f_bkg->Get(treeName));
    }

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
    for(int k = 0; k<var_train_names.size(); k++){
        dataloader->AddVariable(var_train_def.at(k), var_train_names.at(k), "", 'F');
    }


    //signal muons
    TCut cutS = "abs(mu_simPdgId)==13 && abs(mu_simMotherPdgId)==15"; 
    //bkg muons = pions, kaon, muons frmo decay in flight
    TCut cutB_pi = "( abs(mu_simPdgId) == 211 )"; //pion
    TCut cutB_k = "( abs(mu_simPdgId) == 321 )"; //kaon
    TCut cutB_mu_from_pi = "( abs(mu_simPdgId) == 13 && abs(mu_simMotherPdgId) == 211 )"; //true mu from pion decay
    TCut cutB_mu_from_k = "( abs(mu_simPdgId) == 13 && abs(mu_simMotherPdgId) == 321 )"; //true mu from kaon decay
    //running on global muons with associated simInfo
    TCut preselCut = "mu_simType != 0 && mu_isGlobal == 1 && mu_pt > 2 && abs(mu_eta)<2.4";
    //loose preselection cuts on input variables
    TCut cleanInputVar = "mu_combinedQuality_staRelChi2 > 0 &&"
                         "mu_combinedQuality_trkRelChi2 > 0 &&"
                         "mu_GLhitPattern_numberOfValidMuonHits > 0 &&"
                         "mu_validMuonHitComb > 0 &&"
                         "mu_segmentCompatibility > 0.05";

    //training separately on endcap and barrel
    TCut etarange = "abs(mu_eta)<0"; //always false
    TCut barrel = "abs(mu_eta)<1.2";
    TCut endcap = "abs(mu_eta)>=1.2";
    if(categ.Contains("endcap")) etarange = etarange || endcap;
    if(categ.Contains("barrel")) etarange = etarange || barrel;

    dataloader->SetBackgroundWeightExpression("ptetaWeight");

    TString prepareTrainTestOptions = ":SplitMode=Random"
                                       ":NormMode=NumEvents"
                                       ":!V";
    dataloader->PrepareTrainingAndTestTree(preselCut&&cleanInputVar&&etarange&&cutS, preselCut&&cleanInputVar&&etarange&&(cutB_pi||cutB_k||cutB_mu_from_pi||cutB_mu_from_k), prepareTrainTestOptions);
    
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
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}
