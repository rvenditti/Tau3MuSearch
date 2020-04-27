//root -l MVA_code_2017.cpp\(\"A\"\)

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


void MVA_code_2017(TString categ){
    //Check on input argument
    if(!categ.Contains("A") && !categ.Contains("B") && !categ.Contains("C")){
        cout << "Please choose between any combination of 'A', 'B' and 'C'" << endl;
        return;
    }

    // Output file
    TFile *fout = new TFile("TMVA_"+TMVA_outputpath+categ+".root", "RECREATE");

    TString cat_name[] = {"A", "B", "C"};
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
    //Bp
    TFile *f_sig_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
    if (!f_sig_bp || !f_sig_bp->IsOpen()) {
       f_sig_bp = new TFile(inputpath_Bp);
    }
    sigTree.push_back( (TTree*)f_sig_bp->Get(treeName));

    //background
    //Loop on run
    for(int j = 0; j<5; j++){
        TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[j]);
        if (!f_bkg || !f_bkg->IsOpen()) {
            f_bkg = new TFile(inputpath_datarun[j]);
        }
        bkgTree.push_back((TTree*)f_bkg->Get(treeName));
    }

    // Set the event weights per tree
    Double_t sigWeight1  = wNormDs/wNormDs; //1.0 Ds
    Double_t sigWeight2  = wNormB0/wNormDs; //B0
    Double_t sigWeight3  = wNormBp/wNormDs; //Bp
    cout<<"Ds sigWeight1: "<<sigWeight1<<" - B0 sigWeight2: "<<sigWeight2<<" - Bp sigWeight3: "<<sigWeight3<<endl;

    Double_t bkgWeight1 = 1.0;
    Double_t bkgWeight2 = 1.0;
    Double_t bkgWeight3 = 1.0;
    Double_t bkgWeight4 = 1.0;
    Double_t bkgWeight5 = 1.0;
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
    
    dataloader->AddSignalTree(sigTree.at(0), sigWeight1);
    dataloader->AddSignalTree(sigTree.at(1), sigWeight2);
    dataloader->AddSignalTree(sigTree.at(2), sigWeight3);
    dataloader->AddBackgroundTree(bkgTree.at(0), bkgWeight1);
    dataloader->AddBackgroundTree(bkgTree.at(1), bkgWeight2);
    dataloader->AddBackgroundTree(bkgTree.at(2), bkgWeight3);
    dataloader->AddBackgroundTree(bkgTree.at(3), bkgWeight4);
    dataloader->AddBackgroundTree(bkgTree.at(4), bkgWeight5);

    // Spectators declaration
    for(int k = 0; k<var_spec_names.size(); k++){
        dataloader->AddSpectator(var_spec_names.at(k), var_spec_names.at(k), "", 'F');
    }

    // Variables declaration
    for(int k = 0; k<var_names.size(); k++){
        dataloader->AddVariable(var_BDT.at(k), var_names.at(k), "", 'F');
    }

    //dataloader->SetSignalWeightExpression( "puFactor" );

    TCut cutS = "tripletMass<2.0 && tripletMass>1.62"; //Signal -> MC full range 
    TCut cutB = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)"; //Background -> data sidebands
    TCut preselCut = "";

    TCut reso_A = "tripletMassReso < 0.007";
    TCut reso_B = "tripletMassReso >= 0.007 && tripletMassReso <= 0.0105";
    TCut reso_C = "tripletMassReso > 0.0105";
    TCut reso_cat = "tripletMassReso < 0"; //always false    

    if(categ.Contains("A")) reso_cat = reso_cat || reso_A;
    if(categ.Contains("B")) reso_cat = reso_cat || reso_B;
    if(categ.Contains("C")) reso_cat = reso_cat || reso_C;

    TString prepareTrainTestOptions = "";
    prepareTrainTestOptions = ":SplitMode=Random"
                              ":NormMode=NumEvents"
                              ":!V";
    dataloader->PrepareTrainingAndTestTree(reso_cat&&preselCut&&cutS, reso_cat&&preselCut&&cutB, prepareTrainTestOptions);

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
