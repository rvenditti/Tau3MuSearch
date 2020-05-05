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

void MVA_code_2018(TString categ){
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
//    //Bp
//    TFile *f_sig_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
//    if (!f_sig_bp || !f_sig_bp->IsOpen()) {
//       f_sig_bp = new TFile(inputpath_Bp);
//    }
//    sigTree.push_back( (TTree*)f_sig_bp->Get(treeName));

    //add TTreeFriend
    TString friendTreeName[] = {"TreeMu1=TreeMu1", "TreeMu2=TreeMu2", "TreeMu3=TreeMu3"};
    for(int i=0; i<3; i++){
        sigTree.at(0)->AddFriend(friendTreeName[i], inputpath_Ds);
        sigTree.at(1)->AddFriend(friendTreeName[i], inputpath_B0);
    //    sigTree.at(2)->AddFriend(friendTreeName[i], inputpath_Bp);
    }
    //background
    //Loop on run
    int n_bkg = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    for(int j = 0; j<n_bkg; j++){
        TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[j]);
        if (!f_bkg || !f_bkg->IsOpen()) {
            f_bkg = new TFile(inputpath_datarun[j]);
        }
        bkgTree.push_back((TTree*)f_bkg->Get(treeName));
    }
    //add TTreeFriend
    for(int i=0; i<3; i++){
        for(int j=0; j<n_bkg; j++){
            bkgTree.at(j)->AddFriend(friendTreeName[i], inputpath_datarun[j]);
        }
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
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
    
    dataloader->AddSignalTree(sigTree.at(0), sigWeight1);
    dataloader->AddSignalTree(sigTree.at(1), sigWeight2);
    //dataloader->AddSignalTree(sigTree.at(2), sigWeight3);
    dataloader->AddBackgroundTree(bkgTree.at(0), bkgWeight1);
    dataloader->AddBackgroundTree(bkgTree.at(1), bkgWeight2);
    dataloader->AddBackgroundTree(bkgTree.at(2), bkgWeight3);
    //dataloader->AddBackgroundTree(bkgTree.at(3), bkgWeight4);

    std::vector<TString> var_train_name;
    std::vector<TString> var_train_def;
    std::vector<TString> var_spec_name;
    std::vector<TString> var_spec_def;
    TString BDTinVar;
    if(categ.Contains("A")) BDTinVar = BDTinVar_A;
    if(categ.Contains("B")) BDTinVar = BDTinVar_B;
    if(categ.Contains("C")) BDTinVar = BDTinVar_C;
    
    readVarName(var_train_name, var_train_def, BDTinVar);
    readVarName(var_spec_name, var_spec_def, BDTspecVar);
 
    // Spectators declaration
    cout<<"Declaration of spectator variables - category "<<categ<<" from file:"<<BDTinVar<<endl;
    for(int k = 0; k<var_spec_name.size(); k++){
        if(var_spec_name.at(k)=="" || var_spec_def.at(k)=="") continue;
        cout<<k<<" - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
        dataloader->AddSpectator(var_spec_def.at(k), var_spec_name.at(k), "", 'D');
    }

    // Variables declaration
    cout<<"Declaration of variables for training - category "<<categ<<" from file:"<<BDTspecVar<<endl;
    for(int k = 0; k<var_train_name.size(); k++){
        if(var_train_name.at(k)=="" || var_train_def.at(k)=="") continue;
        cout<<k<<" - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
        dataloader->AddVariable(var_train_def.at(k), var_train_name.at(k), "", 'D');
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
    dataloader->PrepareTrainingAndTestTree(reso_cat&&preselCut&&cutS&&chi2cut, reso_cat&&preselCut&&cutB&&chi2cut, prepareTrainTestOptions);

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
