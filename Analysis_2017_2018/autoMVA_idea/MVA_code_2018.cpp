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
using namespace std;

void MVA_code_2018(TString categ){
    //Check on input argument
    if(!categ.Contains("A") && !categ.Contains("B") && !categ.Contains("C")){
        cout << "Please choose between any combination of 'A', 'B' and 'C'" << endl;
        return;
    }

    // Output file
    TFile *fout = new TFile("TMVA_"+TMVA_outputpath+categ+".root", "RECREATE");

    TString cat_name[] = {"A", "B", "C"};
    vector<TTree*> sigTree, bkgTree;

    // Get the signal and background trees from TFile source(s);
    // SIGNAL tree
    for(int i = 0; i<3; i++){
       //signal
       TString treeName_sig = "FinalTree";
       //Ds
       TFile *f_sig_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
       if (!f_sig_ds || !f_sig_ds->IsOpen()) {
            f_sig_ds = new TFile(inputpath_Ds);
       }
       sigTree.push_back( (TTree*)f_sig_ds->Get(treeName_sig)); //A:0 B:3 C:6
       //B0
       TFile *f_sig_b0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
       if (!f_sig_b0 || !f_sig_b0->IsOpen()) {
           f_sig_b0 = new TFile(inputpath_B0);
       }
       sigTree.push_back( (TTree*)f_sig_b0->Get(treeName_sig)); //A:1 B:4 C:7
       //Bp
       TFile *f_sig_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
       if (!f_sig_bp || !f_sig_bp->IsOpen()) {
          f_sig_bp = new TFile(inputpath_Bp);
       }
       sigTree.push_back( (TTree*)f_sig_bp->Get(treeName_sig)); //A:2 B:5 C:8

       //background
       TString treeName_bkg = "FinalTree";
       //Loop on run
       for(int j = 0; j<1; j++){
           TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[j]);
           if (!f_bkg || !f_bkg->IsOpen()) {
               f_bkg = new TFile(inputpath_datarun[j]);
           }
           bkgTree.push_back((TTree*)f_bkg->Get(treeName_bkg)); //A: 0 1 2 3  B: 4 5 6 7  C: 8 9 10 11
       }
    }
    // Set the event weights per tree
    Double_t sigWeight1  = 1.0; //Ds
    Double_t sigWeight2  = 0.76; //B0
    Double_t sigWeight3  = 0.35; //Bp

    Double_t bkgWeight1 = 1.0;
    Double_t bkgWeight2 = 1.0;
    Double_t bkgWeight3 = 1.0;
    Double_t bkgWeight4 = 1.0;
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
    
    if(categ.Contains("A")){
        dataloader->AddSignalTree(sigTree.at(0), sigWeight1);
        dataloader->AddSignalTree(sigTree.at(1), sigWeight2);
        dataloader->AddSignalTree(sigTree.at(2), sigWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(0), bkgWeight1);
//        dataloader->AddBackgroundTree(bkgTree.at(1), bkgWeight2);
//        dataloader->AddBackgroundTree(bkgTree.at(2), bkgWeight3);
//        dataloader->AddBackgroundTree(bkgTree.at(3), bkgWeight4);
    }
    if(categ.Contains("B")){
        dataloader->AddSignalTree(sigTree.at(3), sigWeight1);
        dataloader->AddSignalTree(sigTree.at(4), sigWeight2);
        dataloader->AddSignalTree(sigTree.at(5), sigWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(4), bkgWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(5), bkgWeight2);
        dataloader->AddBackgroundTree(bkgTree.at(6), bkgWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(7), bkgWeight4);
    }
    if(categ.Contains("C")){
        dataloader->AddSignalTree(sigTree.at(6), sigWeight1);
        dataloader->AddSignalTree(sigTree.at(7), sigWeight2);
        dataloader->AddSignalTree(sigTree.at(8), sigWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(8), bkgWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(9), bkgWeight2);
        dataloader->AddBackgroundTree(bkgTree.at(10), bkgWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(11), bkgWeight4);
    }

    // Import BDT input variables
    ImportVar(VarName, VarSel, VarType);
    // Spectators declaration
    dataloader->AddSpectator("tripletMass", 'F'); // triplet invariant mass
    dataloader->AddSpectator("puFactor", 'F'); // event-by-event PU weight //1 for data
    dataloader->AddSpectator("isGlb2", 'I');
    dataloader->AddSpectator("isGlb3", 'I');
    // Variables declaration
    for(int j=0; j<VarName.size()-1; j++) dataloader->AddVariable(VarSel.at(j), VarName.at(j), "", VarType.at(j));

    dataloader->SetSignalWeightExpression( "puFactor" );

    TCut cutS = "tripletMass<2.0 && tripletMass>1.62"; //Signal -> MC full range 
    TCut cutB = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)"; //Background -> data sidebands
    TCut preselCut = "isGlb2==1 && isGlb3==1";

    //TCut preselCut = "";
    TString prepareTrainTestOptions = "";
    if(doCV) prepareTrainTestOptions = "nTest_Signal=1"
                                       ":nTest_Background=1"
                                       ":SplitMode=Random"
                                       ":NormMode=NumEvents"
                                       ":!V";
    else     prepareTrainTestOptions = ":SplitMode=Random"
                                       ":NormMode=NumEvents"
                                       ":!V";
    dataloader->PrepareTrainingAndTestTree(preselCut&&cutS, preselCut&&cutB, prepareTrainTestOptions);

    Int_t numFolds = 5;
    TString analysisType = "Classification";
    TString splitExpr = "";
    //TString outputEnsembling = "Avg";
    TString foldFileOutput = "False";
    TString cvOptions = Form("!V"
                           ":!Silent"
                           ":ModelPersistence"
                           ":AnalysisType=%s"
                           ":NumFolds=%i"
                           ":SplitExpr=%s"
      //                     ":OutputEnsembling=%s"
                           ":FoldFileOutput=%s",
                           analysisType.Data(), numFolds,
                           splitExpr.Data(), 
                           //outputEnsembling.Data(), 
                           foldFileOutput.Data());
    

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
