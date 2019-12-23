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

#include "TMVA/CrossValidation.h"
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
    // SIGNAL tree
    for(int i = 0; i<3; i++){
       //signal
       TString treeName_sig = "FinalTree"+cat_name[i]+"_sgn";
       //Ds
       TFile *f_sig_ds = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Ds);
       if (!f_sig_ds || !f_sig_ds->IsOpen()) {
            f_sig_ds = new TFile(inputpath_Ds);
       }
       sigTree.push_back( (TTree*)f_sig_ds->Get(treeName_sig)); //A:0 B:1 C:2
     //  //B0
     //  TFile *f_sig_b0 = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_B0);
     //  if (!f_sig_b0 || !f_sig_b0->IsOpen()) {
     //      f_sig_b0 = new TFile(inputpath_B0);
     //  }
     //  sigTree.push_back( (TTree*)f_sig_b0->Get(treeName_sig));
     //  //Bp
     //  TFile *f_sig_bp = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_Bp);
     //  if (!f_sig_bp || !f_sig_bp->IsOpen()) {
     //     f_sig_bp = new TFile(inputpath_Bp);
     //  }
     //  sigTree.push_back( (TTree*)f_sig_bp->Get(treeName_sig));

       //background
       TString treeName_bkg = "FinalTree"+cat_name[i]+"_Bkg";
       //Loop on run
       for(int j = 0; j<4; j++){
           TFile *f_bkg = (TFile*)gROOT->GetListOfFiles()->FindObject(inputpath_datarun[j]);
           if (!f_bkg || !f_bkg->IsOpen()) {
               f_bkg = new TFile(inputpath_datarun[j]);
           }
           bkgTree.push_back((TTree*)f_bkg->Get(treeName_bkg)); //A: 0 1 2 3  B: 4 5 6 7  C: 8 9 10 11
       }
    }
    // Set the event weights per tree
    Double_t sigWeight1  = 1.0; //Ds
    //Double_t sigWeight2  = 1.0; //B0
    //Double_t sigWeight3  = 1.0; //Bp
    Double_t bkgWeight1 = 1.0;
    Double_t bkgWeight2 = 1.0;
    Double_t bkgWeight3 = 1.0;
    Double_t bkgWeight4 = 1.0;
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader(TMVA_outputpath+categ);
    
    if(categ.Contains("A")){
        dataloader->AddSignalTree(sigTree.at(0), sigWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(0), bkgWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(1), bkgWeight2);
        dataloader->AddBackgroundTree(bkgTree.at(2), bkgWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(3), bkgWeight4);
    }
    if(categ.Contains("B")){
        dataloader->AddSignalTree(sigTree.at(1), sigWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(4), bkgWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(5), bkgWeight2);
        dataloader->AddBackgroundTree(bkgTree.at(6), bkgWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(7), bkgWeight4);
    }
    if(categ.Contains("C")){
        dataloader->AddSignalTree(sigTree.at(2), sigWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(8), bkgWeight1);
        dataloader->AddBackgroundTree(bkgTree.at(9), bkgWeight2);
        dataloader->AddBackgroundTree(bkgTree.at(10), bkgWeight3);
        dataloader->AddBackgroundTree(bkgTree.at(11), bkgWeight4);
    }

    // Spectators declaration
    dataloader->AddSpectator("tripletMass", 'F'); // triplet invariant mass

    // Variables declaration
    dataloader->AddVariable("Pmu3", 'F');
    dataloader->AddVariable("cLP", 'F');
    dataloader->AddVariable("tKink", 'F');
    dataloader->AddVariable("segmComp", 'F'); // segmComp MIN
    dataloader->AddVariable("fv_nC", 'F'); // Chi2 norm vertex
    dataloader->AddVariable("fv_dphi3D", 'F'); // Angle PVSV
    dataloader->AddVariable("fv_d3Dsig", 'F'); // Flight dist significance
//    dataloader->AddVariable("d0", 'F');
    dataloader->AddVariable("d0sig", 'F'); // Transverse IP significance
    dataloader->AddVariable("mindca_iso", 'F'); // isolation MIN Distance of Closest Approach
//    dataloader->AddVariable("caloComp", 'F'); // calo Comp MIN
    dataloader->AddVariable("trkRel", 'F');
//    dataloader->AddVariable("Pmu1", 'F');
//    dataloader->AddVariable("Ptmu1", 'F');
//    dataloader->AddVariable("Etamu1", 'F');
//    dataloader->AddVariable("Pmu2", 'F');
//    dataloader->AddVariable("Ptmu2", 'F');
//    dataloader->AddVariable("Etamu2", 'F');
//    dataloader->AddVariable("Ptmu3", 'F');
//    dataloader->AddVariable("Etamu3", 'F');
//    dataloader->AddVariable("P_tripl", 'F');
//    dataloader->AddVariable("Pt_tripl", 'F');
//    dataloader->AddVariable("Eta_tripl", 'F');
//    dataloader->AddVariable("nStMu1", 'F');
//    dataloader->AddVariable("nStMu2", 'F');
//    dataloader->AddVariable("nStMu3", 'F'); //number of matched stations
//    dataloader->AddVariable("Iso03Mu1", 'F');
//    dataloader->AddVariable("Iso03Mu2", 'F');
//    dataloader->AddVariable("Iso03Mu3", 'F');
//    dataloader->AddVariable("Iso05Mu1", 'F');
//    dataloader->AddVariable("Iso05Mu2", 'F');
//    dataloader->AddVariable("Iso05Mu3", 'F');
//    dataloader->AddVariable("nMatchesMu1", 'F');
//    dataloader->AddVariable("nMatchesMu2", 'F');
     dataloader->AddVariable("nMatchesMu3", 'F'); //number of matched segments
//    dataloader->AddVariable("timeAtIpInOut1", 'F');
//    dataloader->AddVariable("timeAtIpInOut2", 'F');
//    dataloader->AddVariable("timeAtIpInOut3", 'F');
//    dataloader->AddVariable("cQ_uS", 'F');
//    dataloader->AddVariable("cQ_tK", 'F');
//    dataloader->AddVariable("cQ_gK", 'F');
//    dataloader->AddVariable("cQ_tRChi2", 'F');
//    dataloader->AddVariable("cQ_sRChi2", 'F');
//    dataloader->AddVariable("cQ_Chi2LM", 'F');
//    dataloader->AddVariable("cQ_Chi2lD", 'F');
//    dataloader->AddVariable("cQ_gDEP", 'F');
//    dataloader->AddVariable("cQ_tM", 'F');
//    dataloader->AddVariable("cQ_gTP", 'F');
//    dataloader->AddVariable("calEn_emMu1", 'F');
//    dataloader->AddVariable("calEn_emMu2", 'F');
//    dataloader->AddVariable("calEn_emMu3", 'F');
//    dataloader->AddVariable("calEn_hadMu1", 'F');
//    dataloader->AddVariable("calEn_hadMu2", 'F');
//    dataloader->AddVariable("calEn_hadMu3", 'F');
//    dataloader->AddVariable("fliDistPVSV_Chi2", 'F');
//    dataloader->AddVariable("isGlb3", 'F');
//    dataloader->AddVariable("isTracker3", 'F');
//    dataloader->AddVariable("isLoose3", 'F');
//    dataloader->AddVariable("isSoft3", 'F');
//    dataloader->AddVariable("isPF3", 'F');
//    dataloader->AddVariable("isRPC3", 'F');
//    dataloader->AddVariable("isSA3", 'F');
//    dataloader->AddVariable("isCalo3", 'F');
//    dataloader->AddVariable("Vx1", 'F');
//    dataloader->AddVariable("Vx2", 'F');
//    dataloader->AddVariable("Vx3", 'F');
//    dataloader->AddVariable("Vy1", 'F');
//    dataloader->AddVariable("Vy2", 'F');
//    dataloader->AddVariable("Vy3", 'F');
//    dataloader->AddVariable("Vz1", 'F');
//    dataloader->AddVariable("Vz2", 'F');
//    dataloader->AddVariable("Vz3", 'F');
//    dataloader->AddVariable("RefVx1", 'F');
//    dataloader->AddVariable("RefVy1", 'F');
//    dataloader->AddVariable("RefVz1", 'F');
//    dataloader->AddVariable("SVx", 'F');
//    dataloader->AddVariable("SVy", 'F');
//    dataloader->AddVariable("SVz", 'F');
//    dataloader->AddVariable("had03", 'F');
//    dataloader->AddVariable("had05", 'F');
//    dataloader->AddVariable("nJets03", 'F');
//    dataloader->AddVariable("nJets05", 'F');
//    dataloader->AddVariable("nTracks03", 'F');
//    dataloader->AddVariable("nTracks05", 'F');
//    dataloader->AddVariable("sumPt03", 'F');
//    dataloader->AddVariable("sumPt05", 'F');
//    dataloader->AddVariable("hadVeto03", 'F');
//    dataloader->AddVariable("hadVeto05", 'F');
//    dataloader->AddVariable("emVeto03", 'F');
//    dataloader->AddVariable("emVeto05", 'F');
//    dataloader->AddVariable("trVeto03", 'F');
//    dataloader->AddVariable("trVeto05", 'F');
    

    TCut cutS = "tripletMass<1.82 && tripletMass>1.73"; //Signal -> MC peak region 
    TCut cutB = "(tripletMass<1.73 && tripletMass>1.65) || (tripletMass<1.90 && tripletMass>1.82)"; //Background -> data sidebands
   // TCut preselCut = "Pmu3<30 && cLP>0 && cLP<500 && tKink>0 && tKink<400 && segmComp>0.05 && fv_nC<100 && fv_dphi3D<0.4 && fv_d3Dsig>0 && fv_d3Dsig<150 && d0sig>0 && d0sig<20 && mindca_iso>0 && mindca_iso<0.6 && trkRel>0 && trkRel<5";
    TCut preselCut = "";
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
    
    TMVA::CrossValidation cv{"TMVACrossValidation", dataloader, fout, cvOptions};

    if(doCV){
        // Booking of MVA methods : BDT
        cv.BookMethod( TMVA::Types::kBDT, "BDT",
             "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50");
             //"!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=50" );

        cv.Evaluate();

        size_t iMethod = 0;
        for (auto && result : cv.GetResults()) {
            std::cout << "Summary for method " << cv.GetMethods()[iMethod++].GetValue<TString>("MethodName")
                    << std::endl;
            for (UInt_t iFold = 0; iFold<cv.GetNumFolds(); ++iFold) {
                std::cout << "\tFold " << iFold << ": "
                       << "ROC int: " << result.GetROCValues()[iFold]
                       << ", "
                       << "BkgEff@SigEff=0.3: " << result.GetEff30Values()[iFold]
                       << std::endl;
            }
        }
    }
    else {
        // Booking of MVA methods : BDT
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", 
             "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=50");
             //"!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=50" );

        // Booking of MVA methods : MLP
        //factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
 
        // Training the MVA methods
        factory->TrainAllMethods();
        
        // Testing the MVA methods
        factory->TestAllMethods();
        
        // Evaluating the MVA methods
        factory->EvaluateAllMethods();
    } 

    // Save the output
    fout->Close();
    
    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()){
        TMVAGui("TMVA_"+TMVA_outputpath+categ+".root");
        if(doCV){
            // Draw cv-specific graphs
            // cv.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for BDTG");
            for(int i=0; i<numFolds; i++){
                TString s = std::to_string(i);
               // TMVAGui("BDT_fold"+s+".root");
            }
        }
    }
}
