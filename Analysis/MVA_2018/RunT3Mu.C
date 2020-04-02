#define RunT3Mu_cxx
#include "RunT3Mu.h"
#include "T3M_common.h"
#include "BDT_optimal_cut.h"
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

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/CrossValidation.h"

using namespace std;

void RunT3Mu::Loop(int isMC, TString category){
    
//   In a ROOT session, you can do:
//      root> .L RunT3Mu.C
//      root> RunT3Mu t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   TString fout_path = "datacardT3Mu_"+TMVA_inputpath+category+".root";
   TFile *fout = new TFile(fout_path,"update");
   fout->cd();
   TString datasetName;
   if(!isMC)   datasetName = "data_obs";
   if(isMC==1) datasetName = "signalDs";
   if(isMC==2) datasetName = "signalB0";
   if(isMC==3) datasetName = "signalBp";
   TH1F * hTriplMass1 = new TH1F (datasetName+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
   TH1F * hTriplMass2 = new TH1F (datasetName+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //TMVA Reader initialization
   TMVA::Tools::Instance();
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent:!V" );
   float  forBDTevaluation1; 
   float  forBDTevaluation2; 
   float  forBDTevaluation3; 
   float  forBDTevaluation4; 
   float  forBDTevaluation5; 
   float  forBDTevaluation6; 
   float  forBDTevaluation7; 
   float  forBDTevaluation8; 
   float  forBDTevaluation9; 
   float  forBDTevaluation10;
   float  forBDTevaluation11;
   float  forBDTevaluation12;
   float  forBDTevaluation13;
   Int_t  forBDTevaluation14;

    reader->TMVA::Reader::AddVariable(var_Pmu3,       &forBDTevaluation1 ) ;
    reader->TMVA::Reader::AddVariable(var_cLP,        &forBDTevaluation2 ) ;
    reader->TMVA::Reader::AddVariable(var_tKink,      &forBDTevaluation3 ) ;
    reader->TMVA::Reader::AddVariable(var_segmComp,   &forBDTevaluation4 ) ;
    reader->TMVA::Reader::AddVariable(var_fv_nC,      &forBDTevaluation5 ) ;
    reader->TMVA::Reader::AddVariable(var_fv_dphi3D,  &forBDTevaluation6 ) ;
    reader->TMVA::Reader::AddVariable(var_fv_d3Dsig,  &forBDTevaluation7 ) ;
    reader->TMVA::Reader::AddVariable(var_d0sig,      &forBDTevaluation8 ) ;
    reader->TMVA::Reader::AddVariable(var_mindca_iso, &forBDTevaluation9 ) ;
    reader->TMVA::Reader::AddVariable(var_trkRel,     &forBDTevaluation10 ) ;
    reader->TMVA::Reader::AddVariable( "nMatchesMu3", &forBDTevaluation11 );
    //Spectator variables
    reader->TMVA::Reader::AddSpectator( "tripletMass", &forBDTevaluation12 );
    reader->TMVA::Reader::AddSpectator( "puFactor", &forBDTevaluation13 );
    reader->TMVA::Reader::AddSpectator( "evt := evt % 8192", &forBDTevaluation14 );


   //Book TMVA method
   TString weight_path = TMVA_inputpath+category+TMVA_weightfilename;
   Bool_t weightfileExists = (gSystem->AccessPathName(weight_path) == kFALSE);
   if (weightfileExists) {
      reader->TMVA::Reader::BookMVA(method, weight_path);
      cout<<"Using weights in "<<weight_path<<endl;
   } else {
      std::cout << "Weightfile " <<weight_path<<" for method " << method << " not found."
                   " Did you run TMVACrossValidation with a specified"
                   " splitExpr?" << std::endl;
      exit(0);
   }

/*
   //BDT cut
   TString file_name = "TMVA_file_"+cat_weights+".root";
   TFile *f = new TFile(file_name,"READ");
   TH1F *h_test_signal;
   TH1F *h_test_bkg;
   h_test_signal = (TH1F*)f->Get("dataset_"+cat_weights+"/Method_"+method+"/"+method+"/MVA_"+method+"_S_high");
   h_test_bkg = (TH1F*)f->Get("dataset_"+cat_weights+"/Method_"+method+"/"+method+"/MVA_"+method+"_B_high");
*/

   TString file_name = TMVA_inputpath+category+"/BDTdecision_"+category+".root";
   TFile *f = new TFile(file_name,"READ");
   TH1F *h_test_signal;
   TH1F *h_test_bkg;
   h_test_signal = (TH1F*)f->Get("BDTdecision_signal"+category);
   h_test_bkg = (TH1F*)f->Get("BDTdecision_data_obs"+category);

   BDTcut BDTcutvalues = Get_BDT_cut(category, h_test_signal, h_test_bkg, false);
   cout<<"BDT cut set based on S and B distribution in "<<file_name<<endl;
   cout<<"Category"<< category<<" a cut "<<BDTcutvalues.a<<" b cut "<<BDTcutvalues.b<<endl;

   //Start event Loop
   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      forBDTevaluation1 = Pmu3;
      forBDTevaluation2 = cLP;
      forBDTevaluation3 = tKink;
      forBDTevaluation4 = segmComp;
      forBDTevaluation5 = fv_nC;
      forBDTevaluation6 = fv_dphi3D;
      forBDTevaluation7 = fv_d3Dsig;
      forBDTevaluation8 = d0sig;
      forBDTevaluation9 = mindca_iso;
      forBDTevaluation10 = trkRel;
      forBDTevaluation11 = nMatchesMu3;
      forBDTevaluation12 = tripletMass;
      forBDTevaluation13 = puFactor;
      forBDTevaluation14 = evt % 8192;

      double BDT_decision = reader->EvaluateMVA(method);
   //   cout<<" dec "<< BDT_decision<<" a cut "<<BDTcutvalues.a<<" b cut "<<BDTcutvalues.b<<endl;
    
      if( BDT_decision >= BDTcutvalues.a ) //category 1
           hTriplMass1->Fill(tripletMass);
      else if ( BDT_decision < BDTcutvalues.a && BDT_decision >= BDTcutvalues.b) //category 2
           hTriplMass2->Fill(tripletMass);
   }

   if(isMC){
        //Normalizing Monte Carlo 
        Double_t wNormMC = 1;
        if(isMC == 1) wNormMC = wNormDs * Ds_correction*Dplus_correction; //Ds
        if(isMC == 2) wNormMC = wNormB0 * Ds_correction*Bs_correction*f_correction; //B0
        if(isMC == 3) wNormMC = wNormBp * Ds_correction*Bs_correction*f_correction; //Bp
        cout<<"Normalization factor "<<wNormMC<<endl;
        hTriplMass1->Scale(wNormMC);
        hTriplMass2->Scale(wNormMC);
   } 
   //Write and close the file
   fout->cd();
   hTriplMass1->Write();
   hTriplMass2->Write();
   //duplicate
   if(!isMC) {
       hTriplMass1->Write("background"+category+"1");
       hTriplMass2->Write("background"+category+"2");
   }
   fout->Close();
   delete reader;
   cout<<"Done"<<endl;
}
