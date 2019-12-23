#define BDT_decision_cxx
#include "BDT_decision.h"
#include "T3M_common.h"
#include "BDT_optimal_cut.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/CrossValidation.h"
#include "TMVA/MethodCrossValidation.h"
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"

using namespace std;

void BDT_decision::Loop(int isMC, TString category){
    
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


   TString fout_path = TMVA_inputpath+category+"/BDTdecision_"+category+".root";
   TFile *fout = new TFile(fout_path, "update");
   fout->cd();
   TString datasetName;
   if(!isMC)   datasetName = "data_obs";
   if(isMC==1) datasetName = "signalDs";
   if(isMC==2) datasetName = "signalB0";
   if(isMC==3) datasetName = "signalBp";
   TH1F * hBDTdecision = new TH1F ("BDTdecision_"+datasetName+category,"BDT decision "+category, 240, -0.6, 0.6);

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //TMVA Reader initialization
   TMVA::Tools::Instance();
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
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

   reader->TMVA::Reader::AddVariable( "Pmu3", &forBDTevaluation1 );
   reader->TMVA::Reader::AddVariable( "cLP", &forBDTevaluation2 );
   reader->TMVA::Reader::AddVariable( "tKink", &forBDTevaluation3 );
   reader->TMVA::Reader::AddVariable( "segmComp", &forBDTevaluation4 );
   reader->TMVA::Reader::AddVariable( "fv_nC", &forBDTevaluation5 );
   reader->TMVA::Reader::AddVariable( "fv_dphi3D", &forBDTevaluation6 );
   reader->TMVA::Reader::AddVariable( "fv_d3Dsig", &forBDTevaluation7 );
   reader->TMVA::Reader::AddVariable( "d0sig", &forBDTevaluation8 );
   reader->TMVA::Reader::AddVariable( "mindca_iso", &forBDTevaluation9 );
   reader->TMVA::Reader::AddVariable( "trkRel", &forBDTevaluation10 );
   reader->TMVA::Reader::AddVariable( "nMatchesMu3", &forBDTevaluation11 );
   reader->TMVA::Reader::AddSpectator( "tripletMass", &forBDTevaluation12 );


   //BDT weights
   TString weight_path = TMVA_inputpath+"ABC"+TMVA_weightfilename;
   reader->TMVA::Reader::BookMVA("BDT", weight_path);
   cout<<"Using weights in "<<weight_path<<endl;

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

      //Backgroud
      if( !isMC && (tripletMass >= 1.65 && tripletMass <= 1.73 || tripletMass >= 1.82 && tripletMass <= 1.90) ){
          hBDTdecision->Fill(reader->EvaluateMVA("BDT"));
          //cout<<"data in SB - BDT = "<<reader->EvaluateMVA("BDT")<<endl;
      }
      //Signal
      if( isMC>0 && (tripletMass > 1.73 && tripletMass < 1.82) ){
          hBDTdecision->Fill(reader->EvaluateMVA("BDT"));
          //cout<<"signal - BDT = "<<reader->EvaluateMVA("BDT")<<endl;
      }

   }//end loop on events

   if(isMC){
        //Normalizing Monte Carlo 
        Double_t wNormMC = 1;
        if(isMC == 1) wNormMC = wNormDs * Ds_correction*Dplus_correction; //Ds
        if(isMC == 2) wNormMC = wNormB0 * Ds_correction*Dplus_correction; //B0
        if(isMC == 3) wNormMC = wNormBp * Ds_correction*Dplus_correction; //Bp
        cout<<"Normalization factor "<<wNormMC<<endl;
        hBDTdecision->Scale(wNormMC);
    } 
    //Write and close the file
    fout->cd();
    hBDTdecision->Write();
    fout->Close();
    delete reader;
    cout<<"Done "<<datasetName<<" cat. "<<category<<endl;
}
