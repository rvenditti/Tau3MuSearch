#define RunT3Mu_cxx
#include "RunT3Mu.h"
#include "Get_BDT_cut.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm>

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


   TFile *fout = new TFile("datacardT3Mu_"+category+".root","update");
   fout->cd();
   TString datasetName;
   if(!isMC) datasetName = "data_";
   if(isMC==1) datasetName = "signalDs_";
   if(isMC==2) datasetName = "signalB0_";
   if(isMC==3) datasetName = "signalBp_";
   TH1F * hTriplMass1 = new TH1F ("hTripl"+datasetName+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
   TH1F * hTriplMass2 = new TH1F ("hTripl"+datasetName+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);

   //BDT cut
   BDTcut BDTcutvalues = Get_BDT_cut(category);
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      float  forBDTevaluation1 = Pmu3;
      float  forBDTevaluation2 = cLP;
      float  forBDTevaluation3 = tKink;
      float  forBDTevaluation4 = segmComp;
      float  forBDTevaluation5 = fv_nC;
      float  forBDTevaluation6 = fv_dphi3D;
      float  forBDTevaluation7 = fv_d3Dsig;
      float  forBDTevaluation8 = d0sig;
      float  forBDTevaluation9 = mindca_iso;
      float  forBDTevaluation10 = tripletMass;

      TMVA::Tools::Instance();
      TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

      reader->TMVA::Reader::AddVariable( "Pmu3", &forBDTevaluation1 );
      reader->TMVA::Reader::AddVariable( "cLP", &forBDTevaluation2 );
      reader->TMVA::Reader::AddVariable( "tKink", &forBDTevaluation3 );
      reader->TMVA::Reader::AddVariable( "segmComp", &forBDTevaluation4 );
      reader->TMVA::Reader::AddVariable( "fv_nC", &forBDTevaluation5 );
      reader->TMVA::Reader::AddVariable( "fv_dphi3D", &forBDTevaluation6 );
      reader->TMVA::Reader::AddVariable( "fv_d3Dsig", &forBDTevaluation7 );
      reader->TMVA::Reader::AddVariable( "d0sig", &forBDTevaluation8 );
      reader->TMVA::Reader::AddVariable( "mindca_iso", &forBDTevaluation9 );
      reader->TMVA::Reader::AddSpectator( "tripletMass", &forBDTevaluation10 );

      reader->TMVA::Reader::BookMVA("BDT", "/lustrehome/fsimone/MVA_Cate/dataset_"+category+"/weights/TMVA_new_BDT.weights.xml");

      double BDT_decision = reader->EvaluateMVA("BDT");
     //EvaluateMVA( Pmu3, cLP,tKink, segmComp,fv_nC,fv_dphi3D, fv_d3Dsig, d0sig, mindca_iso);
       cout<<" dec "<< BDT_decision<<" a cut "<<BDTcutvalues.a<<" b cut "<<BDTcutvalues.b<<endl;
    
       if( BDT_decision >= BDTcutvalues.a ) //category 1
           hTriplMass1->Fill(tripletMass);
       else if ( BDT_decision < BDTcutvalues.a && BDT_decision >= BDTcutvalues.b) //category 2
           hTriplMass2->Fill(tripletMass);
   }

   if(isMC){
        //Normalizing Monte Carlo 
        double DsYieldCorrection = 0.74;
        Double_t wNormMC = 1;
        if(isMC == 1) wNormMC = 0.0034*DsYieldCorrection; //Ds
        if(isMC == 2) wNormMC = 0.0034*DsYieldCorrection; //B0
        if(isMC == 3) wNormMC = 0.0034*DsYieldCorrection; //Bp

        hTriplMass1->Scale(wNormMC);  
        hTriplMass2->Scale(wNormMC);  
    } 
    //Write and close the file
    fout->Write();
    fout->Close();
}
