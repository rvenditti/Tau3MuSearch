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

   TString cat_weights = category;
   //if(category == "B" || category == "C") cat_weights = "B+C";
   reader->TMVA::Reader::BookMVA("BDT", "/lustrehome/fsimone/MVA_Cate/dataset_"+cat_weights+"/weights/TMVA_new_BDT.weights.xml");

   //BDT cut
   BDTcut BDTcutvalues = Get_BDT_cut(cat_weights);
   cout<<"Category"<< category<<" a cut "<<BDTcutvalues.a<<" b cut "<<BDTcutvalues.b<<endl;
//   if(category=="C") { BDTcutvalues = {0.182048, 0.1}; cout<<"Category C BDT cut reset as:\n"<<BDTcutvalues.a<<" "<<BDTcutvalues.b<<endl; }

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

      double BDT_decision = reader->EvaluateMVA("BDT");
   //   cout<<" dec "<< BDT_decision<<" a cut "<<BDTcutvalues.a<<" b cut "<<BDTcutvalues.b<<endl;
    
      if( BDT_decision >= BDTcutvalues.a ) //category 1
           hTriplMass1->Fill(tripletMass);
      else if ( BDT_decision < BDTcutvalues.a && BDT_decision >= BDTcutvalues.b) //category 2
           hTriplMass2->Fill(tripletMass);
   }

   if(isMC){
        //Normalizing Monte Carlo 
        double DsYieldCorrection = 0.74;
        Double_t wNormMC = 1;
        if(isMC == 1) wNormMC = 0.00115*DsYieldCorrection; //Ds
        if(isMC == 2) wNormMC = 0.00061*DsYieldCorrection; //B0
        if(isMC == 3) wNormMC = 0.00042*DsYieldCorrection; //Bp

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
