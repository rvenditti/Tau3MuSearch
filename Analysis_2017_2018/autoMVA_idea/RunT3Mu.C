#define RunT3Mu_cxx
#include "RunT3Mu.h"
//#include "T3M_common.h"
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
    vector<Float_t> forBDTevaluation;
    for(int j=0; j<VarName.size()-1; j++) forBDTevaluation.push_back(0);
    // For spectator variables
    float  forBDTevaluation1;
    float  forBDTevaluation2;
    Int_t  forBDTevaluation3;
    Int_t  forBDTevaluation4;
    Int_t  forBDTevaluation5;
    
    for(int j=0; j<VarName.size()-1; j++) reader->Reader::AddVariable( VarSel.at(j), &forBDTevaluation.at(j) ) ;

    //Spectator variables
    reader->TMVA::Reader::AddSpectator( "tripletMass", &forBDTevaluation1 );
    reader->TMVA::Reader::AddSpectator( "puFactor", &forBDTevaluation2 );
    //    reader->TMVA::Reader::AddSpectator( "evt := evt % 8192", &forBDTevaluation3 );
    reader->Reader::AddSpectator( "isGlb2", &forBDTevaluation4 );
    reader->Reader::AddSpectator( "isGlb3", &forBDTevaluation5 );


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
       cout <<endl<< "+++++++++++++++++++"<<endl;
       for(int j=0; j<VarName.size()-1; j++) { forBDTevaluation.at(j) = leaf.at(j); cout << "eval_var n. " << j+1 << " | content: " << forBDTevaluation.at(j) << endl;}
       
       forBDTevaluation1 = tripletMass;
       forBDTevaluation2 = puFactor;
       //      forBDTevaluation3 = evt % 8192;
       forBDTevaluation4 = isGlb2;
       forBDTevaluation5 = isGlb3;
       
       if(!(isGlb2 == 1 && isGlb3 == 1 )) continue;

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
