#define BDT_decision_cxx
#include "BDT_decision.h"
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

void BDT_decision::Loop(int isMC, TString category){

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

      //Backgroud
      if( (isMC==0) && ((tripletMass >= 1.62 && tripletMass <= 1.75) || (tripletMass >= 1.80 && tripletMass <= 2.0 )) && isGlb2 == 1 && isGlb3 == 1 ){
          hBDTdecision->Fill(reader->EvaluateMVA(method));
          cout<<"data in SB - BDT = "<<reader->EvaluateMVA("BDT")<<endl;
      }
      //Signal
      if( isMC>0  && (tripletMass > 1.62 && tripletMass < 2.0) && isGlb2 == 1 && isGlb3 == 1 ){
          hBDTdecision->Fill(reader->EvaluateMVA(method), puFactor);
          cout<<"signal - BDT = "<<reader->EvaluateMVA("BDT")<<endl;
      }

   }//end loop on events

   if(isMC){
        //Normalizing Monte Carlo 
        Double_t wNormMC = 1;
        if(isMC == 1) wNormMC = wNormDs * Ds_correction*Dplus_correction; //Ds
        if(isMC == 2) wNormMC = wNormB0 * Ds_correction*Bs_correction*f_correction; //B0
        if(isMC == 3) wNormMC = wNormBp * Ds_correction*Bs_correction*f_correction; //Bp
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
