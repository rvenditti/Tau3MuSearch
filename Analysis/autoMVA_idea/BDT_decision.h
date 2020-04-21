//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 24 15:04:35 2019 by ROOT version 6.06/09
// from TTree FinalTreeA_Bkg/FinalTreeA_Bkg
// found on file: AnalysedTree_3global_data_2017_tau3mu_03oct.root
//////////////////////////////////////////////////////////

#ifndef BDT_decision_h
#define BDT_decision_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "T3M_common.h"
using namespace std;

// Header file for the classes stored in the TTree if any.

class BDT_decision {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          evt;
   Double_t        isGlb2;
   Double_t        isGlb3;
   Double_t        tripletMass;
   Double_t        puFactor;
   vector<Double_t> leaf;
    
   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_isGlb2;   //!
   TBranch        *b_isGlb3;   //!
   TBranch        *b_tripletMass;   //!
   TBranch        *b_puFactor;   //!

   BDT_decision(TTree *tree=0);
   virtual ~BDT_decision();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int isMC, TString category);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef BDT_decision_cxx
BDT_decision::BDT_decision(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   //if (tree == 0) {
   //   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("AnalysedTree_data_3globalpt2_2017F_tau3mu_24oct.root");
   //   if (!f || !f->IsOpen()) {
   //      f = new TFile("AnalysedTree_data_3globalpt2_2017F_tau3mu_24oct.root");
   //   }
   //   f->GetObject("FinalTree"+categ+"_Bkg",tree);

   //}
   Init(tree);
}

BDT_decision::~BDT_decision()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BDT_decision::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BDT_decision::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BDT_decision::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("isGlb2", &isGlb2, &b_isGlb2);
   fChain->SetBranchAddress("isGlb3", &isGlb3, &b_isGlb3);
   fChain->SetBranchAddress("tripletMass", &tripletMass, &b_tripletMass);
   fChain->SetBranchAddress("puFactor", &puFactor, &b_puFactor);

    for(int j=0; j<VarName.size()-1; j++) leaf.push_back(0);
    for(int j=0; j<VarName.size()-1; j++) fChain->SetBranchAddress( VarName.at(j), &leaf.at(j) );

   Notify();
}

Bool_t BDT_decision::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BDT_decision::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BDT_decision::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BDT_decision_cxx
