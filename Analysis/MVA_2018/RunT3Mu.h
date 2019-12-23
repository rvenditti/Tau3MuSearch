//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 24 15:04:35 2019 by ROOT version 6.06/09
// from TTree FinalTreeA_Bkg/FinalTreeA_Bkg
// found on file: AnalysedTree_3global_data_2017_tau3mu_03oct.root
//////////////////////////////////////////////////////////

#ifndef RunT3Mu_h
#define RunT3Mu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class RunT3Mu {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Pmu3;
   Double_t        cLP;
   Float_t         tKink;
   Double_t        segmComp;
   Double_t        fv_nC;
   Double_t        fv_dphi3D;
   Double_t        fv_d3Dsig;
   Double_t        d0;
   Double_t        d0sig;
   Double_t        mindca_iso;
   Double_t        trkRel;
   Double_t        Pmu1;
   Double_t        Ptmu1;
   Double_t        Etamu1;
   Double_t        Pmu2;
   Double_t        Ptmu2;
   Double_t        Etamu2;
   Double_t        Ptmu3;
   Double_t        Etamu3;
   Double_t        P_tripl;
   Double_t        Pt_tripl;
   Double_t        Eta_tripl;
   Double_t        nStMu1;
   Double_t        nStMu2;
   Double_t        nStMu3;
   Double_t        Iso03Mu1;
   Double_t        Iso03Mu2;
   Double_t        Iso03Mu3;
   Double_t        Iso05Mu1;
   Double_t        Iso05Mu2;
   Double_t        Iso05Mu3;
   Double_t        nMatchesMu1;
   Double_t        nMatchesMu2;
   Double_t        nMatchesMu3;
   Double_t        timeAtIpInOut1;
   Double_t        timeAtIpInOut2;
   Double_t        timeAtIpInOut3;
   Double_t        cQ_uS;
   Double_t        cQ_tK;
   Double_t        cQ_gK;
   Double_t        cQ_tRChi2;
   Double_t        cQ_sRChi2;
   Double_t        cQ_Chi2LM;
   Double_t        cQ_Chi2lD;
   Double_t        cQ_gDEP;
   Double_t        cQ_tM;
   Double_t        cQ_gTP;
   Double_t        calEn_emMu1;
   Double_t        calEn_emMu2;
   Double_t        calEn_emMu3;
   Double_t        calEn_hadMu1;
   Double_t        calEn_hadMu2;
   Double_t        calEn_hadMu3;
   Double_t        caloComp;
   Double_t        fliDistPVSV_Chi2;
   Double_t        isGlb3;
   Double_t        isTracker3;
   Double_t        isLoose3;
   Double_t        isSoft3;
   Double_t        isPF3;
   Double_t        isRPC3;
   Double_t        isSA3;
   Double_t        isCalo3;
   Double_t        Vx1;
   Double_t        Vx2;
   Double_t        Vx3;
   Double_t        Vy1;
   Double_t        Vy2;
   Double_t        Vy3;
   Double_t        Vz1;
   Double_t        Vz2;
   Double_t        Vz3;
   Double_t        RefVx1;
   Double_t        RefVx2;
   Double_t        RefVx3;
   Double_t        RefVy1;
   Double_t        RefVy2;
   Double_t        RefVy3;
   Double_t        RefVz1;
   Double_t        RefVz2;
   Double_t        RefVz3;
   Double_t        SVx;
   Double_t        SVy;
   Double_t        SVz;
   Double_t        had03;
   Double_t        had05;
   Double_t        nJets03;
   Double_t        nJets05;
   Double_t        nTracks03;
   Double_t        nTracks05;
   Double_t        sumPt03;
   Double_t        sumPt05;
   Double_t        hadVeto03;
   Double_t        hadVeto05;
   Double_t        emVeto03;
   Double_t        emVeto05;
   Double_t        trVeto03;
   Double_t        trVeto05;
   Double_t        tripletMass;
   Double_t        tripletMassReso;

   // List of branches
   TBranch        *b_Pmu3;   //!
   TBranch        *b_cLP;   //!
   TBranch        *b_tKink;   //!
   TBranch        *b_segmComp;   //!
   TBranch        *b_fv_nC;   //!
   TBranch        *b_fv_dphi3D;   //!
   TBranch        *b_fv_d3Dsig;   //!
   TBranch        *b_d0;   //!
   TBranch        *b_d0sig;   //!
   TBranch        *b_mindca_iso;   //!
   TBranch        *b_trkRel;   //!
   TBranch        *b_Pmu1;   //!
   TBranch        *b_Ptmu1;   //!
   TBranch        *b_Etamu1;   //!
   TBranch        *b_Pmu2;   //!
   TBranch        *b_Ptmu2;   //!
   TBranch        *b_Etamu2;   //!
   TBranch        *b_Ptmu3;   //!
   TBranch        *b_Etamu3;   //!
   TBranch        *b_P_tripl;   //!
   TBranch        *b_Pt_tripl;   //!
   TBranch        *b_Eta_tripl;   //!
   TBranch        *b_nStMu1;   //!
   TBranch        *b_nStMu2;   //!
   TBranch        *b_nStMu3;   //!
   TBranch        *b_Iso03Mu1;   //!
   TBranch        *b_Iso03Mu2;   //!
   TBranch        *b_Iso03Mu3;   //!
   TBranch        *b_Iso05Mu1;   //!
   TBranch        *b_Iso05Mu2;   //!
   TBranch        *b_Iso05Mu3;   //!
   TBranch        *b_nMatchesMu1;   //!
   TBranch        *b_nMatchesMu2;   //!
   TBranch        *b_nMatchesMu3;   //!
   TBranch        *b_timeAtIpInOut1;   //!
   TBranch        *b_timeAtIpInOut2;   //!
   TBranch        *b_timeAtIpInOut3;   //!
   TBranch        *b_cQ_uS;   //!
   TBranch        *b_cQ_tK;   //!
   TBranch        *b_cQ_gK;   //!
   TBranch        *b_cQ_tRChi2;   //!
   TBranch        *b_cQ_sRChi2;   //!
   TBranch        *b_cQ_Chi2LM;   //!
   TBranch        *b_cQ_Chi2lD;   //!
   TBranch        *b_cQ_gDEP;   //!
   TBranch        *b_cQ_tM;   //!
   TBranch        *b_cQ_gTP;   //!
   TBranch        *b_calEn_emMu1;   //!
   TBranch        *b_calEn_emMu2;   //!
   TBranch        *b_calEn_emMu3;   //!
   TBranch        *b_calEn_hadMu1;   //!
   TBranch        *b_calEn_hadMu2;   //!
   TBranch        *b_calEn_hadMu3;   //!
   TBranch        *b_caloComp;   //!
   TBranch        *b_fliDistPVSV_Chi2;   //!
   TBranch        *b_isGlb3;   //!
   TBranch        *b_isTracker3;   //!
   TBranch        *b_isLoose3;   //!
   TBranch        *b_isSoft3;   //!
   TBranch        *b_isPF3;   //!
   TBranch        *b_isRPC3;   //!
   TBranch        *b_isSA3;   //!
   TBranch        *b_isCalo3;   //!
   TBranch        *b_Vx1;   //!
   TBranch        *b_Vx2;   //!
   TBranch        *b_Vx3;   //!
   TBranch        *b_Vy1;   //!
   TBranch        *b_Vy2;   //!
   TBranch        *b_Vy3;   //!
   TBranch        *b_Vz1;   //!
   TBranch        *b_Vz2;   //!
   TBranch        *b_Vz3;   //!
   TBranch        *b_RefVx1;   //!
   TBranch        *b_RefVx2;   //!
   TBranch        *b_RefVx3;   //!
   TBranch        *b_RefVy1;   //!
   TBranch        *b_RefVy2;   //!
   TBranch        *b_RefVy3;   //!
   TBranch        *b_RefVz1;   //!
   TBranch        *b_RefVz2;   //!
   TBranch        *b_RefVz3;   //!
   TBranch        *b_SVx;   //!
   TBranch        *b_SVy;   //!
   TBranch        *b_SVz;   //!
   TBranch        *b_had03;   //!
   TBranch        *b_had05;   //!
   TBranch        *b_nJets03;   //!
   TBranch        *b_nJets05;   //!
   TBranch        *b_nTracks03;   //!
   TBranch        *b_nTracks05;   //!
   TBranch        *b_sumPt03;   //!
   TBranch        *b_sumPt05;   //!
   TBranch        *b_hadVeto03;   //!
   TBranch        *b_hadVeto05;   //!
   TBranch        *b_emVeto03;   //!
   TBranch        *b_emVeto05;   //!
   TBranch        *b_trVeto03;   //!
   TBranch        *b_trVeto05;   //!
   TBranch        *b_tripletMass;   //!
   TBranch        *b_tripletMassReso;   //!

   RunT3Mu(TTree *tree=0);
   virtual ~RunT3Mu();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int isMC, TString category);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef RunT3Mu_cxx
RunT3Mu::RunT3Mu(TTree *tree) : fChain(0) 
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

RunT3Mu::~RunT3Mu()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RunT3Mu::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RunT3Mu::LoadTree(Long64_t entry)
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

void RunT3Mu::Init(TTree *tree)
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

   fChain->SetBranchAddress("Pmu3", &Pmu3, &b_Pmu3);
   fChain->SetBranchAddress("cLP", &cLP, &b_cLP);
   fChain->SetBranchAddress("tKink", &tKink, &b_tKink);
   fChain->SetBranchAddress("segmComp", &segmComp, &b_segmComp);
   fChain->SetBranchAddress("fv_nC", &fv_nC, &b_fv_nC);
   fChain->SetBranchAddress("fv_dphi3D", &fv_dphi3D, &b_fv_dphi3D);
   fChain->SetBranchAddress("fv_d3Dsig", &fv_d3Dsig, &b_fv_d3Dsig);
   fChain->SetBranchAddress("d0", &d0, &b_d0);
   fChain->SetBranchAddress("d0sig", &d0sig, &b_d0sig);
   fChain->SetBranchAddress("mindca_iso", &mindca_iso, &b_mindca_iso);
   fChain->SetBranchAddress("trkRel", &trkRel, &b_trkRel);
   fChain->SetBranchAddress("Pmu1", &Pmu1, &b_Pmu1);
   fChain->SetBranchAddress("Ptmu1", &Ptmu1, &b_Ptmu1);
   fChain->SetBranchAddress("Etamu1", &Etamu1, &b_Etamu1);
   fChain->SetBranchAddress("Pmu2", &Pmu2, &b_Pmu2);
   fChain->SetBranchAddress("Ptmu2", &Ptmu2, &b_Ptmu2);
   fChain->SetBranchAddress("Etamu2", &Etamu2, &b_Etamu2);
   fChain->SetBranchAddress("Ptmu3", &Ptmu3, &b_Ptmu3);
   fChain->SetBranchAddress("Etamu3", &Etamu3, &b_Etamu3);
   fChain->SetBranchAddress("P_tripl", &P_tripl, &b_P_tripl);
   fChain->SetBranchAddress("Pt_tripl", &Pt_tripl, &b_Pt_tripl);
   fChain->SetBranchAddress("Eta_tripl", &Eta_tripl, &b_Eta_tripl);
   fChain->SetBranchAddress("nStMu1", &nStMu1, &b_nStMu1);
   fChain->SetBranchAddress("nStMu2", &nStMu2, &b_nStMu2);
   fChain->SetBranchAddress("nStMu3", &nStMu3, &b_nStMu3);
   fChain->SetBranchAddress("Iso03Mu1", &Iso03Mu1, &b_Iso03Mu1);
   fChain->SetBranchAddress("Iso03Mu2", &Iso03Mu2, &b_Iso03Mu2);
   fChain->SetBranchAddress("Iso03Mu3", &Iso03Mu3, &b_Iso03Mu3);
   fChain->SetBranchAddress("Iso05Mu1", &Iso05Mu1, &b_Iso05Mu1);
   fChain->SetBranchAddress("Iso05Mu2", &Iso05Mu2, &b_Iso05Mu2);
   fChain->SetBranchAddress("Iso05Mu3", &Iso05Mu3, &b_Iso05Mu3);
   fChain->SetBranchAddress("nMatchesMu1", &nMatchesMu1, &b_nMatchesMu1);
   fChain->SetBranchAddress("nMatchesMu2", &nMatchesMu2, &b_nMatchesMu2);
   fChain->SetBranchAddress("nMatchesMu3", &nMatchesMu3, &b_nMatchesMu3);
   fChain->SetBranchAddress("timeAtIpInOut1", &timeAtIpInOut1, &b_timeAtIpInOut1);
   fChain->SetBranchAddress("timeAtIpInOut2", &timeAtIpInOut2, &b_timeAtIpInOut2);
   fChain->SetBranchAddress("timeAtIpInOut3", &timeAtIpInOut3, &b_timeAtIpInOut3);
   fChain->SetBranchAddress("cQ_uS", &cQ_uS, &b_cQ_uS);
   fChain->SetBranchAddress("cQ_tK", &cQ_tK, &b_cQ_tK);
   fChain->SetBranchAddress("cQ_gK", &cQ_gK, &b_cQ_gK);
   fChain->SetBranchAddress("cQ_tRChi2", &cQ_tRChi2, &b_cQ_tRChi2);
   fChain->SetBranchAddress("cQ_sRChi2", &cQ_sRChi2, &b_cQ_sRChi2);
   fChain->SetBranchAddress("cQ_Chi2LM", &cQ_Chi2LM, &b_cQ_Chi2LM);
   fChain->SetBranchAddress("cQ_Chi2lD", &cQ_Chi2lD, &b_cQ_Chi2lD);
   fChain->SetBranchAddress("cQ_gDEP", &cQ_gDEP, &b_cQ_gDEP);
   fChain->SetBranchAddress("cQ_tM", &cQ_tM, &b_cQ_tM);
   fChain->SetBranchAddress("cQ_gTP", &cQ_gTP, &b_cQ_gTP);
   fChain->SetBranchAddress("calEn_emMu1", &calEn_emMu1, &b_calEn_emMu1);
   fChain->SetBranchAddress("calEn_emMu2", &calEn_emMu2, &b_calEn_emMu2);
   fChain->SetBranchAddress("calEn_emMu3", &calEn_emMu3, &b_calEn_emMu3);
   fChain->SetBranchAddress("calEn_hadMu1", &calEn_hadMu1, &b_calEn_hadMu1);
   fChain->SetBranchAddress("calEn_hadMu2", &calEn_hadMu2, &b_calEn_hadMu2);
   fChain->SetBranchAddress("calEn_hadMu3", &calEn_hadMu3, &b_calEn_hadMu3);
   fChain->SetBranchAddress("caloComp", &caloComp, &b_caloComp);
   fChain->SetBranchAddress("fliDistPVSV_Chi2", &fliDistPVSV_Chi2, &b_fliDistPVSV_Chi2);
   fChain->SetBranchAddress("isGlb3", &isGlb3, &b_isGlb3);
   fChain->SetBranchAddress("isTracker3", &isTracker3, &b_isTracker3);
   fChain->SetBranchAddress("isLoose3", &isLoose3, &b_isLoose3);
   fChain->SetBranchAddress("isSoft3", &isSoft3, &b_isSoft3);
   fChain->SetBranchAddress("isPF3", &isPF3, &b_isPF3);
   fChain->SetBranchAddress("isRPC3", &isRPC3, &b_isRPC3);
   fChain->SetBranchAddress("isSA3", &isSA3, &b_isSA3);
   fChain->SetBranchAddress("isCalo3", &isCalo3, &b_isCalo3);
   fChain->SetBranchAddress("Vx1", &Vx1, &b_Vx1);
   fChain->SetBranchAddress("Vx2", &Vx2, &b_Vx2);
   fChain->SetBranchAddress("Vx3", &Vx3, &b_Vx3);
   fChain->SetBranchAddress("Vy1", &Vy1, &b_Vy1);
   fChain->SetBranchAddress("Vy2", &Vy2, &b_Vy2);
   fChain->SetBranchAddress("Vy3", &Vy3, &b_Vy3);
   fChain->SetBranchAddress("Vz1", &Vz1, &b_Vz1);
   fChain->SetBranchAddress("Vz2", &Vz2, &b_Vz2);
   fChain->SetBranchAddress("Vz3", &Vz3, &b_Vz3);
   fChain->SetBranchAddress("RefVx1", &RefVx1, &b_RefVx1);
   fChain->SetBranchAddress("RefVx2", &RefVx2, &b_RefVx2);
   fChain->SetBranchAddress("RefVx3", &RefVx3, &b_RefVx3);
   fChain->SetBranchAddress("RefVy1", &RefVy1, &b_RefVy1);
   fChain->SetBranchAddress("RefVy2", &RefVy2, &b_RefVy2);
   fChain->SetBranchAddress("RefVy3", &RefVy3, &b_RefVy3);
   fChain->SetBranchAddress("RefVz1", &RefVz1, &b_RefVz1);
   fChain->SetBranchAddress("RefVz2", &RefVz2, &b_RefVz2);
   fChain->SetBranchAddress("RefVz3", &RefVz3, &b_RefVz3);
   fChain->SetBranchAddress("SVx", &SVx, &b_SVx);
   fChain->SetBranchAddress("SVy", &SVy, &b_SVy);
   fChain->SetBranchAddress("SVz", &SVz, &b_SVz);
   fChain->SetBranchAddress("had03", &had03, &b_had03);
   fChain->SetBranchAddress("had05", &had05, &b_had05);
   fChain->SetBranchAddress("nJets03", &nJets03, &b_nJets03);
   fChain->SetBranchAddress("nJets05", &nJets05, &b_nJets05);
   fChain->SetBranchAddress("nTracks03", &nTracks03, &b_nTracks03);
   fChain->SetBranchAddress("nTracks05", &nTracks05, &b_nTracks05);
   fChain->SetBranchAddress("sumPt03", &sumPt03, &b_sumPt03);
   fChain->SetBranchAddress("sumPt05", &sumPt05, &b_sumPt05);
   fChain->SetBranchAddress("hadVeto03", &hadVeto03, &b_hadVeto03);
   fChain->SetBranchAddress("hadVeto05", &hadVeto05, &b_hadVeto05);
   fChain->SetBranchAddress("emVeto03", &emVeto03, &b_emVeto03);
   fChain->SetBranchAddress("emVeto05", &emVeto05, &b_emVeto05);
   fChain->SetBranchAddress("trVeto03", &trVeto03, &b_trVeto03);
   fChain->SetBranchAddress("trVeto05", &trVeto05, &b_trVeto05);
   fChain->SetBranchAddress("tripletMass", &tripletMass, &b_tripletMass);
   Notify();
}

Bool_t RunT3Mu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RunT3Mu::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RunT3Mu::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RunT3Mu_cxx
