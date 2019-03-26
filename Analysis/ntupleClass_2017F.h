//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 25 17:19:53 2019 by ROOT version 6.12/06
// from TTree ntuple/LFVTau ntuple
// found on file: TreeData_Run2017F.root
//////////////////////////////////////////////////////////

#ifndef ntupleClass_2017F_h
#define ntupleClass_2017F_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class ntupleClass_2017F {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          evt;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          nPileUpInt;
   vector<int>     *GenParticle_PdgId;
   vector<double>  *GenParticle_Pt;
   vector<double>  *GenParticle_Eta;
   vector<double>  *GenParticle_Phi;
   vector<int>     *GenParticle_MotherPdgId;
   Int_t           MuonCollectionSize;
   vector<float>   *MuonPt;
   vector<double>  *MuonEnergy;
   vector<double>  *MuonCharge;
   vector<float>   *MuonEta;
   vector<float>   *MuonPhi;
   vector<int>     *Muon_PdgId;
   vector<int>     *Muon_MotherPdgId;
   vector<int>     *Muon_simFlavour;
   vector<float>   *MuonChi2P;
   vector<float>   *MuonChi2LocalPosition;
   vector<float>   *MuonGlbTrackProbability;
   vector<float>   *MuonTrkRelChi2;
   vector<float>   *MuonTrkKink;
   vector<double>  *Muon_vx;
   vector<double>  *Muon_vy;
   vector<double>  *Muon_vz;
   vector<double>  *Muon_isGlobal;
   vector<double>  *Muon_isSoft;
   vector<double>  *Muon_isLoose;
   vector<double>  *Muon_isPF;
   vector<double>  *Muon_isRPCMuon;
   vector<double>  *Muon_isStandAloneMuon;
   vector<double>  *Muon_isTrackerMuon;
   vector<double>  *Muon_isCaloMuon;
   vector<double>  *Muon_isQualityValid;
   vector<double>  *Muon_isTimeValid;
   vector<double>  *Muon_isIsolationValid;
   vector<double>  *Muon_numberOfMatchedStations;
   vector<double>  *Muon_numberOfMatches;
   vector<double>  *Muon_timeAtIpInOut;
   vector<double>  *Muon_GLnormChi2;
   vector<double>  *Muon_GLhitPattern_numberOfValidMuonHits;
   vector<double>  *Muon_trackerLayersWithMeasurement;
   vector<double>  *Muon_Numberofvalidpixelhits;
   vector<double>  *Muon_outerTrack_p;
   vector<double>  *Muon_outerTrack_eta;
   vector<double>  *Muon_outerTrack_phi;
   vector<double>  *Muon_outerTrack_normalizedChi2;
   vector<double>  *Muon_outerTrack_muonStationsWithValidHits;
   vector<double>  *Muon_innerTrack_p;
   vector<double>  *Muon_innerTrack_eta;
   vector<double>  *Muon_innerTrack_phi;
   vector<double>  *Muon_innerTrack_normalizedChi2;
   vector<double>  *Muon_QInnerOuter;
   vector<double>  *Muon_combinedQuality_updatedSta;
   vector<double>  *Muon_combinedQuality_trkKink;
   vector<double>  *Muon_combinedQuality_glbKink;
   vector<double>  *Muon_combinedQuality_trkRelChi2;
   vector<double>  *Muon_combinedQuality_staRelChi2;
   vector<double>  *Muon_combinedQuality_chi2LocalPosition;
   vector<double>  *Muon_combinedQuality_chi2LocalMomentum;
   vector<double>  *Muon_combinedQuality_localDistance;
   vector<double>  *Muon_combinedQuality_globalDeltaEtaPhi;
   vector<double>  *Muon_combinedQuality_tightMatch;
   vector<double>  *Muon_combinedQuality_glbTrackProbability;
   vector<double>  *Muon_calEnergy_em;
   vector<double>  *Muon_calEnergy_emS9;
   vector<double>  *Muon_calEnergy_emS25;
   vector<double>  *Muon_calEnergy_had;
   vector<double>  *Muon_calEnergy_hadS9;
   vector<double>  *Muon_segmentCompatibility;
   vector<double>  *Muon_caloCompatibility;
   vector<double>  *Muon_ptErrOverPt;
   vector<double>  *Muon_emEt03;
   vector<double>  *Muon_hadEt03;
   vector<double>  *Muon_nJets03;
   vector<double>  *Muon_nTracks03;
   vector<double>  *Muon_sumPt03;
   vector<double>  *Muon_hadVetoEt03;
   vector<double>  *Muon_emVetoEt03;
   vector<double>  *Muon_trackerVetoPt03;
   vector<double>  *Muon_emEt05;
   vector<double>  *Muon_hadEt05;
   vector<double>  *Muon_nJets05;
   vector<double>  *Muon_nTracks05;
   vector<double>  *Muon_sumPt05;
   vector<double>  *Muon_hadVetoEt05;
   vector<double>  *Muon_emVetoEt05;
   vector<double>  *Muon_trackerVetoPt05;
   Int_t           PVCollection_Size;
   Double_t        PV_x;
   Double_t        PV_y;
   Double_t        PV_z;
   Double_t        PV_NTracks;
   Int_t           TripletCollectionSize;
   vector<double>  *Mu1_Pt;
   vector<double>  *Mu1_Eta;
   vector<double>  *Mu1_Phi;
   vector<int>     *Mu1_TripletIndex;
   vector<double>  *Mu2_Pt;
   vector<double>  *Mu2_Eta;
   vector<double>  *Mu2_Phi;
   vector<int>     *Mu2_TripletIndex;
   vector<double>  *Mu3_Pt;
   vector<double>  *Mu3_Eta;
   vector<double>  *Mu3_Phi;
   vector<int>     *Mu3_TripletIndex;
   vector<double>  *GenMatchMu1_SimPt;
   vector<double>  *GenMatchMu2_SimPt;
   vector<double>  *GenMatchMu3_SimPt;
   vector<double>  *GenMatchMu1_SimEta;
   vector<double>  *GenMatchMu2_SimEta;
   vector<double>  *GenMatchMu3_SimEta;
   vector<double>  *GenMatchMu1_SimPhi;
   vector<double>  *GenMatchMu2_SimPhi;
   vector<double>  *GenMatchMu3_SimPhi;
   vector<double>  *GenMatchMu1_Pt;
   vector<double>  *GenMatchMu2_Pt;
   vector<double>  *GenMatchMu3_Pt;
   vector<double>  *GenMatchMu1_Eta;
   vector<double>  *GenMatchMu2_Eta;
   vector<double>  *GenMatchMu3_Eta;
   vector<double>  *GenMatchMu1_Phi;
   vector<double>  *GenMatchMu2_Phi;
   vector<double>  *GenMatchMu3_Phi;
   vector<double>  *TripletVtx_x;
   vector<double>  *TripletVtx_y;
   vector<double>  *TripletVtx_z;
   vector<double>  *TripletVtx_Chi2;
   vector<double>  *TripletVtx_NDOF;
   vector<double>  *Triplet_Mass;
   vector<double>  *Triplet_Pt;
   vector<double>  *Triplet_Eta;
   vector<double>  *Triplet_Phi;
   vector<double>  *Triplet_Charge;
   vector<double>  *RefittedPV_x;
   vector<double>  *RefittedPV_y;
   vector<double>  *RefittedPV_z;
   vector<double>  *RefittedPV_NTracks;
   vector<int>     *RefittedPV_isValid;
   vector<double>  *FlightDistPVSV;
   vector<double>  *FlightDistPVSV_Err;
   vector<double>  *FlightDistPVSV_Significance;
   vector<double>  *FlightDistPVSV_chi2;

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_nPileUpInt;   //!
   TBranch        *b_GenParticle_PdgId;   //!
   TBranch        *b_GenParticle_Pt;   //!
   TBranch        *b_GenParticle_Eta;   //!
   TBranch        *b_GenParticle_Phi;   //!
   TBranch        *b_GenParticle_MotherPdgId;   //!
   TBranch        *b_MuonCollectionSize;   //!
   TBranch        *b_MuonPt;   //!
   TBranch        *b_MuonEnergy;   //!
   TBranch        *b_MuonCharge;   //!
   TBranch        *b_MuonEta;   //!
   TBranch        *b_MuonPhi;   //!
   TBranch        *b_Muon_PdgId;   //!
   TBranch        *b_Muon_MotherPdgId;   //!
   TBranch        *b_Muon_simFlavour;   //!
   TBranch        *b_MuonChi2P;   //!
   TBranch        *b_MuonChi2LocalPosition;   //!
   TBranch        *b_MuonGlbTrackProbability;   //!
   TBranch        *b_MuonTrkRelChi2;   //!
   TBranch        *b_MuonTrkKink;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_isGlobal;   //!
   TBranch        *b_Muon_isSoft;   //!
   TBranch        *b_Muon_isLoose;   //!
   TBranch        *b_Muon_isPF;   //!
   TBranch        *b_Muon_isRPCMuon;   //!
   TBranch        *b_Muon_isStandAloneMuon;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isCaloMuon;   //!
   TBranch        *b_Muon_isQualityValid;   //!
   TBranch        *b_Muon_isTimeValid;   //!
   TBranch        *b_Muon_isIsolationValid;   //!
   TBranch        *b_Muon_numberOfMatchedStations;   //!
   TBranch        *b_Muon_numberOfMatches;   //!
   TBranch        *b_Muon_timeAtIpInOut;   //!
   TBranch        *b_Muon_GLnormChi2;   //!
   TBranch        *b_Muon_GLhitPattern_numberOfValidMuonHits;   //!
   TBranch        *b_Muon_trackerLayersWithMeasurement;   //!
   TBranch        *b_Muon_Numberofvalidpixelhits;   //!
   TBranch        *b_Muon_outerTrack_p;   //!
   TBranch        *b_Muon_outerTrack_eta;   //!
   TBranch        *b_Muon_outerTrack_phi;   //!
   TBranch        *b_Muon_outerTrack_normalizedChi2;   //!
   TBranch        *b_Muon_outerTrack_muonStationsWithValidHits;   //!
   TBranch        *b_Muon_innerTrack_p;   //!
   TBranch        *b_Muon_innerTrack_eta;   //!
   TBranch        *b_Muon_innerTrack_phi;   //!
   TBranch        *b_Muon_innerTrack_normalizedChi2;   //!
   TBranch        *b_Muon_QInnerOuter;   //!
   TBranch        *b_Muon_combinedQuality_updatedSta;   //!
   TBranch        *b_Muon_combinedQuality_trkKink;   //!
   TBranch        *b_Muon_combinedQuality_glbKink;   //!
   TBranch        *b_Muon_combinedQuality_trkRelChi2;   //!
   TBranch        *b_Muon_combinedQuality_staRelChi2;   //!
   TBranch        *b_Muon_combinedQuality_chi2LocalPosition;   //!
   TBranch        *b_Muon_combinedQuality_chi2LocalMomentum;   //!
   TBranch        *b_Muon_combinedQuality_localDistance;   //!
   TBranch        *b_Muon_combinedQuality_globalDeltaEtaPhi;   //!
   TBranch        *b_Muon_combinedQuality_tightMatch;   //!
   TBranch        *b_Muon_combinedQuality_glbTrackProbability;   //!
   TBranch        *b_Muon_calEnergy_em;   //!
   TBranch        *b_Muon_calEnergy_emS9;   //!
   TBranch        *b_Muon_calEnergy_emS25;   //!
   TBranch        *b_Muon_calEnergy_had;   //!
   TBranch        *b_Muon_calEnergy_hadS9;   //!
   TBranch        *b_Muon_segmentCompatibility;   //!
   TBranch        *b_Muon_caloCompatibility;   //!
   TBranch        *b_Muon_ptErrOverPt;   //!
   TBranch        *b_Muon_emEt03;   //!
   TBranch        *b_Muon_hadEt03;   //!
   TBranch        *b_Muon_nJets03;   //!
   TBranch        *b_Muon_nTracks03;   //!
   TBranch        *b_Muon_sumPt03;   //!
   TBranch        *b_Muon_hadVetoEt03;   //!
   TBranch        *b_Muon_emVetoEt03;   //!
   TBranch        *b_Muon_trackerVetoPt03;   //!
   TBranch        *b_Muon_emEt05;   //!
   TBranch        *b_Muon_hadEt05;   //!
   TBranch        *b_Muon_nJets05;   //!
   TBranch        *b_Muon_nTracks05;   //!
   TBranch        *b_Muon_sumPt05;   //!
   TBranch        *b_Muon_hadVetoEt05;   //!
   TBranch        *b_Muon_emVetoEt05;   //!
   TBranch        *b_Muon_trackerVetoPt05;   //!
   TBranch        *b_PVCollection_Size;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PV_NTracks;   //!
   TBranch        *b_TripletCollectionSize;   //!
   TBranch        *b_Mu1_Pt;   //!
   TBranch        *b_Mu1_Eta;   //!
   TBranch        *b_Mu1_Phi;   //!
   TBranch        *b_Mu1_TripletIndex;   //!
   TBranch        *b_Mu2_Pt;   //!
   TBranch        *b_Mu2_Eta;   //!
   TBranch        *b_Mu2_Phi;   //!
   TBranch        *b_Mu2_TripletIndex;   //!
   TBranch        *b_Mu3_Pt;   //!
   TBranch        *b_Mu3_Eta;   //!
   TBranch        *b_Mu3_Phi;   //!
   TBranch        *b_Mu3_TripletIndex;   //!
   TBranch        *b_GenMatchMu1_SimPt;   //!
   TBranch        *b_GenMatchMu2_SimPt;   //!
   TBranch        *b_GenMatchMu3_SimPt;   //!
   TBranch        *b_GenMatchMu1_SimEta;   //!
   TBranch        *b_GenMatchMu2_SimEta;   //!
   TBranch        *b_GenMatchMu3_SimEta;   //!
   TBranch        *b_GenMatchMu1_SimPhi;   //!
   TBranch        *b_GenMatchMu2_SimPhi;   //!
   TBranch        *b_GenMatchMu3_SimPhi;   //!
   TBranch        *b_GenMatchMu1_Pt;   //!
   TBranch        *b_GenMatchMu2_Pt;   //!
   TBranch        *b_GenMatchMu3_Pt;   //!
   TBranch        *b_GenMatchMu1_Eta;   //!
   TBranch        *b_GenMatchMu2_Eta;   //!
   TBranch        *b_GenMatchMu3_Eta;   //!
   TBranch        *b_GenMatchMu1_Phi;   //!
   TBranch        *b_GenMatchMu2_Phi;   //!
   TBranch        *b_GenMatchMu3_Phi;   //!
   TBranch        *b_TripletVtx_x;   //!
   TBranch        *b_TripletVtx_y;   //!
   TBranch        *b_TripletVtx_z;   //!
   TBranch        *b_TripletVtx_Chi2;   //!
   TBranch        *b_TripletVtx_NDOF;   //!
   TBranch        *b_Triplet_Mass;   //!
   TBranch        *b_Triplet_Pt;   //!
   TBranch        *b_Triplet_Eta;   //!
   TBranch        *b_Triplet_Phi;   //!
   TBranch        *b_Triplet_Charge;   //!
   TBranch        *b_RefittedPV_x;   //!
   TBranch        *b_RefittedPV_y;   //!
   TBranch        *b_RefittedPV_z;   //!
   TBranch        *b_RefittedPV_NTracks;   //!
   TBranch        *b_RefittedPV_isValid;   //!
   TBranch        *b_FlightDistPVSV;   //!
   TBranch        *b_FlightDistPVSV_Err;   //!
   TBranch        *b_FlightDistPVSV_Significance;   //!
   TBranch        *b_FlightDistPVSV_chi2;   //!

   ntupleClass_2017F(TTree *tree=0);
    virtual ~ntupleClass_2017F();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop_New();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual Double_t MuonFinder(Double_t pt, Double_t eta, Double_t phi);
    virtual Bool_t   isDeltaRGood(Float_t eta1, Float_t eta2, Float_t phi1, Float_t phi2, Float_t DeltaRmax);
    virtual Bool_t   isDeltaZGood(Float_t vz1, Float_t vz2, Float_t DeltaZmax);
    virtual Bool_t   isPairDeltaRGood(Int_t ntriplet, Float_t DeltaRmax);
    virtual Bool_t   isPairDeltaZGood(Float_t DeltaZ1, Float_t DeltaZ2, Float_t DeltaZ3, Float_t DeltaZmax);
    virtual Double_t MuonP(Double_t pt, Double_t eta, Double_t phi, Double_t mass);
    virtual Double_t DimuonMass(Double_t charge1, Double_t charge2, Double_t pt1, Double_t pt2, Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2, Double_t mass);
    virtual Float_t  QuadMuonMass(Float_t pt1, Float_t pt2, Float_t pt3, Float_t pt4, Float_t eta1, Float_t eta2, Float_t eta3, Float_t eta4, Float_t phi1, Float_t phi2, Float_t phi3, Float_t phi4, Double_t mass);
    virtual Bool_t  isNotAPhi(Double_t dimumass, Double_t sigma);
    virtual Bool_t  isPairNotAPhi(Double_t dimu1_2, Double_t dimu2_3, Double_t dimu1_3, Double_t sigma);
    virtual void    particleName(TString pId[37]);
    virtual void    particleId(Int_t pdgId, Int_t Idsummary[37]);
    virtual void    cutName(TString listCut[7]);
    virtual void    cutName_New(TString listCut[7]);
    virtual void    cutName2(TString listCut[7]);
    virtual void    PdgIdHisto(TCanvas *canv, TH1I *hist, Int_t Idsummary[37], TString pIdList[37]);
    virtual void    CutEffHisto(TCanvas *canv, TH1I *hist, Int_t cut[7], TString listCut[7]);
    virtual void    DiMuonHisto(TH1D *hist, Double_t dimu1_2, Double_t dimu2_3, Double_t dimu1_3);
};

#endif

#ifdef ntupleClass_2017F_cxx
ntupleClass_2017F::ntupleClass_2017F(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TreeData_Run2017F.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("TreeData_Run2017F.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("TreeData_Run2017F.root:/TreeMakerBkg");
      dir->GetObject("ntuple",tree);

   }
   Init(tree);
}

ntupleClass_2017F::~ntupleClass_2017F()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntupleClass_2017F::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntupleClass_2017F::LoadTree(Long64_t entry)
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

void ntupleClass_2017F::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   GenParticle_PdgId = 0;
   GenParticle_Pt = 0;
   GenParticle_Eta = 0;
   GenParticle_Phi = 0;
   GenParticle_MotherPdgId = 0;
   MuonPt = 0;
   MuonEnergy = 0;
   MuonCharge = 0;
   MuonEta = 0;
   MuonPhi = 0;
   Muon_PdgId = 0;
   Muon_MotherPdgId = 0;
   Muon_simFlavour = 0;
   MuonChi2P = 0;
   MuonChi2LocalPosition = 0;
   MuonGlbTrackProbability = 0;
   MuonTrkRelChi2 = 0;
   MuonTrkKink = 0;
   Muon_vx = 0;
   Muon_vy = 0;
   Muon_vz = 0;
   Muon_isGlobal = 0;
   Muon_isSoft = 0;
   Muon_isLoose = 0;
   Muon_isPF = 0;
   Muon_isRPCMuon = 0;
   Muon_isStandAloneMuon = 0;
   Muon_isTrackerMuon = 0;
   Muon_isCaloMuon = 0;
   Muon_isQualityValid = 0;
   Muon_isTimeValid = 0;
   Muon_isIsolationValid = 0;
   Muon_numberOfMatchedStations = 0;
   Muon_numberOfMatches = 0;
   Muon_timeAtIpInOut = 0;
   Muon_GLnormChi2 = 0;
   Muon_GLhitPattern_numberOfValidMuonHits = 0;
   Muon_trackerLayersWithMeasurement = 0;
   Muon_Numberofvalidpixelhits = 0;
   Muon_outerTrack_p = 0;
   Muon_outerTrack_eta = 0;
   Muon_outerTrack_phi = 0;
   Muon_outerTrack_normalizedChi2 = 0;
   Muon_outerTrack_muonStationsWithValidHits = 0;
   Muon_innerTrack_p = 0;
   Muon_innerTrack_eta = 0;
   Muon_innerTrack_phi = 0;
   Muon_innerTrack_normalizedChi2 = 0;
   Muon_QInnerOuter = 0;
   Muon_combinedQuality_updatedSta = 0;
   Muon_combinedQuality_trkKink = 0;
   Muon_combinedQuality_glbKink = 0;
   Muon_combinedQuality_trkRelChi2 = 0;
   Muon_combinedQuality_staRelChi2 = 0;
   Muon_combinedQuality_chi2LocalPosition = 0;
   Muon_combinedQuality_chi2LocalMomentum = 0;
   Muon_combinedQuality_localDistance = 0;
   Muon_combinedQuality_globalDeltaEtaPhi = 0;
   Muon_combinedQuality_tightMatch = 0;
   Muon_combinedQuality_glbTrackProbability = 0;
   Muon_calEnergy_em = 0;
   Muon_calEnergy_emS9 = 0;
   Muon_calEnergy_emS25 = 0;
   Muon_calEnergy_had = 0;
   Muon_calEnergy_hadS9 = 0;
   Muon_segmentCompatibility = 0;
   Muon_caloCompatibility = 0;
   Muon_ptErrOverPt = 0;
   Muon_emEt03 = 0;
   Muon_hadEt03 = 0;
   Muon_nJets03 = 0;
   Muon_nTracks03 = 0;
   Muon_sumPt03 = 0;
   Muon_hadVetoEt03 = 0;
   Muon_emVetoEt03 = 0;
   Muon_trackerVetoPt03 = 0;
   Muon_emEt05 = 0;
   Muon_hadEt05 = 0;
   Muon_nJets05 = 0;
   Muon_nTracks05 = 0;
   Muon_sumPt05 = 0;
   Muon_hadVetoEt05 = 0;
   Muon_emVetoEt05 = 0;
   Muon_trackerVetoPt05 = 0;
   Mu1_Pt = 0;
   Mu1_Eta = 0;
   Mu1_Phi = 0;
   Mu1_TripletIndex = 0;
   Mu2_Pt = 0;
   Mu2_Eta = 0;
   Mu2_Phi = 0;
   Mu2_TripletIndex = 0;
   Mu3_Pt = 0;
   Mu3_Eta = 0;
   Mu3_Phi = 0;
   Mu3_TripletIndex = 0;
   GenMatchMu1_SimPt = 0;
   GenMatchMu2_SimPt = 0;
   GenMatchMu3_SimPt = 0;
   GenMatchMu1_SimEta = 0;
   GenMatchMu2_SimEta = 0;
   GenMatchMu3_SimEta = 0;
   GenMatchMu1_SimPhi = 0;
   GenMatchMu2_SimPhi = 0;
   GenMatchMu3_SimPhi = 0;
   GenMatchMu1_Pt = 0;
   GenMatchMu2_Pt = 0;
   GenMatchMu3_Pt = 0;
   GenMatchMu1_Eta = 0;
   GenMatchMu2_Eta = 0;
   GenMatchMu3_Eta = 0;
   GenMatchMu1_Phi = 0;
   GenMatchMu2_Phi = 0;
   GenMatchMu3_Phi = 0;
   TripletVtx_x = 0;
   TripletVtx_y = 0;
   TripletVtx_z = 0;
   TripletVtx_Chi2 = 0;
   TripletVtx_NDOF = 0;
   Triplet_Mass = 0;
   Triplet_Pt = 0;
   Triplet_Eta = 0;
   Triplet_Phi = 0;
   Triplet_Charge = 0;
   RefittedPV_x = 0;
   RefittedPV_y = 0;
   RefittedPV_z = 0;
   RefittedPV_NTracks = 0;
   RefittedPV_isValid = 0;
   FlightDistPVSV = 0;
   FlightDistPVSV_Err = 0;
   FlightDistPVSV_Significance = 0;
   FlightDistPVSV_chi2 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("nPileUpInt", &nPileUpInt, &b_nPileUpInt);
   fChain->SetBranchAddress("GenParticle_PdgId", &GenParticle_PdgId, &b_GenParticle_PdgId);
   fChain->SetBranchAddress("GenParticle_Pt", &GenParticle_Pt, &b_GenParticle_Pt);
   fChain->SetBranchAddress("GenParticle_Eta", &GenParticle_Eta, &b_GenParticle_Eta);
   fChain->SetBranchAddress("GenParticle_Phi", &GenParticle_Phi, &b_GenParticle_Phi);
   fChain->SetBranchAddress("GenParticle_MotherPdgId", &GenParticle_MotherPdgId, &b_GenParticle_MotherPdgId);
   fChain->SetBranchAddress("MuonCollectionSize", &MuonCollectionSize, &b_MuonCollectionSize);
   fChain->SetBranchAddress("MuonPt", &MuonPt, &b_MuonPt);
   fChain->SetBranchAddress("MuonEnergy", &MuonEnergy, &b_MuonEnergy);
   fChain->SetBranchAddress("MuonCharge", &MuonCharge, &b_MuonCharge);
   fChain->SetBranchAddress("MuonEta", &MuonEta, &b_MuonEta);
   fChain->SetBranchAddress("MuonPhi", &MuonPhi, &b_MuonPhi);
   fChain->SetBranchAddress("Muon_PdgId", &Muon_PdgId, &b_Muon_PdgId);
   fChain->SetBranchAddress("Muon_MotherPdgId", &Muon_MotherPdgId, &b_Muon_MotherPdgId);
   fChain->SetBranchAddress("Muon_simFlavour", &Muon_simFlavour, &b_Muon_simFlavour);
   fChain->SetBranchAddress("MuonChi2P", &MuonChi2P, &b_MuonChi2P);
   fChain->SetBranchAddress("MuonChi2LocalPosition", &MuonChi2LocalPosition, &b_MuonChi2LocalPosition);
   fChain->SetBranchAddress("MuonGlbTrackProbability", &MuonGlbTrackProbability, &b_MuonGlbTrackProbability);
   fChain->SetBranchAddress("MuonTrkRelChi2", &MuonTrkRelChi2, &b_MuonTrkRelChi2);
   fChain->SetBranchAddress("MuonTrkKink", &MuonTrkKink, &b_MuonTrkKink);
   fChain->SetBranchAddress("Muon_vx", &Muon_vx, &b_Muon_vx);
   fChain->SetBranchAddress("Muon_vy", &Muon_vy, &b_Muon_vy);
   fChain->SetBranchAddress("Muon_vz", &Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon_isGlobal", &Muon_isGlobal, &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isSoft", &Muon_isSoft, &b_Muon_isSoft);
   fChain->SetBranchAddress("Muon_isLoose", &Muon_isLoose, &b_Muon_isLoose);
   fChain->SetBranchAddress("Muon_isPF", &Muon_isPF, &b_Muon_isPF);
   fChain->SetBranchAddress("Muon_isRPCMuon", &Muon_isRPCMuon, &b_Muon_isRPCMuon);
   fChain->SetBranchAddress("Muon_isStandAloneMuon", &Muon_isStandAloneMuon, &b_Muon_isStandAloneMuon);
   fChain->SetBranchAddress("Muon_isTrackerMuon", &Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon_isCaloMuon", &Muon_isCaloMuon, &b_Muon_isCaloMuon);
   fChain->SetBranchAddress("Muon_isQualityValid", &Muon_isQualityValid, &b_Muon_isQualityValid);
   fChain->SetBranchAddress("Muon_isTimeValid", &Muon_isTimeValid, &b_Muon_isTimeValid);
   fChain->SetBranchAddress("Muon_isIsolationValid", &Muon_isIsolationValid, &b_Muon_isIsolationValid);
   fChain->SetBranchAddress("Muon_numberOfMatchedStations", &Muon_numberOfMatchedStations, &b_Muon_numberOfMatchedStations);
   fChain->SetBranchAddress("Muon_numberOfMatches", &Muon_numberOfMatches, &b_Muon_numberOfMatches);
   fChain->SetBranchAddress("Muon_timeAtIpInOut", &Muon_timeAtIpInOut, &b_Muon_timeAtIpInOut);
   fChain->SetBranchAddress("Muon_GLnormChi2", &Muon_GLnormChi2, &b_Muon_GLnormChi2);
   fChain->SetBranchAddress("Muon_GLhitPattern_numberOfValidMuonHits", &Muon_GLhitPattern_numberOfValidMuonHits, &b_Muon_GLhitPattern_numberOfValidMuonHits);
   fChain->SetBranchAddress("Muon_trackerLayersWithMeasurement", &Muon_trackerLayersWithMeasurement, &b_Muon_trackerLayersWithMeasurement);
   fChain->SetBranchAddress("Muon_Numberofvalidpixelhits", &Muon_Numberofvalidpixelhits, &b_Muon_Numberofvalidpixelhits);
   fChain->SetBranchAddress("Muon_outerTrack_p", &Muon_outerTrack_p, &b_Muon_outerTrack_p);
   fChain->SetBranchAddress("Muon_outerTrack_eta", &Muon_outerTrack_eta, &b_Muon_outerTrack_eta);
   fChain->SetBranchAddress("Muon_outerTrack_phi", &Muon_outerTrack_phi, &b_Muon_outerTrack_phi);
   fChain->SetBranchAddress("Muon_outerTrack_normalizedChi2", &Muon_outerTrack_normalizedChi2, &b_Muon_outerTrack_normalizedChi2);
   fChain->SetBranchAddress("Muon_outerTrack_muonStationsWithValidHits", &Muon_outerTrack_muonStationsWithValidHits, &b_Muon_outerTrack_muonStationsWithValidHits);
   fChain->SetBranchAddress("Muon_innerTrack_p", &Muon_innerTrack_p, &b_Muon_innerTrack_p);
   fChain->SetBranchAddress("Muon_innerTrack_eta", &Muon_innerTrack_eta, &b_Muon_innerTrack_eta);
   fChain->SetBranchAddress("Muon_innerTrack_phi", &Muon_innerTrack_phi, &b_Muon_innerTrack_phi);
   fChain->SetBranchAddress("Muon_innerTrack_normalizedChi2", &Muon_innerTrack_normalizedChi2, &b_Muon_innerTrack_normalizedChi2);
   fChain->SetBranchAddress("Muon_QInnerOuter", &Muon_QInnerOuter, &b_Muon_QInnerOuter);
   fChain->SetBranchAddress("Muon_combinedQuality_updatedSta", &Muon_combinedQuality_updatedSta, &b_Muon_combinedQuality_updatedSta);
   fChain->SetBranchAddress("Muon_combinedQuality_trkKink", &Muon_combinedQuality_trkKink, &b_Muon_combinedQuality_trkKink);
   fChain->SetBranchAddress("Muon_combinedQuality_glbKink", &Muon_combinedQuality_glbKink, &b_Muon_combinedQuality_glbKink);
   fChain->SetBranchAddress("Muon_combinedQuality_trkRelChi2", &Muon_combinedQuality_trkRelChi2, &b_Muon_combinedQuality_trkRelChi2);
   fChain->SetBranchAddress("Muon_combinedQuality_staRelChi2", &Muon_combinedQuality_staRelChi2, &b_Muon_combinedQuality_staRelChi2);
   fChain->SetBranchAddress("Muon_combinedQuality_chi2LocalPosition", &Muon_combinedQuality_chi2LocalPosition, &b_Muon_combinedQuality_chi2LocalPosition);
   fChain->SetBranchAddress("Muon_combinedQuality_chi2LocalMomentum", &Muon_combinedQuality_chi2LocalMomentum, &b_Muon_combinedQuality_chi2LocalMomentum);
   fChain->SetBranchAddress("Muon_combinedQuality_localDistance", &Muon_combinedQuality_localDistance, &b_Muon_combinedQuality_localDistance);
   fChain->SetBranchAddress("Muon_combinedQuality_globalDeltaEtaPhi", &Muon_combinedQuality_globalDeltaEtaPhi, &b_Muon_combinedQuality_globalDeltaEtaPhi);
   fChain->SetBranchAddress("Muon_combinedQuality_tightMatch", &Muon_combinedQuality_tightMatch, &b_Muon_combinedQuality_tightMatch);
   fChain->SetBranchAddress("Muon_combinedQuality_glbTrackProbability", &Muon_combinedQuality_glbTrackProbability, &b_Muon_combinedQuality_glbTrackProbability);
   fChain->SetBranchAddress("Muon_calEnergy_em", &Muon_calEnergy_em, &b_Muon_calEnergy_em);
   fChain->SetBranchAddress("Muon_calEnergy_emS9", &Muon_calEnergy_emS9, &b_Muon_calEnergy_emS9);
   fChain->SetBranchAddress("Muon_calEnergy_emS25", &Muon_calEnergy_emS25, &b_Muon_calEnergy_emS25);
   fChain->SetBranchAddress("Muon_calEnergy_had", &Muon_calEnergy_had, &b_Muon_calEnergy_had);
   fChain->SetBranchAddress("Muon_calEnergy_hadS9", &Muon_calEnergy_hadS9, &b_Muon_calEnergy_hadS9);
   fChain->SetBranchAddress("Muon_segmentCompatibility", &Muon_segmentCompatibility, &b_Muon_segmentCompatibility);
   fChain->SetBranchAddress("Muon_caloCompatibility", &Muon_caloCompatibility, &b_Muon_caloCompatibility);
   fChain->SetBranchAddress("Muon_ptErrOverPt", &Muon_ptErrOverPt, &b_Muon_ptErrOverPt);
   fChain->SetBranchAddress("Muon_emEt03", &Muon_emEt03, &b_Muon_emEt03);
   fChain->SetBranchAddress("Muon_hadEt03", &Muon_hadEt03, &b_Muon_hadEt03);
   fChain->SetBranchAddress("Muon_nJets03", &Muon_nJets03, &b_Muon_nJets03);
   fChain->SetBranchAddress("Muon_nTracks03", &Muon_nTracks03, &b_Muon_nTracks03);
   fChain->SetBranchAddress("Muon_sumPt03", &Muon_sumPt03, &b_Muon_sumPt03);
   fChain->SetBranchAddress("Muon_hadVetoEt03", &Muon_hadVetoEt03, &b_Muon_hadVetoEt03);
   fChain->SetBranchAddress("Muon_emVetoEt03", &Muon_emVetoEt03, &b_Muon_emVetoEt03);
   fChain->SetBranchAddress("Muon_trackerVetoPt03", &Muon_trackerVetoPt03, &b_Muon_trackerVetoPt03);
   fChain->SetBranchAddress("Muon_emEt05", &Muon_emEt05, &b_Muon_emEt05);
   fChain->SetBranchAddress("Muon_hadEt05", &Muon_hadEt05, &b_Muon_hadEt05);
   fChain->SetBranchAddress("Muon_nJets05", &Muon_nJets05, &b_Muon_nJets05);
   fChain->SetBranchAddress("Muon_nTracks05", &Muon_nTracks05, &b_Muon_nTracks05);
   fChain->SetBranchAddress("Muon_sumPt05", &Muon_sumPt05, &b_Muon_sumPt05);
   fChain->SetBranchAddress("Muon_hadVetoEt05", &Muon_hadVetoEt05, &b_Muon_hadVetoEt05);
   fChain->SetBranchAddress("Muon_emVetoEt05", &Muon_emVetoEt05, &b_Muon_emVetoEt05);
   fChain->SetBranchAddress("Muon_trackerVetoPt05", &Muon_trackerVetoPt05, &b_Muon_trackerVetoPt05);
   fChain->SetBranchAddress("PVCollection_Size", &PVCollection_Size, &b_PVCollection_Size);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_NTracks", &PV_NTracks, &b_PV_NTracks);
   fChain->SetBranchAddress("TripletCollectionSize", &TripletCollectionSize, &b_TripletCollectionSize);
   fChain->SetBranchAddress("Mu1_Pt", &Mu1_Pt, &b_Mu1_Pt);
   fChain->SetBranchAddress("Mu1_Eta", &Mu1_Eta, &b_Mu1_Eta);
   fChain->SetBranchAddress("Mu1_Phi", &Mu1_Phi, &b_Mu1_Phi);
   fChain->SetBranchAddress("Mu1_TripletIndex", &Mu1_TripletIndex, &b_Mu1_TripletIndex);
   fChain->SetBranchAddress("Mu2_Pt", &Mu2_Pt, &b_Mu2_Pt);
   fChain->SetBranchAddress("Mu2_Eta", &Mu2_Eta, &b_Mu2_Eta);
   fChain->SetBranchAddress("Mu2_Phi", &Mu2_Phi, &b_Mu2_Phi);
   fChain->SetBranchAddress("Mu2_TripletIndex", &Mu2_TripletIndex, &b_Mu2_TripletIndex);
   fChain->SetBranchAddress("Mu3_Pt", &Mu3_Pt, &b_Mu3_Pt);
   fChain->SetBranchAddress("Mu3_Eta", &Mu3_Eta, &b_Mu3_Eta);
   fChain->SetBranchAddress("Mu3_Phi", &Mu3_Phi, &b_Mu3_Phi);
   fChain->SetBranchAddress("Mu3_TripletIndex", &Mu3_TripletIndex, &b_Mu3_TripletIndex);
   fChain->SetBranchAddress("GenMatchMu1_SimPt", &GenMatchMu1_SimPt, &b_GenMatchMu1_SimPt);
   fChain->SetBranchAddress("GenMatchMu2_SimPt", &GenMatchMu2_SimPt, &b_GenMatchMu2_SimPt);
   fChain->SetBranchAddress("GenMatchMu3_SimPt", &GenMatchMu3_SimPt, &b_GenMatchMu3_SimPt);
   fChain->SetBranchAddress("GenMatchMu1_SimEta", &GenMatchMu1_SimEta, &b_GenMatchMu1_SimEta);
   fChain->SetBranchAddress("GenMatchMu2_SimEta", &GenMatchMu2_SimEta, &b_GenMatchMu2_SimEta);
   fChain->SetBranchAddress("GenMatchMu3_SimEta", &GenMatchMu3_SimEta, &b_GenMatchMu3_SimEta);
   fChain->SetBranchAddress("GenMatchMu1_SimPhi", &GenMatchMu1_SimPhi, &b_GenMatchMu1_SimPhi);
   fChain->SetBranchAddress("GenMatchMu2_SimPhi", &GenMatchMu2_SimPhi, &b_GenMatchMu2_SimPhi);
   fChain->SetBranchAddress("GenMatchMu3_SimPhi", &GenMatchMu3_SimPhi, &b_GenMatchMu3_SimPhi);
   fChain->SetBranchAddress("GenMatchMu1_Pt", &GenMatchMu1_Pt, &b_GenMatchMu1_Pt);
   fChain->SetBranchAddress("GenMatchMu2_Pt", &GenMatchMu2_Pt, &b_GenMatchMu2_Pt);
   fChain->SetBranchAddress("GenMatchMu3_Pt", &GenMatchMu3_Pt, &b_GenMatchMu3_Pt);
   fChain->SetBranchAddress("GenMatchMu1_Eta", &GenMatchMu1_Eta, &b_GenMatchMu1_Eta);
   fChain->SetBranchAddress("GenMatchMu2_Eta", &GenMatchMu2_Eta, &b_GenMatchMu2_Eta);
   fChain->SetBranchAddress("GenMatchMu3_Eta", &GenMatchMu3_Eta, &b_GenMatchMu3_Eta);
   fChain->SetBranchAddress("GenMatchMu1_Phi", &GenMatchMu1_Phi, &b_GenMatchMu1_Phi);
   fChain->SetBranchAddress("GenMatchMu2_Phi", &GenMatchMu2_Phi, &b_GenMatchMu2_Phi);
   fChain->SetBranchAddress("GenMatchMu3_Phi", &GenMatchMu3_Phi, &b_GenMatchMu3_Phi);
   fChain->SetBranchAddress("TripletVtx_x", &TripletVtx_x, &b_TripletVtx_x);
   fChain->SetBranchAddress("TripletVtx_y", &TripletVtx_y, &b_TripletVtx_y);
   fChain->SetBranchAddress("TripletVtx_z", &TripletVtx_z, &b_TripletVtx_z);
   fChain->SetBranchAddress("TripletVtx_Chi2", &TripletVtx_Chi2, &b_TripletVtx_Chi2);
   fChain->SetBranchAddress("TripletVtx_NDOF", &TripletVtx_NDOF, &b_TripletVtx_NDOF);
   fChain->SetBranchAddress("Triplet_Mass", &Triplet_Mass, &b_Triplet_Mass);
   fChain->SetBranchAddress("Triplet_Pt", &Triplet_Pt, &b_Triplet_Pt);
   fChain->SetBranchAddress("Triplet_Eta", &Triplet_Eta, &b_Triplet_Eta);
   fChain->SetBranchAddress("Triplet_Phi", &Triplet_Phi, &b_Triplet_Phi);
   fChain->SetBranchAddress("Triplet_Charge", &Triplet_Charge, &b_Triplet_Charge);
   fChain->SetBranchAddress("RefittedPV_x", &RefittedPV_x, &b_RefittedPV_x);
   fChain->SetBranchAddress("RefittedPV_y", &RefittedPV_y, &b_RefittedPV_y);
   fChain->SetBranchAddress("RefittedPV_z", &RefittedPV_z, &b_RefittedPV_z);
   fChain->SetBranchAddress("RefittedPV_NTracks", &RefittedPV_NTracks, &b_RefittedPV_NTracks);
   fChain->SetBranchAddress("RefittedPV_isValid", &RefittedPV_isValid, &b_RefittedPV_isValid);
   fChain->SetBranchAddress("FlightDistPVSV", &FlightDistPVSV, &b_FlightDistPVSV);
   fChain->SetBranchAddress("FlightDistPVSV_Err", &FlightDistPVSV_Err, &b_FlightDistPVSV_Err);
   fChain->SetBranchAddress("FlightDistPVSV_Significance", &FlightDistPVSV_Significance, &b_FlightDistPVSV_Significance);
   fChain->SetBranchAddress("FlightDistPVSV_chi2", &FlightDistPVSV_chi2, &b_FlightDistPVSV_chi2);
   Notify();
}

Bool_t ntupleClass_2017F::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntupleClass_2017F::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntupleClass_2017F::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

//When run inside an event, given the characteristics of a muon (pt, eta, phi), the function return the # of the corresponding muon in the event

//                cout << "Muone trovato n. " << l << endl;
//                cout << "Pt richiesto " << pt << " | Pt trovato " << MuonPt->at(k) << endl;
//                cout << "Eta richiesto " << eta << " | Eta trovato " << MuonEta->at(k) << endl;
//                cout << "Phi richiesto " << phi << " | Phi trovato " << MuonPhi->at(k) << endl;
//                cout << "GLnormChi2 -> " << Muon_GLnormChi2->at(k) << endl;
//                cout << "OuterTrackNormChi2 -> " << Muon_outerTrack_normalizedChi2->at(k) << endl;
//                cout << "InnerTrackNormChi2 -> " << Muon_innerTrack_normalizedChi2->at(k) << endl;
//                cout << "nMatches -> " << Muon_numberOfMatches->at(k) << endl << endl;

// * PRIORITY:
// *** if (outerTrack != -999)
//      -> Vedi chi ha l'outerTrackChi2 minore
//      -> A parità di outerTrackChi2:
//          -> guarda chi ha l'innerTrackChi2 minore
//              -> guarda chi ha il numero di matches maggiore

// *** if (there is NOT at least one 1 track != -999)
//      -> Vedi chi ha l'innerTrackChi2 minore
//      -> A parità di Chi2: guarda chi ha il numero di matches maggiore

Double_t ntupleClass_2017F::MuonFinder(Double_t pt, Double_t eta, Double_t phi){
    int n=0, m=0, badOuterChi2=0;
    for(int g=0; g<MuonPt->size(); g++){
        if(pt == MuonPt->at(g) && eta == MuonEta->at(g) && phi == MuonPhi->at(g)){
            n++;
            m = g;
            if(Muon_outerTrack_normalizedChi2->at(g) == -999)
                badOuterChi2++;
        }
    }
    
    // MULTIPLE muon section
    if(n>1) {  //cout << "Error: There is more than 1 muon that matches the conditions!" << endl;
        
        
        // ####### CASO 0:
        // NON c'è NESSUN muone tra quelli uguali che ha outerTrachChi2 != -999
        //implementa controllo sul Chi2 dell'inner (per gli n muoni):
        // * Se min unico -> ritorna il suo indice
        // * Se min uguali  -> implementa controllo sul numero matches e restituisci l'indice di quello con più matches
        // * Se n. matches uguale -> stampa un messaggio di errore
        
        if((n-badOuterChi2) == 0){
            double Chi2InnerTrackmin[2] = {0};
            // Chi2InnerTrackmin[1] è il valore del Chi2,
            // Chi2InnerTrackmin[0] è l'indice corrispondente
            Chi2InnerTrackmin[1] = 999;
            // Cerca il valore MIN di Chi2 dell'innerTrack
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_innerTrack_normalizedChi2->at(k) < Chi2InnerTrackmin[1]){
                    Chi2InnerTrackmin[1] = Muon_innerTrack_normalizedChi2->at(k);
                    Chi2InnerTrackmin[0] = k;
                }
            }
            // Conta il # di muoni con questo valor MIN
            int nMuonChi2InnerMin = 0;
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k)){
                    nMuonChi2InnerMin++;
                }
            }
            // Se è uno solo, restituisci il suo indice
            if(nMuonChi2InnerMin == 1){
                return Chi2InnerTrackmin[0];
            }
            // Altrimenti cerca il valore MAX del numero dei matches
            else {
                double NumberOfMatches[2] = {0}; // stessa idea di Chi2InnerTrackmin[]
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k)){
                        if(Muon_numberOfMatches->at(k) > NumberOfMatches[1]){
                            NumberOfMatches[1] = Muon_numberOfMatches->at(k);
                            NumberOfMatches[0] = k;
                        }
                    }
                }
                // Conta il # di muoni con questo valor MAX
                int nMuonMatchesMax = 0;
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k) && NumberOfMatches[1] == Muon_numberOfMatches->at(k)){
                        nMuonMatchesMax++;
                    }
                }
                // Se è uno solo, restituisci il suo indice
                if(nMuonMatchesMax == 1){
                    return NumberOfMatches[0];
                }
                // Altrimenti stampa un messaggio di errore e restituisci l'indice di uno qualsiasi (in questo caso, l'ultimo)
                else{
                    cout << "Non so quale scegliere tra muoni multipli con outerChi2 = -999! " << endl;
                    return NumberOfMatches[0];
                }
            }
        }
        
        
        // ####### CASO 1:
        // C'è solo 1 muone tra tutti quelli uguali che ha outerTrachChi2 != -999
        //controlla tra i muoni uguali qual è quello con outerTrachChi2 != -999 e restituisci il suo indice
        else if((n-badOuterChi2) == 1){
            int indexGoodMu = -1;
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) != -999)
                    indexGoodMu = k;
            }
            if (indexGoodMu == -1)
                cout << "C'è sicuramente un BUG !!!" << endl;
            
            return indexGoodMu;
        }
        
        
        
        // ####### CASO 2:
        // C'è più di 1 muone con outerTrachChi2 != -999
        //confronta Chi2 dell'outer su tutti i muoni uguali (n-badOuterChi2) e ritorna l'indice di quello minore.
        // Se ce ne sono più di uno uguale al minimo, confronta il # di matches e restituisci l'indice di quello con numero maggiore
        // Se continuano ancora ad essere uguali, stampa un messaggio di errore
        
        else if((n-badOuterChi2) > 1){
            double Chi2OuterTrackmin[2] = {0};
            Chi2OuterTrackmin[1] = 999;
            
            // Cerca il valore MIN di Chi2 dell'outerTrack
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) < Chi2OuterTrackmin[1]){
                    Chi2OuterTrackmin[1] = Muon_outerTrack_normalizedChi2->at(k);
                    Chi2OuterTrackmin[0] = k;
                }
            }
            // Conta il # di muoni con questo valor MIN
            int nMuonChi2OuterMin = 0;
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2OuterTrackmin[1] == Muon_outerTrack_normalizedChi2->at(k)){
                    nMuonChi2OuterMin++;
                }
            }
            // Se è uno solo, restituisci il suo indice
            if(nMuonChi2OuterMin == 1){
                return Chi2OuterTrackmin[0];
            }
            // Altrimenti :
            else{
                // cerca il valore MIN di Chi2 dell'innerTrack
                double Chi2InnerTrackmin[2] = {0};
                Chi2InnerTrackmin[1] = 999;
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) == Chi2OuterTrackmin[1] && Muon_innerTrack_normalizedChi2->at(k) < Chi2InnerTrackmin[1]){
                        Chi2InnerTrackmin[1] = Muon_innerTrack_normalizedChi2->at(k);
                        Chi2InnerTrackmin[0] = k;
                    }
                }
                // Conta il # di muoni con questo valor MIN
                int nMuonChi2InnerMin = 0;
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) == Chi2OuterTrackmin[1] && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k)){
                        nMuonChi2InnerMin++;
                    }
                }
                // Se è uno solo, restituisci il suo indice
                if(nMuonChi2InnerMin == 1){
                    return Chi2InnerTrackmin[0];
                }
                // Altrimenti cerca il valore MAX del numero dei matches
                else {
                    double NumberOfMatches[2] = {0};
                    for (int k=0; k<MuonPt->size(); k++){
                        if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) == Chi2OuterTrackmin[1] && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k)){
                            if(Muon_numberOfMatches->at(k) > NumberOfMatches[1]){
                                NumberOfMatches[1] = Muon_numberOfMatches->at(k);
                                NumberOfMatches[0] = k;
                            }
                        }
                    }
                    // Conta il # di muoni con questo valor MAX
                    int nMuonMatchesMax = 0;
                    for (int k=0; k<MuonPt->size(); k++){
                        if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) == Chi2OuterTrackmin[1] && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k) && NumberOfMatches[1] == Muon_numberOfMatches->at(k)){
                            nMuonMatchesMax++;
                        }
                    }
                    // Se è uno solo, restituisci il suo indice
                    if(nMuonMatchesMax == 1){
                        return NumberOfMatches[0];
                    }
                    // Altrimenti stampa un messaggio di errore e restituisci l'indice di uno qualsiasi (in questo caso, l'ultimo)
                    else{
                        cout << "Non so quale scegliere tra muoni multipli con outerChi2 != -999! " << endl;
                        return NumberOfMatches[0];
                    }
                }
            }
            
        }
        
        // ####### CASO 3: Se non rientra in nessuno dei casi precedenti
        else {
            cout << "ERROR in multiple muon number!" << endl;
            return m;
        }
        
    }
    // end multiple muons section
    
    else
        return m;
}

// Given 2 muons it returns 'true' if DeltaR < DeltaRmax
Bool_t ntupleClass_2017F::isDeltaRGood(Float_t eta1, Float_t eta2, Float_t phi1, Float_t phi2, Float_t DeltaRmax){
    float n=0;
    n=TMath::Sqrt(pow((eta1-eta2),2)+pow((phi1-phi2),2));
    if(n<DeltaRmax) return true;
    else return false;
}

// Given 2 muons it returns 'true' if |DeltaZ| < DeltaZmax
Bool_t ntupleClass_2017F::isDeltaZGood(Float_t vz1, Float_t vz2, Float_t DeltaZmax){
    float n = TMath::Abs(vz2 - vz1);
    if(n<DeltaZmax) return true;
    else return false;
}

// Given the 'DeltaZ' value of each of the 3 muons of the triplet & the DeltaZmax value, the function returns true if all of the 3 possible pairs of muons satisfy isDeltaZGood
Bool_t ntupleClass_2017F::isPairDeltaZGood(Float_t DeltaZ1, Float_t DeltaZ2, Float_t DeltaZ3, Float_t DeltaZmax){
    if(isDeltaZGood(DeltaZ1, DeltaZ2, DeltaZmax) == true && isDeltaZGood(DeltaZ2, DeltaZ3, DeltaZmax) == true && isDeltaZGood(DeltaZ1, DeltaZ3, DeltaZmax) == true)
        return true;
    else return false;
}
//Given the "refering number" of the triplet & the DeltaZmax value, the function returns true if all of the 3 possible pairs of muons satisfy isDeltaRGood
Bool_t ntupleClass_2017F::isPairDeltaRGood(Int_t ntriplet, Float_t DeltaRmax){
    if(isDeltaRGood(Mu1_Eta->at(ntriplet), Mu2_Eta->at(ntriplet), Mu1_Phi->at(ntriplet), Mu2_Phi->at(ntriplet), DeltaRmax) == true && isDeltaRGood(Mu2_Eta->at(ntriplet), Mu3_Eta->at(ntriplet), Mu2_Phi->at(ntriplet), Mu3_Phi->at(ntriplet), DeltaRmax) == true && isDeltaRGood(Mu1_Eta->at(ntriplet), Mu3_Eta->at(ntriplet), Mu1_Phi->at(ntriplet), Mu3_Phi->at(ntriplet), DeltaRmax) == true)
        return true;
    else return false;
}

//Given the energy, eta, phi of a muon, the function returns the value of the impulse of the muon
Double_t ntupleClass_2017F::MuonP(Double_t pt, Double_t eta, Double_t phi, Double_t mass){
    TLorentzVector mu;
    mu.SetPtEtaPhiM(pt, eta, phi, mass);
    return mu.P();
}

//Given the energy, eta, phi, pt of 2 muons, if their charge is opposite the function returns their invariant mass, otherwise it returns 0
Double_t ntupleClass_2017F::DimuonMass(Double_t charge1, Double_t charge2, Double_t pt1, Double_t pt2, Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2, Double_t mass){
    double inv = 0;
    if(charge1 - charge2 != 0){
        return inv;
    }
    else {
        TLorentzVector mu1, mu2, mutot;
        mu1.SetPtEtaPhiM(pt1, eta1, phi1, mass);
        mu2.SetPtEtaPhiM(pt2, eta2, phi2, mass);
        mutot = mu1 + mu2;
        return mutot.M();
    }
}

//Given the energy, eta, phi, pt of 4 muons it returns their invariant mass
Float_t ntupleClass_2017F::QuadMuonMass(Float_t pt1, Float_t pt2, Float_t pt3, Float_t pt4, Float_t eta1, Float_t eta2, Float_t eta3, Float_t eta4, Float_t phi1, Float_t phi2, Float_t phi3, Float_t phi4, Double_t mass){
    TLorentzVector mu1, mu2, mu3, mu4, mutot;
    mu1.SetPtEtaPhiM(pt1, eta1, phi1, mass);
    mu2.SetPtEtaPhiM(pt2, eta2, phi2, mass);
    mu3.SetPtEtaPhiM(pt3, eta3, phi3, mass);
    mu4.SetPtEtaPhiM(pt4, eta4, phi4, mass);
    mutot = mu1 + mu2 + mu3 + mu4;
    return mutot.M();
}

// Given the dimuon mass and a value for sigma, the function return true if the dimuon mass is NOT compatible w/ the mass of a Phi(1020) in 2 sigmas
Bool_t ntupleClass_2017F::isNotAPhi(Double_t dimumass, Double_t sigma){
    const double PhiMass = 1.019461; //Phi Mass in GeV
    if (dimumass < (PhiMass-2*sigma) || dimumass > (PhiMass+2*sigma))
        return true;
    else return false;
}

//Given 3 muons it checks, for all the pairs having opposite charges, if the corresponding dimuon mass is NOT compatible w/ the Phi mass(1020)
Bool_t ntupleClass_2017F::isPairNotAPhi(Double_t dimu1_2, Double_t dimu2_3, Double_t dimu1_3, Double_t sigma){
    if(dimu1_2!=0 && dimu2_3!=0 && dimu1_3!=0){
        cout << "ERROR!!! There are 3 pairs of muons w/ opposite charges !!!" << endl;
        return false;
    }
    else if (dimu1_2!=0 && dimu2_3!=0 && dimu1_3==0){
        if (isNotAPhi(dimu1_2,sigma) == true && isNotAPhi(dimu2_3,sigma) == true)
            return true;
        else return false;
    }
    else if(dimu1_2!=0 && dimu2_3==0 && dimu1_3!=0){
        if (isNotAPhi(dimu1_2,sigma) == true && isNotAPhi(dimu1_3,sigma) == true)
            return true;
        else return false;
    }
    else if(dimu1_2==0 && dimu2_3!=0 && dimu1_3!=0){
        if (isNotAPhi(dimu2_3,sigma) == true && isNotAPhi(dimu1_3,sigma) == true)
            return true;
        else return false;
    }
    else if(dimu1_2!=0 && dimu2_3==0 && dimu1_3==0){
        if (isNotAPhi(dimu1_2,sigma) == true)
            return true;
        else return false;
    }
    else if(dimu1_2==0 && dimu2_3!=0 && dimu1_3==0){
        if (isNotAPhi(dimu2_3,sigma) == true)
            return true;
        else return false;
    }
    else if(dimu1_2==0 && dimu2_3==0 && dimu1_3!=0){
        if (isNotAPhi(dimu1_3,sigma) == true)
            return true;
        else return false;
    }
    cout << "C'è qualquadra che non cosa!" << endl;
    return false;
}

// Given a vector of strings the function initializes it with the names of particles
void ntupleClass_2017F::particleName(TString pId[37]){
    //Particle name list
    pId[0] = "#mu^{-}";
    pId[1] = "#mu^{+}";
    pId[2] = "#tau^{-}";
    pId[3] = "#tau^{+}";
    pId[4] = "e^{-}";
    pId[5] = "e^{+}";
    pId[6] = "D^{+}";
    pId[7] = "D^{-}";
    pId[8] = "D^{0}";
    pId[9] = "#bar{D}^{0}";
    pId[10] = "K^{+}";
    pId[11] = "K^{-}";
    pId[12] = "B^{+}_{c}";
    pId[13] = "B^{-}_{c}";
    pId[14] = "B^{0}_{s}";
    pId[15] = "#bar{B}^{0}_{s}";
    pId[16] = "#Pi^{+}";
    pId[17] = "#Pi^{-}";
    pId[18] = "D^{+}_{s}";
    pId[19] = "D^{-}_{s}";
    pId[20] = "B^{0}";
    pId[21] = "#bar{B}^{0}";
    pId[22] = "B^{+}";
    pId[23] = "B^{-}";
    pId[24] = "p";
    pId[25] = "#Xi^{-}_{b}";
    pId[26] = "#Xi^{+}_{b}";
    pId[27] = "K^{0}_{L}";
    pId[28] = "#omega"; // omega(782)
    pId[29] = "n";
    pId[30] = "#Lambda^{+}_{c}";
    pId[31] = "#Lambda^{-}_{c}";
    pId[32] = "#Lambda^{0}_{b}";
    pId[33] = "#bar{#Lambda}^{0}_{b}";
    pId[34] = "#bar{p}";
    pId[35] = "#Phi"; //Phi(1020)
    pId[36] = "NN";
}

// Given the particle Id, it checks at which particle corresponds and increases the relative counter
void ntupleClass_2017F::particleId(Int_t pdgId, Int_t IdSummary[37]){
    
    switch (pdgId) {
        case 13: // is a mu-
            IdSummary[0]++;
            break;
        case -13: // is a mu+
            IdSummary[1]++;
            break;
        case 15: // is a tau-
            IdSummary[2]++;
            break;
        case -15: // is a tau+
            IdSummary[3]++;
            break;
        case 11: // is a e-
            IdSummary[4]++;
            break;
        case -11: // is a e+
            IdSummary[5]++;
            break;
        case 411: // is a D+
            IdSummary[6]++;
            break;
        case -411: // is a D-
            IdSummary[7]++;
            break;
        case 421: // is a D0
            IdSummary[8]++;
            break;
        case -421: // is an antiD0
            IdSummary[9]++;
            break;
        case 321: // is a K+
            IdSummary[10]++;
            break;
        case -321: // is a K-
            IdSummary[11]++;
            break;
        case 541: // is a B+_c
            IdSummary[12]++;
            break;
        case -541: // is a B-_c
            IdSummary[13]++;
            break;
        case 531: // is a B0_s
            IdSummary[14]++;
            break;
        case -531: // is a antiB0_s
            IdSummary[15]++;
            break;
        case 211: // is a Pi+
            IdSummary[16]++;
            break;
        case -211: // is a Pi-
            IdSummary[17]++;
            break;
        case 431: // is a Ds+
            IdSummary[18]++;
            break;
        case -431: // is a Ds-
            IdSummary[19]++;
            break;
        case 511: // is a B0
            IdSummary[20]++;
            break;
        case -511: // is an antiB0
            IdSummary[21]++;
            break;
        case 521: // is a B+
            IdSummary[22]++;
            break;
        case -521: // is a B-
            IdSummary[23]++;
            break;
        case 2212: // is a p
            IdSummary[24]++;
            break;
        case 5132: // is a Csi-b
            IdSummary[25]++;
            break;
        case -5132: // is a Csi+b
            IdSummary[26]++;
            break;
        case 130: // is a K0_L
            IdSummary[27]++;
            break;
        case 223: // is a omega(782)
            IdSummary[28]++;
            break;
        case 2112: // is a n
            IdSummary[29]++;
            break;
        case 4122: // is a Lambda+_c
            IdSummary[30]++;
            break;
        case -4122: // is a Lambda-_c
            IdSummary[31]++;
            break;
        case 5122: // is a Lambda0_b
            IdSummary[32]++;
            break;
        case -5122: // is an antiLambda0_b
            IdSummary[33]++;
            break;
        case -2212: // is an antip
            IdSummary[34]++;
            break;
        case 333: // is a Phi(1020)
            IdSummary[35]++;
            break;
        case 0: // is a UNKNOWN
            IdSummary[36]++;
            break;
            
        default: cout << "Unknown particle! PdgId = " << pdgId << endl;
            break;
    }
}

// Given a vector of strings the function initializes it with the names of the cuts
void ntupleClass_2017F::cutName(TString listCut[7]){
    //    List of the names of the cuts
    listCut[0] = "BeforeCuts";
    listCut[1] = "Chi2Triplet";
    listCut[2] = "MassTriplet";
    listCut[3] = "2glb+1trk";
    listCut[4] = "DeltaR";
    listCut[5] = "DeltaZ";
    listCut[6] = "DimuonMass-Phi";
}

// Given a vector of strings the function initializes it with the names of the cuts
void ntupleClass_2017F::cutName_New(TString listCut[7]){
    //    List of the names of the cuts
    listCut[0] = "BeforeCuts";
    listCut[1] = "Chi2Triplet";
    listCut[2] = "2glb+1trk";
    listCut[3] = "MassTriplet";
    listCut[4] = "DeltaR";
    listCut[5] = "DeltaZ";
    listCut[6] = "DimuonMass-Phi";
}

// Given a vector of strings the function initializes it with the names of the cuts
void ntupleClass_2017F::cutName2(TString listCut[7]){
    //    List of the names of the cuts
    listCut[0] = "BeforeCuts";
    listCut[1] = "";
    listCut[2] = "";
    listCut[3] = "";
    listCut[4] = "";
    listCut[5] = "";
    listCut[6] = "";
}

// This function writes on the canvas the histo w/ the particle names
void ntupleClass_2017F::PdgIdHisto(TCanvas *canv, TH1I *hist, Int_t Idsummary[37], TString pIdList[37]){
    for(int k=0; k<37; k++){
        hist->Fill(k, Idsummary[k]);
        hist->GetXaxis()->SetBinLabel(k+1, pIdList[k]);
    }
    canv->SetLogy();
    hist->DrawCopy();
    hist->DrawCopy("HIST TEXT0");
    canv->Write();
    canv->Close();
}

// This function writes on the canvas the histo of the cuts efficiency
void ntupleClass_2017F::CutEffHisto(TCanvas *canv, TH1I *hist, Int_t cut[7], TString listCut[7]){
    for(int k=0; k<7; k++){
        hist->Fill(k+1, cut[k]);
        hist->GetXaxis()->SetBinLabel(k+1, listCut[k]);
    }
    canv->SetLogy();
    hist->DrawCopy();
    hist->DrawCopy("HIST TEXT0");
    canv->Write();
    canv->Close();
}


// This function writes on the canvas the histo of the cuts efficiency
void ntupleClass_2017F::DiMuonHisto(TH1D *hist, Double_t dimu1_2, Double_t dimu2_3, Double_t dimu1_3){
    if(dimu1_2 != 0) hist->Fill(dimu1_2);
    if(dimu2_3 != 0) hist->Fill(dimu2_3);
    if(dimu1_3 != 0) hist->Fill(dimu1_3);
}




#endif // #ifdef ntupleClass_2017F_cxx
