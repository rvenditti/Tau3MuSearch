//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar 29 16:32:34 2019 by ROOT version 6.12/06
// from TTree ntuple/LFVTau ntuple
// found on file: TreeDsPhiPi.root
//////////////////////////////////////////////////////////

#ifndef ntupleClass_Control_h
#define ntupleClass_Control_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

using namespace std;

class ntupleClass_Control {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TString    fileName;//output filename given in Analysis.cpp

    // Fixed size dimensions of array or collections stored in the TTree if any.
    
    // Declaration of leaf types
    UInt_t          evt;
    UInt_t          run;
    UInt_t          lumi;
    vector<int>     *NGoodTriplets;
    UInt_t          nPileUpInt;
    vector<int>     *GenParticle_PdgId;
    vector<double>  *GenParticle_Pt;
    vector<double>  *GenParticle_Eta;
    vector<double>  *GenParticle_Phi;
    vector<int>     *GenParticle_isDs;
    vector<int>     *GenParticle_isB;
    vector<int>     *GenParticle_isBdecay;
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
    vector<double>  *Muon_timeAtIpInOutErr;
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
    vector<double>  *Track_pt;
    vector<double>  *Track_eta;
    vector<double>  *Track_phi;
    vector<double>  *Track_normalizedChi2;
    vector<double>  *Track_numberOfValidHits;
    vector<double>  *Track_charge;
    vector<double>  *Track_dxy;
    vector<double>  *Track_dxyError;
    vector<double>  *Track_dz;
    vector<double>  *Track_dzError;
    vector<double>  *Track_vx;
    vector<double>  *Track_vy;
    vector<double>  *Track_vz;
    Int_t           PVCollection_Size;
    vector<double>  PV_x;
    vector<double>  PV_y;
    vector<double>  PV_z;
    vector<double>  PV_NTracks;
    vector<string>  *Trigger_l1name;
    vector<int>     *Trigger_l1decision;
    vector<int>     *Trigger_l1prescale;
    vector<string>  *Trigger_hltname;
    vector<int>     *Trigger_hltdecision;
    Int_t           TripletCollectionSize2;
    Int_t           SelectedTripletsSize;
    vector<double>  *Mu01_Pt;
    vector<double>  *Mu01_Eta;
    vector<double>  *Mu01_Phi;
    vector<float>   *Mu01_dRtriggerMatch;
    vector<float>   *Mu1_dRtriggerMatch_Mu7;
    vector<float>   *Mu1_dRtriggerMatch_Mu8;
    vector<float>   *Mu1_dRtriggerMatch_Mu8_IP5;
    vector<float>   *Mu1_dRtriggerMatch_Mu8_IP6;
    vector<float>   *Mu1_dRtriggerMatch_Mu9_IP0;
    vector<float>   *Mu1_dRtriggerMatch_Mu9_IP3;
    vector<float>   *Mu1_dRtriggerMatch_Mu9_IP4;
    vector<float>   *Mu1_dRtriggerMatch_Mu9_IP5;
    vector<float>   *Mu1_dRtriggerMatch_Mu9_IP6;
    vector<float>   *Mu1_dRtriggerMatch_Mu12_IP6;
    vector<int>     *Mu01_TripletIndex;
    vector<double>  *Mu02_Pt;
    vector<double>  *Mu02_Eta;
    vector<double>  *Mu02_Phi;
    vector<float>   *Mu02_dRtriggerMatch;
    vector<float>   *Mu2_dRtriggerMatch_Mu7;
    vector<float>   *Mu2_dRtriggerMatch_Mu8;
    vector<int>     *Mu02_TripletIndex;
    vector<double>  *Tr_Pt;
    vector<double>  *Tr_Eta;
    vector<double>  *Tr_Phi;
    vector<float>   *Tr_dRtriggerMatch;
    vector<int>     *Tr_TripletIndex;
    vector<int>     *selectedTripletsIndex;
    vector<double>  *GenMatchMu01_SimPt;
    vector<double>  *GenMatchMu02_SimPt;
    vector<double>  *GenMatchMu01_SimEta;
    vector<double>  *GenMatchMu02_SimEta;
    vector<double>  *GenMatchMu01_SimPhi;
    vector<double>  *GenMatchMu02_SimPhi;
    vector<double>  *GenMatchMu01_Pt;
    vector<double>  *GenMatchMu02_Pt;
    vector<double>  *GenMatchMu01_Eta;
    vector<double>  *GenMatchMu02_Eta;
    vector<double>  *GenMatchMu01_Phi;
    vector<double>  *GenMatchMu02_Phi;
    vector<double>  *Triplet_mindca_iso;
    vector<double>  *Triplet_relativeiso;
    vector<double>  *TripletVtx2_x;
    vector<double>  *TripletVtx2_y;
    vector<double>  *TripletVtx2_z;
    vector<double>  *TripletVtx2_Chi2;
    vector<double>  *TripletVtx2_NDOF;
    vector<double>  *Triplet2_Mass;
    vector<double>  *Triplet2_Pt;
    vector<double>  *Triplet2_Eta;
    vector<double>  *Triplet2_Phi;
    vector<double>  *Triplet2_Charge;
    vector<double>  *RefittedPV2_x;
    vector<double>  *RefittedPV2_y;
    vector<double>  *RefittedPV2_z;
    vector<double>  *RefittedPV2_NTracks;
    vector<int>     *RefittedPV2_isValid;
    vector<double>  *FlightDistPVSV2;
    vector<double>  *FlightDistPVSV2_Err;
    vector<double>  *FlightDistPVSV2_Significance;
    vector<double>  *FlightDistPVSV2_chi2;
    vector<double>  *dxy_mu1;
    vector<double>  *dxy_mu2;
    vector<double>  *dxy_mu3;
    vector<double>  *dxyErr_mu1;
    vector<double>  *dxyErr_mu2;
    vector<double>  *dxyErr_mu3;
    
    // List of branches
    TBranch        *b_evt;   //!
    TBranch        *b_run;   //!
    TBranch        *b_lumi;   //!
    TBranch        *b_NGoodTriplets;   //!
    TBranch        *b_nPileUpInt;   //!
    TBranch        *b_GenParticle_PdgId;   //!
    TBranch        *b_GenParticle_Pt;   //!
    TBranch        *b_GenParticle_Eta;   //!
    TBranch        *b_GenParticle_Phi;   //!
    TBranch        *b_GenParticle_isDs;   //!
    TBranch        *b_GenParticle_isB;   //!
    TBranch        *b_GenParticle_isBdecay;   //!
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
    TBranch        *b_Muon_timeAtIpInOutErr;   //!
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
    TBranch        *b_Track_pt;   //!
    TBranch        *b_Track_eta;   //!
    TBranch        *b_Track_phi;   //!
    TBranch        *b_Track_normalizedChi2;   //!
    TBranch        *b_Track_numberOfValidHits;   //!
    TBranch        *b_Track_charge;   //!
    TBranch        *b_Track_dxy;   //!
    TBranch        *b_Track_dxyError;   //!
    TBranch        *b_Track_dz;   //!
    TBranch        *b_Track_dzError;   //!
    TBranch        *b_Track_vx;   //!
    TBranch        *b_Track_vy;   //!
    TBranch        *b_Track_vz;   //!
    TBranch        *b_PVCollection_Size;   //!
    TBranch        *b_PV_x;   //!
    TBranch        *b_PV_y;   //!
    TBranch        *b_PV_z;   //!
    TBranch        *b_PV_NTracks;   //!
    TBranch        *b_Trigger_l1name;   //!
    TBranch        *b_Trigger_l1decision;   //!
    TBranch        *b_Trigger_l1prescale;   //!
    TBranch        *b_Trigger_hltname;   //!
    TBranch        *b_Trigger_hltdecision;   //!
    TBranch        *b_TripletCollectionSize2;   //!
    TBranch        *b_SelectedTripletsSize;   //!
    TBranch        *b_Mu01_Pt;   //!
    TBranch        *b_Mu01_Eta;   //!
    TBranch        *b_Mu01_Phi;   //!
    TBranch        *b_Mu01_dRtriggerMatch;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu7;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu8;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu8_IP5;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu8_IP6;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu9_IP0;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu9_IP3;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu9_IP4;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu9_IP5;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu9_IP6;   //!
    TBranch        *b_Mu1_dRtriggerMatch_Mu12_IP6;   //!
    TBranch        *b_Mu01_TripletIndex;   //!
    TBranch        *b_Mu02_Pt;   //!
    TBranch        *b_Mu02_Eta;   //!
    TBranch        *b_Mu02_Phi;   //!
    TBranch        *b_Mu02_dRtriggerMatch;   //!
    TBranch        *b_Mu2_dRtriggerMatch_Mu7;   //!
    TBranch        *b_Mu2_dRtriggerMatch_Mu8;   //!
    TBranch        *b_Mu02_TripletIndex;   //!
    TBranch        *b_Tr_Pt;   //!
    TBranch        *b_Tr_Eta;   //!
    TBranch        *b_Tr_Phi;   //!
    TBranch        *b_Tr_dRtriggerMatch;   //!
    TBranch        *b_Tr_TripletIndex;   //!
    TBranch        *b_selectedTripletsIndex;   //!
    TBranch        *b_GenMatchMu01_SimPt;   //!
    TBranch        *b_GenMatchMu02_SimPt;   //!
    TBranch        *b_GenMatchMu01_SimEta;   //!
    TBranch        *b_GenMatchMu02_SimEta;   //!
    TBranch        *b_GenMatchMu01_SimPhi;   //!
    TBranch        *b_GenMatchMu02_SimPhi;   //!
    TBranch        *b_GenMatchMu01_Pt;   //!
    TBranch        *b_GenMatchMu02_Pt;   //!
    TBranch        *b_GenMatchMu01_Eta;   //!
    TBranch        *b_GenMatchMu02_Eta;   //!
    TBranch        *b_GenMatchMu01_Phi;   //!
    TBranch        *b_GenMatchMu02_Phi;   //!
    TBranch        *b_Triplet_mindca_iso;   //!
    TBranch        *b_Triplet_relativeiso;   //!
    TBranch        *b_TripletVtx2_x;   //!
    TBranch        *b_TripletVtx2_y;   //!
    TBranch        *b_TripletVtx2_z;   //!
    TBranch        *b_TripletVtx2_Chi2;   //!
    TBranch        *b_TripletVtx2_NDOF;   //!
    TBranch        *b_Triplet2_Mass;   //!
    TBranch        *b_Triplet2_Pt;   //!
    TBranch        *b_Triplet2_Eta;   //!
    TBranch        *b_Triplet2_Phi;   //!
    TBranch        *b_Triplet2_Charge;   //!
    TBranch        *b_RefittedPV2_x;   //!
    TBranch        *b_RefittedPV2_y;   //!
    TBranch        *b_RefittedPV2_z;   //!
    TBranch        *b_RefittedPV2_NTracks;   //!
    TBranch        *b_RefittedPV2_isValid;   //!
    TBranch        *b_FlightDistPVSV2;   //!
    TBranch        *b_FlightDistPVSV2_Err;   //!
    TBranch        *b_FlightDistPVSV2_Significance;   //!
    TBranch        *b_FlightDistPVSV2_chi2;   //!
    TBranch        *b_dxy_mu1;   //!
    TBranch        *b_dxy_mu2;   //!
    TBranch        *b_dxy_mu3;   //!
    TBranch        *b_dxyErr_mu1;   //!
    TBranch        *b_dxyErr_mu2;   //!
    TBranch        *b_dxyErr_mu3;   //!
 
   ntupleClass_Control(TTree *tree, TString fname);
    virtual ~ntupleClass_Control();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    // in "ntupleClass_Control.C"
    virtual void     LoopControl(TString type, TString datasetName);
    // in "Utilities_Control.C"
    virtual Double_t DimuonMass(Double_t charge1, Double_t charge2, Double_t pt1, Double_t pt2, Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2, Double_t en1, Double_t en2);
    virtual void     Fill_MuonAndTrackVariables(Int_t mu_Ind[NTOT], Double_t pt[NTOT], Double_t eta[NTOT], Double_t phi[NTOT]);
    virtual void     Fill_MuonVariablesGen(Int_t muGen[NTOT], Double_t ptGEN[NMU_C], Double_t etaGEN[NMU_C], Double_t phiGEN[NMU_C]);
    virtual void     Fill_MuonVariablesGen_Sim(Int_t muGen[NTOT], Double_t ptSimGEN[NMU_C], Double_t etaSimGEN[NMU_C], Double_t phiSimGEN[NMU_C]);
    virtual Bool_t   DuplicateFinder(Double_t eta1, Double_t eta2, Double_t etaTr, Double_t phi1, Double_t phi2, Double_t phiTr, Double_t pt1, Double_t pt2, Double_t ptTr);
    virtual void     MatchIndex(TString type, Int_t ind, Int_t mu_Ind[NTOT], Int_t mu[NTOT]);
    virtual Double_t MuonFinder(Double_t pt, Double_t eta, Double_t phi);
    virtual Double_t MuonFinderGen(Int_t muind, Double_t pt, Double_t eta, Double_t phi);
    virtual Double_t MuonP(Double_t pt, Double_t eta, Double_t phi);
    virtual Double_t TrackFinder(Double_t pt, Double_t eta, Double_t phi);
    // Functions for the final tree
    virtual Double_t TreeFin_Angle(Int_t ind);
    virtual void     TreeFin_Fill(TTree *tree, Int_t ind, Int_t mu_Ind[NMU], Int_t mu[NMU], Double_t &lumi, Double_t &run, Double_t &evt, Double_t &puFactor, Double_t &nHLT, Double_t &DeltaR_max, Double_t &DeltaZ_max, Double_t &Pmu3, Double_t &cLP, Double_t &tKink, Double_t &segmComp, Double_t &tripletMass, Double_t &tripletMassReso, Double_t &fv_nC, Double_t &fv_dphi3D, Double_t &fv_d3D,  Double_t &fv_d3Dsig, Double_t &d0, Double_t &d0sig, Double_t &mindca_iso, Double_t &trkRel, Double_t &Pmu1, Double_t &Ptmu1, Double_t &etamu1, Double_t &Pmu2, Double_t &Ptmu2, Double_t &etamu2, Double_t &Ptmu3, Double_t &etamu3, Double_t &P_trip, Double_t &Pt_trip, Double_t &eta_trip, Double_t &nStationsMu1, Double_t &nStationsMu2, Double_t &nStationsMu3, Double_t &Iso03Mu1, Double_t &Iso03Mu2, Double_t &Iso03Mu3, Double_t &Iso05Mu1, Double_t &Iso05Mu2, Double_t &Iso05Mu3, Double_t &nMatchesMu1, Double_t &nMatchesMu2, Double_t &nMatchesMu3, Double_t &timeAtIpInOutMu_sig1, Double_t &timeAtIpInOutMu_sig2, Double_t &timeAtIpInOutMu_sig3, Double_t &cQ_uS, Double_t &cQ_tK, Double_t &cQ_gK, Double_t &cQ_tRChi2, Double_t &cQ_sRChi2, Double_t &cQ_Chi2LP, Double_t &cQ_Chi2LM, Double_t &cQ_lD, Double_t &cQ_gDEP, Double_t &cQ_tM, Double_t &cQ_gTP, Double_t &calEn_emMu1, Double_t &calEn_emMu2, Double_t &calEn_emMu3, Double_t &calEn_hadMu1, Double_t &calEn_hadMu2, Double_t &calEn_hadMu3, Double_t &caloComp, Double_t &fliDistPVSV_Chi2, Double_t &isGlb1, Double_t &isTracker1, Double_t &isLoose1, Double_t &isSoft1, Double_t &isPF1, Double_t &isRPC1, Double_t &isSA1, Double_t &isCalo1, Double_t &isGlb2, Double_t &isTracker2, Double_t &isLoose2, Double_t &isSoft2, Double_t &isPF2, Double_t &isRPC2, Double_t &isSA2, Double_t &isCalo2, Double_t &isGlb3, Double_t &isTracker3, Double_t &isLoose3, Double_t &isSoft3, Double_t &isPF3, Double_t &isRPC3, Double_t &isSA3, Double_t &isCalo3, Double_t &Vx1, Double_t &Vx2, Double_t &Vx3, Double_t &Vy1, Double_t &Vy2, Double_t &Vy3, Double_t &Vz1, Double_t &Vz2, Double_t &Vz3, Double_t &RefVx1, Double_t &RefVx2, Double_t &RefVx3, Double_t &RefVy1, Double_t &RefVy2, Double_t &RefVy3, Double_t &RefVz1, Double_t &RefVz2, Double_t &RefVz3, Double_t &SVx, Double_t &SVy, Double_t &SVz, Double_t &had03, Double_t &had05, Double_t &nJets03, Double_t &nJets05, Double_t &nTracks03, Double_t &nTracks05, Double_t &sumPt03, Double_t &sumPt05, Double_t &hadVeto03, Double_t &hadVeto05, Double_t &emVeto03, Double_t &emVeto05, Double_t &trVeto03, Double_t &trVeto05, Double_t &EnMu1, Double_t &EnMu2, Double_t &EnMu3, Double_t &ChargeMu1, Double_t &ChargeMu2, Double_t &ChargeMu3, Double_t &isQValid1, Double_t &isTValid1, Double_t &isIsoValid1, Double_t &GLnormChi2_mu1, Double_t &GL_nValidMuHits1, Double_t &trkLayersWMeas1, Double_t &nValidPixelHits1, Double_t &outerTrk_P_1, Double_t &outerTrk_Eta_1, Double_t &outerTrk_normChi2_1, Double_t &outerTrk_muStValidHits_1, Double_t &innerTrk_P_1, Double_t &innerTrk_Eta_1, Double_t &innerTrk_normChi2_1, Double_t &QInnerOuter_1, Double_t &cQ_uS_1, Double_t &cQ_tK_1, Double_t &cQ_gK_1, Double_t &cQ_tRChi2_1, Double_t &cQ_sRChi2_1, Double_t &cQ_Chi2LP_1, Double_t &cQ_Chi2LM_1, Double_t &cQ_lD_1, Double_t &cQ_gDEP_1, Double_t &cQ_tM_1, Double_t &cQ_gTP_1, Double_t &segmComp_1, Double_t &caloComp_1, Double_t &isQValid2, Double_t &isTValid2, Double_t &isIsoValid2, Double_t &GLnormChi2_mu2, Double_t &GL_nValidMuHits2, Double_t &trkLayersWMeas2, Double_t &nValidPixelHits2, Double_t &outerTrk_P_2, Double_t &outerTrk_Eta_2, Double_t &outerTrk_normChi2_2, Double_t &outerTrk_muStValidHits_2, Double_t &innerTrk_P_2, Double_t &innerTrk_Eta_2, Double_t &innerTrk_normChi2_2, Double_t &QInnerOuter_2, Double_t &cQ_uS_2, Double_t &cQ_tK_2, Double_t &cQ_gK_2, Double_t &cQ_tRChi2_2, Double_t &cQ_sRChi2_2, Double_t &cQ_Chi2LP_2, Double_t &cQ_Chi2LM_2, Double_t &cQ_lD_2, Double_t &cQ_gDEP_2, Double_t &cQ_tM_2, Double_t &cQ_gTP_2, Double_t &segmComp_2, Double_t &caloComp_2, Double_t &isQValid3, Double_t &isTValid3, Double_t &isIsoValid3, Double_t &GLnormChi2_mu3, Double_t &GL_nValidMuHits3, Double_t &trkLayersWMeas3, Double_t &nValidPixelHits3, Double_t &outerTrk_P_3, Double_t &outerTrk_Eta_3, Double_t &outerTrk_normChi2_3, Double_t &outerTrk_muStValidHits_3, Double_t &innerTrk_P_3, Double_t &innerTrk_Eta_3, Double_t &innerTrk_normChi2_3, Double_t &QInnerOuter_3, Double_t &cQ_uS_3, Double_t &cQ_tK_3, Double_t &cQ_gK_3, Double_t &cQ_tRChi2_3, Double_t &cQ_sRChi2_3, Double_t &cQ_Chi2LP_3, Double_t &cQ_Chi2LM_3, Double_t &cQ_lD_3, Double_t &cQ_gDEP_3, Double_t &cQ_tM_3, Double_t &cQ_gTP_3, Double_t &segmComp_3, Double_t &caloComp_3, Double_t &trk_dZ, Double_t &trk_dXY);
    virtual void     TreeFin_Init(TTree *&tree, Double_t &lumi, Double_t &run, Double_t &evt, Double_t &puFactor, Double_t &nHLT, Double_t &DeltaR_max, Double_t &DeltaZ_max, Double_t &Pmu3, Double_t &cLP, Double_t &tKink, Double_t &segmComp, Double_t &tripletMass, Double_t &tripletMassReso, Double_t &fv_nC, Double_t &fv_dphi3D, Double_t &fv_d3D,  Double_t &fv_d3Dsig, Double_t &d0, Double_t &d0sig, Double_t &mindca_iso, Double_t &trkRel, Double_t &Pmu1, Double_t &Ptmu1, Double_t &etamu1, Double_t &Pmu2, Double_t &Ptmu2, Double_t &etamu2, Double_t &Ptmu3, Double_t &etamu3, Double_t &P_trip, Double_t &Pt_trip, Double_t &eta_trip, Double_t &nStationsMu1, Double_t &nStationsMu2, Double_t &nStationsMu3, Double_t &Iso03Mu1, Double_t &Iso03Mu2, Double_t &Iso03Mu3, Double_t &Iso05Mu1, Double_t &Iso05Mu2, Double_t &Iso05Mu3, Double_t &nMatchesMu1, Double_t &nMatchesMu2, Double_t &nMatchesMu3, Double_t &timeAtIpInOutMu_sig1, Double_t &timeAtIpInOutMu_sig2, Double_t &timeAtIpInOutMu_sig3, Double_t &cQ_uS, Double_t &cQ_tK, Double_t &cQ_gK, Double_t &cQ_tRChi2, Double_t &cQ_sRChi2, Double_t &cQ_Chi2LP, Double_t &cQ_Chi2LM, Double_t &cQ_lD, Double_t &cQ_gDEP, Double_t &cQ_tM, Double_t &cQ_gTP, Double_t &calEn_emMu1, Double_t &calEn_emMu2, Double_t &calEn_emMu3, Double_t &calEn_hadMu1, Double_t &calEn_hadMu2, Double_t &calEn_hadMu3, Double_t &caloComp, Double_t &fliDistPVSV_Chi2, Double_t &isGlb1, Double_t &isTracker1, Double_t &isLoose1, Double_t &isSoft1, Double_t &isPF1, Double_t &isRPC1, Double_t &isSA1, Double_t &isCalo1, Double_t &isGlb2, Double_t &isTracker2, Double_t &isLoose2, Double_t &isSoft2, Double_t &isPF2, Double_t &isRPC2, Double_t &isSA2, Double_t &isCalo2, Double_t &isGlb3, Double_t &isTracker3, Double_t &isLoose3, Double_t &isSoft3, Double_t &isPF3, Double_t &isRPC3, Double_t &isSA3, Double_t &isCalo3, Double_t &Vx1, Double_t &Vx2, Double_t &Vx3, Double_t &Vy1, Double_t &Vy2, Double_t &Vy3, Double_t &Vz1, Double_t &Vz2, Double_t &Vz3, Double_t &RefVx1, Double_t &RefVx2, Double_t &RefVx3, Double_t &RefVy1, Double_t &RefVy2, Double_t &RefVy3, Double_t &RefVz1, Double_t &RefVz2, Double_t &RefVz3, Double_t &SVx, Double_t &SVy, Double_t &SVz, Double_t &had03, Double_t &had05, Double_t &nJets03, Double_t &nJets05, Double_t &nTracks03, Double_t &nTracks05, Double_t &sumPt03, Double_t &sumPt05, Double_t &hadVeto03, Double_t &hadVeto05, Double_t &emVeto03, Double_t &emVeto05, Double_t &trVeto03, Double_t &trVeto05, Double_t &EnMu1, Double_t &EnMu2, Double_t &EnMu3, Double_t &ChargeMu1, Double_t &ChargeMu2, Double_t &ChargeMu3, Double_t &isQValid1, Double_t &isTValid1, Double_t &isIsoValid1, Double_t &GLnormChi2_mu1, Double_t &GL_nValidMuHits1, Double_t &trkLayersWMeas1, Double_t &nValidPixelHits1, Double_t &outerTrk_P_1, Double_t &outerTrk_Eta_1, Double_t &outerTrk_normChi2_1, Double_t &outerTrk_muStValidHits_1, Double_t &innerTrk_P_1, Double_t &innerTrk_Eta_1, Double_t &innerTrk_normChi2_1, Double_t &QInnerOuter_1, Double_t &cQ_uS_1, Double_t &cQ_tK_1, Double_t &cQ_gK_1, Double_t &cQ_tRChi2_1, Double_t &cQ_sRChi2_1, Double_t &cQ_Chi2LP_1, Double_t &cQ_Chi2LM_1, Double_t &cQ_lD_1, Double_t &cQ_gDEP_1, Double_t &cQ_tM_1, Double_t &cQ_gTP_1, Double_t &segmComp_1, Double_t &caloComp_1, Double_t &isQValid2, Double_t &isTValid2, Double_t &isIsoValid2, Double_t &GLnormChi2_mu2, Double_t &GL_nValidMuHits2, Double_t &trkLayersWMeas2, Double_t &nValidPixelHits2, Double_t &outerTrk_P_2, Double_t &outerTrk_Eta_2, Double_t &outerTrk_normChi2_2, Double_t &outerTrk_muStValidHits_2, Double_t &innerTrk_P_2, Double_t &innerTrk_Eta_2, Double_t &innerTrk_normChi2_2, Double_t &QInnerOuter_2, Double_t &cQ_uS_2, Double_t &cQ_tK_2, Double_t &cQ_gK_2, Double_t &cQ_tRChi2_2, Double_t &cQ_sRChi2_2, Double_t &cQ_Chi2LP_2, Double_t &cQ_Chi2LM_2, Double_t &cQ_lD_2, Double_t &cQ_gDEP_2, Double_t &cQ_tM_2, Double_t &cQ_gTP_2, Double_t &segmComp_2, Double_t &caloComp_2, Double_t &isQValid3, Double_t &isTValid3, Double_t &isIsoValid3, Double_t &GLnormChi2_mu3, Double_t &GL_nValidMuHits3, Double_t &trkLayersWMeas3, Double_t &nValidPixelHits3, Double_t &outerTrk_P_3, Double_t &outerTrk_Eta_3, Double_t &outerTrk_normChi2_3, Double_t &outerTrk_muStValidHits_3, Double_t &innerTrk_P_3, Double_t &innerTrk_Eta_3, Double_t &innerTrk_normChi2_3, Double_t &QInnerOuter_3, Double_t &cQ_uS_3, Double_t &cQ_tK_3, Double_t &cQ_gK_3, Double_t &cQ_tRChi2_3, Double_t &cQ_sRChi2_3, Double_t &cQ_Chi2LP_3, Double_t &cQ_Chi2LM_3, Double_t &cQ_lD_3, Double_t &cQ_gDEP_3, Double_t &cQ_tM_3, Double_t &cQ_gTP_3, Double_t &segmComp_3, Double_t &caloComp_3, Double_t &trk_dZ, Double_t &trk_dXY);
};

#endif

#ifdef ntupleClass_Control_cxx
ntupleClass_Control::ntupleClass_Control(TTree *tree, TString fname) : fChain(0)
{
   Init(tree);
   fileName=fname;
}

ntupleClass_Control::~ntupleClass_Control()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntupleClass_Control::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntupleClass_Control::LoadTree(Long64_t entry)
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

void ntupleClass_Control::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

    // Set object pointer
    evt = 0;
    lumi = 0;
    run = 0;
    NGoodTriplets = 0;
    nPileUpInt = 0;
    GenParticle_PdgId = 0;
    GenParticle_Pt = 0;
    GenParticle_Eta = 0;
    GenParticle_Phi = 0;
    GenParticle_isDs = 0;
    GenParticle_isB = 0;
    GenParticle_isBdecay = 0;
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
    Muon_timeAtIpInOutErr = 0;
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
    Track_pt = 0;
    Track_eta = 0;
    Track_phi = 0;
    Track_normalizedChi2 = 0;
    Track_numberOfValidHits = 0;
    Track_charge = 0;
    Track_dxy = 0;
    Track_dxyError = 0;
    Track_dz = 0;
    Track_dzError = 0;
    Track_vx = 0;
    Track_vy = 0;
    Track_vz = 0;
    Trigger_l1name = 0;
    Trigger_l1decision = 0;
    Trigger_l1prescale = 0;
    Trigger_hltname = 0;
    Trigger_hltdecision = 0;
    Mu01_Pt = 0;
    Mu01_Eta = 0;
    Mu01_Phi = 0;
    Mu01_dRtriggerMatch = 0;
    Mu01_TripletIndex = 0;
    Mu1_dRtriggerMatch_Mu7 = 0;
    Mu1_dRtriggerMatch_Mu8 = 0;
    Mu1_dRtriggerMatch_Mu8_IP5 = 0;
    Mu1_dRtriggerMatch_Mu8_IP6 = 0;
    Mu1_dRtriggerMatch_Mu9_IP0 = 0;
    Mu1_dRtriggerMatch_Mu9_IP3 = 0;
    Mu1_dRtriggerMatch_Mu9_IP4 = 0;
    Mu1_dRtriggerMatch_Mu9_IP5 = 0;
    Mu1_dRtriggerMatch_Mu9_IP6 = 0;
    Mu1_dRtriggerMatch_Mu12_IP6 = 0;
    Mu02_Pt = 0;
    Mu02_Eta = 0;
    Mu02_Phi = 0;
    Mu02_dRtriggerMatch = 0;
    Mu2_dRtriggerMatch_Mu7 = 0;
    Mu2_dRtriggerMatch_Mu8 = 0;
    Mu02_TripletIndex = 0;
    Tr_Pt = 0;
    Tr_Eta = 0;
    Tr_Phi = 0;
    Tr_dRtriggerMatch = 0;
    Tr_TripletIndex = 0;
    selectedTripletsIndex = 0;
    GenMatchMu01_SimPt = 0;
    GenMatchMu02_SimPt = 0;
    GenMatchMu01_SimEta = 0;
    GenMatchMu02_SimEta = 0;
    GenMatchMu01_SimPhi = 0;
    GenMatchMu02_SimPhi = 0;
    GenMatchMu01_Pt = 0;
    GenMatchMu02_Pt = 0;
    GenMatchMu01_Eta = 0;
    GenMatchMu02_Eta = 0;
    GenMatchMu01_Phi = 0;
    GenMatchMu02_Phi = 0;
    Triplet_mindca_iso = 0;
    Triplet_relativeiso = 0;
    TripletVtx2_x = 0;
    TripletVtx2_y = 0;
    TripletVtx2_z = 0;
    TripletVtx2_Chi2 = 0;
    TripletVtx2_NDOF = 0;
    Triplet2_Mass = 0;
    Triplet2_Pt = 0;
    Triplet2_Eta = 0;
    Triplet2_Phi = 0;
    Triplet2_Charge = 0;
    RefittedPV2_x = 0;
    RefittedPV2_y = 0;
    RefittedPV2_z = 0;
    RefittedPV2_NTracks = 0;
    RefittedPV2_isValid = 0;
    FlightDistPVSV2 = 0;
    FlightDistPVSV2_Err = 0;
    FlightDistPVSV2_Significance = 0;
    FlightDistPVSV2_chi2 = 0;
    dxy_mu1 = 0;
    dxy_mu2 = 0;
    dxy_mu3 = 0;
    dxyErr_mu1 = 0;
    dxyErr_mu2 = 0;
    dxyErr_mu3 = 0;
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("evt", &evt, &b_evt);
    fChain->SetBranchAddress("run", &run, &b_run);
    fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
    fChain->SetBranchAddress("NGoodTriplets", &NGoodTriplets, &b_NGoodTriplets);
    fChain->SetBranchAddress("nPileUpInt", &nPileUpInt, &b_nPileUpInt);
    fChain->SetBranchAddress("GenParticle_PdgId", &GenParticle_PdgId, &b_GenParticle_PdgId);
    fChain->SetBranchAddress("GenParticle_Pt", &GenParticle_Pt, &b_GenParticle_Pt);
    fChain->SetBranchAddress("GenParticle_Eta", &GenParticle_Eta, &b_GenParticle_Eta);
    fChain->SetBranchAddress("GenParticle_Phi", &GenParticle_Phi, &b_GenParticle_Phi);
    fChain->SetBranchAddress("GenParticle_isDs", &GenParticle_isDs, &b_GenParticle_isDs);
    fChain->SetBranchAddress("GenParticle_isB", &GenParticle_isB, &b_GenParticle_isB);
    fChain->SetBranchAddress("GenParticle_isBdecay", &GenParticle_isBdecay, &b_GenParticle_isBdecay);
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
    fChain->SetBranchAddress("Muon_timeAtIpInOutErr", &Muon_timeAtIpInOutErr, &b_Muon_timeAtIpInOutErr);
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
    fChain->SetBranchAddress("Track_pt", &Track_pt, &b_Track_pt);
    fChain->SetBranchAddress("Track_eta", &Track_eta, &b_Track_eta);
    fChain->SetBranchAddress("Track_phi", &Track_phi, &b_Track_phi);
    fChain->SetBranchAddress("Track_normalizedChi2", &Track_normalizedChi2, &b_Track_normalizedChi2);
    fChain->SetBranchAddress("Track_numberOfValidHits", &Track_numberOfValidHits, &b_Track_numberOfValidHits);
    fChain->SetBranchAddress("Track_charge", &Track_charge, &b_Track_charge);
    fChain->SetBranchAddress("Track_dxy", &Track_dxy, &b_Track_dxy);
    fChain->SetBranchAddress("Track_dxyError", &Track_dxyError, &b_Track_dxyError);
    fChain->SetBranchAddress("Track_dz", &Track_dz, &b_Track_dz);
    fChain->SetBranchAddress("Track_dzError", &Track_dzError, &b_Track_dzError);
    fChain->SetBranchAddress("Track_vx", &Track_vx, &b_Track_vx);
    fChain->SetBranchAddress("Track_vy", &Track_vy, &b_Track_vy);
    fChain->SetBranchAddress("Track_vz", &Track_vz, &b_Track_vz);
    fChain->SetBranchAddress("PVCollection_Size", &PVCollection_Size, &b_PVCollection_Size);
    fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
    fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
    fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
    fChain->SetBranchAddress("PV_NTracks", &PV_NTracks, &b_PV_NTracks);
    fChain->SetBranchAddress("Trigger_l1name", &Trigger_l1name, &b_Trigger_l1name);
    fChain->SetBranchAddress("Trigger_l1decision", &Trigger_l1decision, &b_Trigger_l1decision);
    fChain->SetBranchAddress("Trigger_l1prescale", &Trigger_l1prescale, &b_Trigger_l1prescale);
    fChain->SetBranchAddress("Trigger_hltname", &Trigger_hltname, &b_Trigger_hltname);
    fChain->SetBranchAddress("Trigger_hltdecision", &Trigger_hltdecision, &b_Trigger_hltdecision);
    fChain->SetBranchAddress("TripletCollectionSize2", &TripletCollectionSize2, &b_TripletCollectionSize2);
    fChain->SetBranchAddress("SelectedTripletsSize", &SelectedTripletsSize, &b_SelectedTripletsSize);
    fChain->SetBranchAddress("Mu01_Pt", &Mu01_Pt, &b_Mu01_Pt);
    fChain->SetBranchAddress("Mu01_Eta", &Mu01_Eta, &b_Mu01_Eta);
    fChain->SetBranchAddress("Mu01_Phi", &Mu01_Phi, &b_Mu01_Phi);
    fChain->SetBranchAddress("Mu01_dRtriggerMatch", &Mu01_dRtriggerMatch, &b_Mu01_dRtriggerMatch);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu7", &Mu1_dRtriggerMatch_Mu7, &b_Mu1_dRtriggerMatch_Mu7);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu8", &Mu1_dRtriggerMatch_Mu8, &b_Mu1_dRtriggerMatch_Mu8);
    fChain->SetBranchAddress("Mu01_TripletIndex", &Mu01_TripletIndex, &b_Mu01_TripletIndex);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu8_IP5", &Mu1_dRtriggerMatch_Mu8_IP5, &b_Mu1_dRtriggerMatch_Mu8_IP5);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu8_IP6", &Mu1_dRtriggerMatch_Mu8_IP6, &b_Mu1_dRtriggerMatch_Mu8_IP6);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu9_IP0", &Mu1_dRtriggerMatch_Mu9_IP0, &b_Mu1_dRtriggerMatch_Mu9_IP0);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu9_IP3", &Mu1_dRtriggerMatch_Mu9_IP3, &b_Mu1_dRtriggerMatch_Mu9_IP3);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu9_IP4", &Mu1_dRtriggerMatch_Mu9_IP4, &b_Mu1_dRtriggerMatch_Mu9_IP4);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu9_IP5", &Mu1_dRtriggerMatch_Mu9_IP5, &b_Mu1_dRtriggerMatch_Mu9_IP5);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu9_IP6", &Mu1_dRtriggerMatch_Mu9_IP6, &b_Mu1_dRtriggerMatch_Mu9_IP6);
    fChain->SetBranchAddress("Mu1_dRtriggerMatch_Mu12_IP6", &Mu1_dRtriggerMatch_Mu12_IP6, &b_Mu1_dRtriggerMatch_Mu12_IP6);
    fChain->SetBranchAddress("Mu02_Pt", &Mu02_Pt, &b_Mu02_Pt);
    fChain->SetBranchAddress("Mu02_Eta", &Mu02_Eta, &b_Mu02_Eta);
    fChain->SetBranchAddress("Mu02_Phi", &Mu02_Phi, &b_Mu02_Phi);
    fChain->SetBranchAddress("Mu02_dRtriggerMatch", &Mu02_dRtriggerMatch, &b_Mu02_dRtriggerMatch);
    fChain->SetBranchAddress("Mu2_dRtriggerMatch_Mu7", &Mu2_dRtriggerMatch_Mu7, &b_Mu2_dRtriggerMatch_Mu7);
    fChain->SetBranchAddress("Mu2_dRtriggerMatch_Mu8", &Mu2_dRtriggerMatch_Mu8, &b_Mu2_dRtriggerMatch_Mu8);
    fChain->SetBranchAddress("Mu02_TripletIndex", &Mu02_TripletIndex, &b_Mu02_TripletIndex);
    fChain->SetBranchAddress("Tr_Pt", &Tr_Pt, &b_Tr_Pt);
    fChain->SetBranchAddress("Tr_Eta", &Tr_Eta, &b_Tr_Eta);
    fChain->SetBranchAddress("Tr_Phi", &Tr_Phi, &b_Tr_Phi);
    fChain->SetBranchAddress("Tr_dRtriggerMatch", &Tr_dRtriggerMatch, &b_Tr_dRtriggerMatch);
    fChain->SetBranchAddress("Tr_TripletIndex", &Tr_TripletIndex, &b_Tr_TripletIndex);
    fChain->SetBranchAddress("selectedTripletsIndex", &selectedTripletsIndex, &b_selectedTripletsIndex);
    fChain->SetBranchAddress("GenMatchMu01_SimPt", &GenMatchMu01_SimPt, &b_GenMatchMu01_SimPt);
    fChain->SetBranchAddress("GenMatchMu02_SimPt", &GenMatchMu02_SimPt, &b_GenMatchMu02_SimPt);
    fChain->SetBranchAddress("GenMatchMu01_SimEta", &GenMatchMu01_SimEta, &b_GenMatchMu01_SimEta);
    fChain->SetBranchAddress("GenMatchMu02_SimEta", &GenMatchMu02_SimEta, &b_GenMatchMu02_SimEta);
    fChain->SetBranchAddress("GenMatchMu01_SimPhi", &GenMatchMu01_SimPhi, &b_GenMatchMu01_SimPhi);
    fChain->SetBranchAddress("GenMatchMu02_SimPhi", &GenMatchMu02_SimPhi, &b_GenMatchMu02_SimPhi);
    fChain->SetBranchAddress("GenMatchMu01_Pt", &GenMatchMu01_Pt, &b_GenMatchMu01_Pt);
    fChain->SetBranchAddress("GenMatchMu02_Pt", &GenMatchMu02_Pt, &b_GenMatchMu02_Pt);
    fChain->SetBranchAddress("GenMatchMu01_Eta", &GenMatchMu01_Eta, &b_GenMatchMu01_Eta);
    fChain->SetBranchAddress("GenMatchMu02_Eta", &GenMatchMu02_Eta, &b_GenMatchMu02_Eta);
    fChain->SetBranchAddress("GenMatchMu01_Phi", &GenMatchMu01_Phi, &b_GenMatchMu01_Phi);
    fChain->SetBranchAddress("GenMatchMu02_Phi", &GenMatchMu02_Phi, &b_GenMatchMu02_Phi);
    fChain->SetBranchAddress("Triplet_mindca_iso", &Triplet_mindca_iso, &b_Triplet_mindca_iso);
    fChain->SetBranchAddress("Triplet_relativeiso", &Triplet_relativeiso, &b_Triplet_relativeiso);
    fChain->SetBranchAddress("TripletVtx2_x", &TripletVtx2_x, &b_TripletVtx2_x);
    fChain->SetBranchAddress("TripletVtx2_y", &TripletVtx2_y, &b_TripletVtx2_y);
    fChain->SetBranchAddress("TripletVtx2_z", &TripletVtx2_z, &b_TripletVtx2_z);
    fChain->SetBranchAddress("TripletVtx2_Chi2", &TripletVtx2_Chi2, &b_TripletVtx2_Chi2);
    fChain->SetBranchAddress("TripletVtx2_NDOF", &TripletVtx2_NDOF, &b_TripletVtx2_NDOF);
    fChain->SetBranchAddress("Triplet2_Mass", &Triplet2_Mass, &b_Triplet2_Mass);
    fChain->SetBranchAddress("Triplet2_Pt", &Triplet2_Pt, &b_Triplet2_Pt);
    fChain->SetBranchAddress("Triplet2_Eta", &Triplet2_Eta, &b_Triplet2_Eta);
    fChain->SetBranchAddress("Triplet2_Phi", &Triplet2_Phi, &b_Triplet2_Phi);
    fChain->SetBranchAddress("Triplet2_Charge", &Triplet2_Charge, &b_Triplet2_Charge);
    fChain->SetBranchAddress("RefittedPV2_x", &RefittedPV2_x, &b_RefittedPV2_x);
    fChain->SetBranchAddress("RefittedPV2_y", &RefittedPV2_y, &b_RefittedPV2_y);
    fChain->SetBranchAddress("RefittedPV2_z", &RefittedPV2_z, &b_RefittedPV2_z);
    fChain->SetBranchAddress("RefittedPV2_NTracks", &RefittedPV2_NTracks, &b_RefittedPV2_NTracks);
    fChain->SetBranchAddress("RefittedPV2_isValid", &RefittedPV2_isValid, &b_RefittedPV2_isValid);
    fChain->SetBranchAddress("FlightDistPVSV2", &FlightDistPVSV2, &b_FlightDistPVSV2);
    fChain->SetBranchAddress("FlightDistPVSV2_Err", &FlightDistPVSV2_Err, &b_FlightDistPVSV2_Err);
    fChain->SetBranchAddress("FlightDistPVSV2_Significance", &FlightDistPVSV2_Significance, &b_FlightDistPVSV2_Significance);
    fChain->SetBranchAddress("FlightDistPVSV2_chi2", &FlightDistPVSV2_chi2, &b_FlightDistPVSV2_chi2);
    fChain->SetBranchAddress("dxy_mu1", &dxy_mu1, &b_dxy_mu1);
    fChain->SetBranchAddress("dxy_mu2", &dxy_mu2, &b_dxy_mu2);
    fChain->SetBranchAddress("dxy_mu3", &dxy_mu3, &b_dxy_mu3);
    fChain->SetBranchAddress("dxyErr_mu1", &dxyErr_mu1, &b_dxyErr_mu1);
    fChain->SetBranchAddress("dxyErr_mu2", &dxyErr_mu2, &b_dxyErr_mu2);
    fChain->SetBranchAddress("dxyErr_mu3", &dxyErr_mu3, &b_dxyErr_mu3);
    Notify();
}

Bool_t ntupleClass_Control::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntupleClass_Control::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntupleClass_Control::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


#endif // #ifdef ntupleClass_Control_cxx
