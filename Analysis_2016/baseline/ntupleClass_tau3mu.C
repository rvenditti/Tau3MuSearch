#define ntupleClass_tau3mu_cxx
#define NCUTS 18
#define NPARTICLES 560
#define NMU 3
#define mumass 0.1056583715
#define PhiMass 1.019461 // Phi mass in GeV
#define OmegaMass 0.78265 // Omega mass in GeV
#define ptmin 2.0
#define catA 0.007
#define catB 0.0105

#include "ntupleClass_tau3mu.h"
#include "Utilities.C"
#include "PdgId_list.C"
#include <stdio.h>
#include <iostream>

int Idsummary2D[NCUTS][NPARTICLES][NPARTICLES] = {0};
int Idsummary2D_Gen[NPARTICLES][NPARTICLES] = {0};

// ######################################### MC: SIGNAL ANALYSIS CUTFLOW

// Cuts: (over events)
// // * cut[0] -> Before cuts
// // * cut[1] -> Event fires L1 seed Double OR Triple
// // * cut[2] -> Event fires HLT
// Cuts: (over triplets)
// * cut[3] -> not negative Chi2
// * cut[4] -> SV-BS displacement >= 2*std dev
// * cut[5] -> All mu are PF
// * cut[6] -> There are 3 mu glb w/ pt>ptmin (=2) & |eta|<2.4
// * cut[7] -> Triplet mass (in 1.62 - 2 GeV)
// * cut[8] -> The 3 possible pairs of mu of the triplet have proper DeltaR (<0.8)
// * cut[9] -> The 3 possible pairs of mu of the triplet have proper |DeltaZ| (<0.5)
// * cut[10] -> Cut on the dimuon mass w.r.t. Phi(1020) per pairs of mu of the triplet w/ opposite sign
// * cut[11] -> Cut on the dimuon mass w.r.t. Omega(782) per pairs of mu of the triplet w/ opposite sign
// * cut[12] -> Mu1 Trigger Matching dR<0.03, abs(dpt)/pt<0.1
// * cut[13] -> Mu2 Trigger Matching dR<0.03, abs(dpt)/pt<0.1
// * cut[14] -> Mu3 Trigger Matching dR<0.03, abs(dpt)/pt<0.1
// Cuts: (over events containing at least 1 triplet)
// * cut[15] -> Filling HistoStepByStep if HLT fired and L1Double fired
// * cut[16] -> Filling HistoStepByStep if HLT fired and L1Triple fired
// * cut[17] -> Filling HistoStepByStep if HLT and L1Double OR L1Triple fired
//
// N.B.: cut[NCUTS] total number of triplets passing each selection
//       cutevt[NCUTS] total number of events passing each selection
//       cuttripl[NCUTS] number of triplets passing each selection in current event


void ntupleClass_tau3mu::Loop(TString type, TString datasetName){
   
    bool isMC = false;
    if(strcmp(type, "MC") == 0 ) isMC = true;

    bool isVerbose = false;
 
    // Pile-up reweighting
    // data pileup from full 2017
    if( isMC && datasetName.Contains("2017") ) {
        TFile *fPileUp = new TFile("/lustrehome/fsimone/Analysis/Pile_up_reweighing_tools/PileUp_ReweightingStudy_"+datasetName+"Tau3Mu.root");
        TH1 *hPileUpRew = (TH1*)fPileUp->Get("PileUp_Reweighting");
        int Nbins = hPileUpRew->GetNbinsX();
        for(int m=0; m<Nbins; m++){
            pileup_weight.push_back(hPileUpRew->GetBinContent(m));
        }
    }
    // End pile-up reweighting
    
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntries();
    // Variables definition
    int ntripl, trInd = 0, ind = 0, mu_Ind[NMU] = {0}, mu[NMU] = {0}, muGen[NMU] = {0}, NgoodTripl = 0, NbadTripl = 0, cut[NCUTS] = {0}, cutevt[NCUTS] = {0}, cutSignal[NCUTS] = {0}, cutevtSignal[NCUTS] = {0}, Ncut = 0, IdsummaryDaughter[NCUTS][NPARTICLES] = {0}, IdsummaryMother[NCUTS][NPARTICLES] = {0}, IdsummaryDaughter_Gen[NPARTICLES] = {0}, IdsummaryMother_Gen[NPARTICLES] = {0};
    float ptminTrack = 0.5, DeltaRmax = 0.8, DeltaZmax = 0.5;
    double sigmaPhi = 0.011, sigmaOmega = 0.0085, EtaMax = 2.4;
    TString listCut[NCUTS], pId[NPARTICLES];
    std::vector< Double_t > dimu;

    // Variables for the final tree
    double Pmu3 = 0, cLP = 0, segmComp = 0, tripletMass = 0, tripletMassReso = 0, fv_nC = 0, fv_dphi3D = 0, fv_d3D = 0, fv_d3Dsig = 0, d0 = 0, d0sig = 0, mindca_iso = 0, trkRel = 0, Pmu1 = 0, Ptmu1 = 0, etamu1 = 0, Pmu2 = 0, Ptmu2 = 0, etamu2 = 0, Ptmu3 = 0, etamu3 = 0, P_trip = 0, Pt_trip = 0, eta_trip = 0, nStationsMu1 = 0, nStationsMu2 = 0, nStationsMu3 = 0, Iso03Mu1 = 0, Iso03Mu2 = 0, Iso03Mu3 = 0, Iso05Mu1 = 0, Iso05Mu2 = 0, Iso05Mu3 = 0, nMatchesMu1 = 0, nMatchesMu2 = 0, nMatchesMu3 = 0, timeAtIpInOutMu1 = 0, timeAtIpInOutMu2 = 0, timeAtIpInOutMu3 = 0, cQ_uS = 0, cQ_tK, cQ_gK = 0, cQ_tRChi2 = 0, cQ_sRChi2 = 0, cQ_Chi2LM = 0, cQ_Chi2lD = 0, cQ_gDEP = 0, cQ_tM = 0, cQ_gTP = 0, calEn_emMu1 = 0, calEn_emMu2 = 0, calEn_emMu3 = 0, calEn_hadMu1 = 0, calEn_hadMu2 = 0, calEn_hadMu3 = 0, caloComp = 0, fliDistPVSV_Chi2 = 0, isGlb3 = 0, isTracker3 = 0, isLoose3 = 0, isSoft3 = 0, isPF3 = 0, isRPC3 = 0, isSA3 = 0, isCalo3 = 0, Vx1 = 0, Vx2 = 0, Vx3 = 0, Vy1 = 0, Vy2 = 0, Vy3 = 0, Vz1 = 0, Vz2 = 0, Vz3 = 0, Refvx1 = 0, Refvx2 = 0, Refvx3 = 0, Refvy1 = 0, Refvy2 = 0, Refvy3 = 0, Refvz1 = 0, Refvz2 = 0, Refvz3 = 0, SVx = 0, SVy = 0, SVz = 0, had03 = 0, had05 = 0, nJets03 = 0, nJets05 = 0, nTracks03 = 0, nTracks05 = 0, sumPt03 = 0, sumPt05 = 0, hadVeto03 = 0, hadVeto05 = 0, emVeto03 = 0, emVeto05 = 0, trVeto03 = 0, trVeto05 = 0;
    float tKink = 0;
    double run_n = 0, lumi_n = 0, evt_n = 0;
    double dPtriggerMatching = 0;
    //Variables inizialization
    Fill_particleName(pId);
    //cutevt[0] = nentries;
    Fill_CutName(listCut);
    
    gStyle->SetOptStat(1111);
    // Creation of output file & final tree
    TString root_fileName = fileName;
    TFile *fout = new TFile(root_fileName, "RECREATE");
    fout->cd();
    TTree *treeA = new TTree("FinalTreeA_sgn","FinalTreeA_sgn");
    TTree *treeB = new TTree("FinalTreeB_sgn","FinalTreeB_sgn");
    TTree *treeC = new TTree("FinalTreeC_sgn","FinalTreeC_sgn");
    TreeFin_Init(treeA, lumi_n, run_n, evt_n, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, Vx1, Vx2, Vx3, Vy1, Vy2, Vy3, Vz1, Vz2, Vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
    TreeFin_Init(treeB, lumi_n, run_n, evt_n, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, Vx1, Vx2, Vx3, Vy1, Vy2, Vy3, Vz1, Vz2, Vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
    TreeFin_Init(treeC, lumi_n, run_n, evt_n, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso,  fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, Vx1, Vx2, Vx3, Vy1, Vy2, Vy3, Vz1, Vz2, Vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
    // Creation of histograms for variables BEFORE cuts
    TDirectory *dirBeforeCuts = fout->mkdir("BeforeCuts");
    dirBeforeCuts->cd();
    TH1I *hNtrigObj;
    hNtrigObj = new TH1I("hNtrigObj", "hNtrigObj", 100, -0.5, 99.5);
    TH2I *hNtrigObJ_L1;
    hNtrigObJ_L1 = new TH2I("hNtrigObJ_L1", "hNtrigObJ_L1", 100, -0.5, 99.5, 100, -0.5, 99.5);
    TH2I *hNtrigObJ_HLT;
    hNtrigObJ_HLT = new TH2I("hNtrigObJ_HLT", "hNtrigObJ_HLT", 100, -0.5, 99.5, 100, -0.5, 99.5);
    TH1I *hNtrigObj_skimmed;
    hNtrigObj_skimmed = new TH1I("hNtrigObj_skimmed", "hNtrigObj_skimmed", 100, -0.5, 99.5);
    TH1D *hPileUp_BC, *hNPrVert_BC, *hMass_tripl_BC, *hChi2Vertex, *hMass_di_Zero_BC, *hMass_di_Zero2_BC, *hPtRes_BC, *hPtRes_BC_mu[NMU], *hPtResBarrel_BC, *hPtResBarrel_BC_mu[NMU], *hPtResEndcap_BC, *hPtResEndcap_BC_mu[NMU]; TH2D *hMassvsChi2;
    TH1F *hMass_quad_BC, *hMass_quad_Zero_BC;
    TH1I *hPdgId_Gen, *hMotherPdgId_Gen; TH2I *hPdgId2D_Gen;
    hPileUp_BC = new TH1D("hNPileUp", "hNPileUp", 100, -0.5, 99.5);
    hPileUp_BC->Sumw2();
    hNPrVert_BC = new TH1D("hNPrimaryVertices", "hNPrimaryVertices", 100, -0.5, 99.5);
    hNPrVert_BC->Sumw2();
    TCanvas *PdgIdCanvas_Gen, *PdgIdMotherCanvas_Gen, *PdgIdCanvas2D_Gen;
    InitHistoBC(hMass_tripl_BC, hChi2Vertex, hMassvsChi2, hMass_quad_BC, hMass_quad_Zero_BC, hMass_di_Zero_BC, hMass_di_Zero2_BC, hPtRes_BC, hPtRes_BC_mu, hPtResBarrel_BC, hPtResBarrel_BC_mu, hPtResEndcap_BC, hPtResEndcap_BC_mu, hPdgId_Gen, hMotherPdgId_Gen, hPdgId2D_Gen);
    fout->cd();
    // Creation of histograms for variables AFTER cuts
    TDirectory *dirAfterCuts = fout->mkdir("AfterCuts");
    dirAfterCuts->cd();
    TH1I *hNtrigObj_AC;
    hNtrigObj_AC = new TH1I("hNtrigObj_AC", "hNtrigObj_AC", 100, -0.5, 99.5);
    TH2I *hNtrigObJ_AC_L1;
    hNtrigObJ_AC_L1 = new TH2I("hNtrigObJ_AC_L1", "hNtrigObJ_AC_L1", 100, -0.5, 99.5, 100, -0.5, 99.5);
    TH1I *hNtrigObj_AC_skimmed;
    hNtrigObj_AC_skimmed = new TH1I("hNtrigObj_AC_skimmed", "hNtrigObj_AC_skimmed", 100, -0.5, 99.5);
    TH1I *hNtripl; TH1F *hChi2Track, *hmassQuad, *hmassQuad_Zero, *hSegmComp, *hfv_d3Dsig;
    TH1D *hPileUp_AC, *hNPrVert_AC, *hTripTriggerMatched, *hMassTriRes, *hMassTriResBarrel, *hMassTriResEndcap, *hmassdi, *hPtRes_AC, *hPtRes_AC_mu[NMU], *hPtResBarrel_AC, *hPtResBarrel_AC_mu[NMU], *hPtResEndcap_AC, *hPtResEndcap_AC_mu[NMU], *hNMatchedStat, *hFlightDist, *hFlightDist_Signif, *hPtErrOverPt, *hPt_tripl_good, *hPt_tripl_fake, *hDeltaX, *hDeltaY, *hDeltaZ, *hDeltaX_fake, *hDeltaY_fake, *hDeltaZ_fake, *hTripMassA, *hTripMassB, *hTripMassC, *hEtaA, *hEtaB, *hEtaC; TH2D *hFlightDistvsP;
    TH1D *h_FinalInvMass_A1, *h_FinalInvMass_A2, *h_FinalInvMass_B1, *h_FinalInvMass_B2, *h_FinalInvMass_C1, *h_FinalInvMass_C2;
    hPileUp_AC = new TH1D("hNPileUp", "hNPileUp", 100, -0.5, 99.5);
    hPileUp_AC->Sumw2();
    hNPrVert_AC = new TH1D("hNPrimaryVertices", "hNPrimaryVertices", 100, -0.5, 99.5);
    hNPrVert_AC->Sumw2();
    hTripTriggerMatched = new TH1D("hTriplMassTrigMatched", "hTriplMassTrigMatched", 42, 1.60, 2.02); // binning 10 MeV
    hTripTriggerMatched->Sumw2();
    hTripMassA = new TH1D("hTripMassA", "hTripMassA", 42, 1.60, 2.02); // binning 10 MeV
    hTripMassA->Sumw2();
    hTripMassB = new TH1D("hTripMassB", "hTripMassB", 42, 1.60, 2.02); // binning 10 MeV
    hTripMassB->Sumw2();
    hTripMassC = new TH1D("hTripMassC", "hTripMassC", 42, 1.60, 2.02); // binning 10 MeV
    hTripMassC->Sumw2();
    hEtaA = new TH1D("hEtaA", "hEtaA", 50, -0.05, 2.45); // binning 0.05
    hEtaA->Sumw2();
    hEtaB = new TH1D("hEtaB", "hEtaB", 50, -0.05, 2.45); // binning 0.05
    hEtaB->Sumw2();
    hEtaC = new TH1D("hEtaC", "hEtaC", 50, -0.05, 2.45); // binning 0.05
    hEtaC->Sumw2();
    h_FinalInvMass_A1 = new TH1D("h_FinalInvMass_A1", "h_FinalInvMass_A1", 42, 1.60, 2.02); // binning 10 MeV
    h_FinalInvMass_A2 = new TH1D("h_FinalInvMass_A2", "h_FinalInvMass_A2", 42, 1.60, 2.02); // binning 10 MeV
    h_FinalInvMass_B1 = new TH1D("h_FinalInvMass_B1", "h_FinalInvMass_B1", 42, 1.60, 2.02); // binning 10 MeV
    h_FinalInvMass_B2 = new TH1D("h_FinalInvMass_B2", "h_FinalInvMass_B2", 42, 1.60, 2.02); // binning 10 MeV
    h_FinalInvMass_C1 = new TH1D("h_FinalInvMass_C1", "h_FinalInvMass_C1", 42, 1.60, 2.02); // binning 10 MeV
    h_FinalInvMass_C2 = new TH1D("h_FinalInvMass_C2", "h_FinalInvMass_C2", 42, 1.60, 2.02); // binning 10 MeV
    InitHistoAC(hNtripl, hSegmComp, hfv_d3Dsig, hChi2Track, hMassTriRes, hMassTriResBarrel, hMassTriResEndcap, hmassdi, hmassQuad, hmassQuad_Zero, hPtRes_AC, hPtRes_AC_mu, hPtResBarrel_AC, hPtResBarrel_AC_mu, hPtResEndcap_AC, hPtResEndcap_AC_mu, hNMatchedStat, hFlightDist, hFlightDist_Signif, hFlightDistvsP, hPtErrOverPt, hPt_tripl_good, hPt_tripl_fake, hDeltaX, hDeltaY, hDeltaZ, hDeltaX_fake, hDeltaY_fake, hDeltaZ_fake);
    TH1D *hIsolation_03 = new TH1D("Isolation03_AC", "Isolation03_AC", 30, -0.5, 29.5); // binning di 1
    hIsolation_03->Sumw2();
    TH1D *hIsolation_05 = new TH1D("Isolation05_AC", "Isolation05_AC", 30, -0.5, 29.5); // binning di 1
    hIsolation_05->Sumw2();
    fout->cd();
    // Creation of histo StepByStep
    TDirectory *dirStepByStep = fout->mkdir("StepByStep");
    dirStepByStep->cd();
    TDirectory *dirSingleMu = dirStepByStep->mkdir("SingleMu"); // Single mu variables histo
    dirSingleMu->cd();
    TH1D *hPt[NCUTS], *hPt_mu[NCUTS][NMU], *hEta[NCUTS], *hEta_mu[NCUTS][NMU], *hPhi[NCUTS], *hVx[NCUTS], *hVy[NCUTS], *hVz[NCUTS];
    InitHistoStepByStep_SingleMu(hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz);
    dirStepByStep->cd();
    TDirectory *dirPair = dirStepByStep->mkdir("Pair");  // Pair variables histo
    dirPair->cd();
    TH1D *hMass_pair[NCUTS], *hDeltaR_pair[NCUTS], *hDeltaZ_pair[NCUTS];
    InitHistoStepByStep_Pair(hMass_pair, hDeltaR_pair, hDeltaZ_pair);
    dirStepByStep->cd();
    TDirectory *dirTriplet = dirStepByStep->mkdir("Triplet");  // Triplet variables histo
    dirTriplet->cd();
    TH1D *hL1[NCUTS], *hPt_tripl[NCUTS], *hEta_tripl[NCUTS], *hPhi_tripl[NCUTS], *hMass_tripl[NCUTS], *hChi2_tripl[NCUTS];
    InitHistoStepByStep_Triplet(hL1, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl);
    dirStepByStep->cd();
    TDirectory *dirPdgId = dirStepByStep->mkdir("PdgId"); // PdgId histo
    dirPdgId->cd();
    TH1I *hPdgId_cut[NCUTS], *hMotherPdgId_cut[NCUTS]; TH2I *hPdgId2D_cut[NCUTS];
    TCanvas *PdgIdCanvas_cut[NCUTS], *PdgIdMotherCanvas_cut[NCUTS], *PdgIdCanvas2D_cut[NCUTS];
    InitHistoStepByStep_PdgId(hPdgId_cut, hMotherPdgId_cut, hPdgId2D_cut);
    dirStepByStep->cd();
    fout->cd();
    // Creation of total histograms
    TH1I *hCutEff = new TH1I("CutEff_Ntriplets", "CutEff_Ntriplets", NCUTS, 0.5, (NCUTS+0.5));
    TH1I *hCutEffEvt = new TH1I("CutEff_NEvents", "CutEff_NEvents", NCUTS, 0.5, (NCUTS+0.5));
    TH1I *hCutEffSignal = new TH1I("CutEff_Ntriplets_isSignal", "CutEff_Ntriplets_isSignal", NCUTS, 0.5, (NCUTS+0.5));
    TH1I *hCutEffEvtSignal = new TH1I("CutEff_NEvents_isSignal", "CutEff_NEvents_isSignal", NCUTS, 0.5, (NCUTS+0.5));
    
    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        //cout << "Event n. " << jentry << endl;
        std::vector<Int_t> triplIndex;
	ntripl = 0, trInd = 0; int cuttripl[NCUTS] = {0}; int cuttriplSignal[NCUTS] = {0};
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetTree()->GetEntry(ientry);

        //Skip event if no good triplets
        if(NGoodTriplets->at(0) == 0) continue;
        cout<<"=================================\nevt "<<evt<<" run "<<run<<" lumi "<<lumi<<endl;

        
        //Define vector of arrays (Pt, Eta, Phi) for trigger objects in the event 
        std::vector< std::array<double, 3> > Muon_HLT;
        for(std::size_t k=0; k<MuonPt_HLT->size(); k++){
            std::array<double, 3> temp;
            temp[0] = MuonPt_HLT->at(k); temp[1] = MuonEta_HLT->at(k); temp[2] = MuonPhi_HLT->at(k);
            Muon_HLT.push_back(temp);
        }
        if(isVerbose){
            cout<<"\n"<<MuonPt_HLT->size()<<" TriggerObjects in the event:\n | Pt | Eta | Phi "<<endl;
            for( auto const& mu: Muon_HLT ){
                for( auto const& var: mu ){
                     cout<<" | "<<var;
                } cout<<"\n";
            } cout<<"\n";
        }
        //Sort the vector and removes the duplicates
        std::sort(Muon_HLT.begin(), Muon_HLT.end());
        Muon_HLT.erase(std::unique(Muon_HLT.begin(), Muon_HLT.end()), Muon_HLT.end());
        if(isVerbose){
            cout<<"After sorting and removing duplicates "<<endl;
            for( auto const& mu: Muon_HLT ){
                for( auto const& var: mu ){
                     cout<<" | "<<var;
                } cout<<"\n";
            } cout<<"\n";
        }

        //CUT 0 : Before cuts
        Ncut = 0; cutevt[Ncut]++; cut[Ncut] += TripletVtx_Chi2->size();

        //Check HLT and L1 decision
        bool hlt_fired = 0;
        bool l1double_fired = 0;
        bool l1triple_fired = 0;
        int n_l1fired = 0;
        int n_hltfired = 0;
        for(int h=0; h<Trigger_hltname->size(); h++) {
            TString hltName = Trigger_hltname->at(h);
            //HLT 2016
            if( (datasetName.Contains("2016") != std::string::npos) && strncmp(hltName, "HLT_DoubleMu3_Trk_Tau3mu_v", 26) == 0 && Trigger_hltdecision->at(h) == 1) {
                hlt_fired = 1;
                n_hltfired++;
            }
            //HLT 2017
            if( (datasetName.Contains("2017") != std::string::npos) && strncmp(hltName, "HLT_DoubleMu3_Trk_Tau3mu_v", 26) == 0 && Trigger_hltdecision->at(h) == 1) {
                hlt_fired = 1;
                n_hltfired++;
            }
            //HLT 2018
            if( (datasetName.Contains("2018") != std::string::npos) && strncmp(hltName, "HLT_DoubleMu3_TkMu_DsTau3Mu_v", 29) == 0 && Trigger_hltdecision->at(h) == 1) {
                hlt_fired = 1;
                n_hltfired++;
            }
        }
        for(int k=0; k<Trigger_l1name->size(); k++) {
            TString l1Name = Trigger_l1name->at(k);
            //2016 + 2017 + 2018 L1_DoubleMu
            if( ( 
                  strcmp(l1Name, "L1_DoubleMu0er1p6_dEta_Max1p8_OS") == 0 || //2016 Main
                  strcmp(l1Name, "L1_DoubleMu_10_0_dEta_Max1p8") == 0 ||     //2016 Main
                  strcmp(l1Name, "L1_DoubleMu_11_4") == 0 ||                 //2016 Backup
                  strcmp(l1Name, "L1_DoubleMu0er1p4_dEta_Max1p8_OS") == 0 || //2016 Backup
                  strcmp(l1Name, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") == 0 ||//2017 + 2018
                  strcmp(l1Name, "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4") == 0 ||//2017 + 2018 Backup
                  strcmp(l1Name, "L1_DoubleMu4_SQ_OS_dR_Max1p2") == 0 ||     //2017 + 2018
                  strcmp(l1Name, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2") == 0 )    //2017 + 2018 Backup
                  && Trigger_l1decision->at(k) == 1){
                l1double_fired = 1;
                n_l1fired++;
            //    if ( hlt_fired == 1 ) h_l1prescaleDouble->Fill(Trigger_l1prescale->at(k));
            }
            //2016 + 2017 + 2018  L1_TripleMu
            if( ( 
                 strcmp(l1Name, "L1_TripleMu0") == 0 ||                                      //2016 Main
                 strcmp(l1Name, "L1_TripleMu_5_0_0") == 0 ||                                 //2016 Backup
                 strcmp(l1Name, "L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17") == 0 ||      //2017 + 2018
                 strcmp(l1Name, "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") == 0 ) //2017 + 2018
                 && Trigger_l1decision->at(k) == 1){
                l1triple_fired = 1;
                n_l1fired++;
            //    if ( hlt_fired == 1 ) h_l1prescaleTriple->Fill(Trigger_l1prescale->at(k));
            }
        }
        bool isTrigger_forAna = 0;
        if( hlt_fired == 1 && ( l1double_fired == 1 || l1triple_fired == 1 ) ) isTrigger_forAna = 1;
        if (l1double_fired == 1) hL1[0]->Fill(1);
        if (l1triple_fired == 1) hL1[0]->Fill(2);
        if (l1double_fired == 1 || l1triple_fired == 1) hL1[0]->Fill(3);
        if(l1double_fired || l1triple_fired){ 
            if (l1double_fired == 1) hL1[1]->Fill(1);
            if (l1triple_fired == 1) hL1[1]->Fill(2);
            if (l1double_fired == 1 || l1triple_fired == 1) hL1[1]->Fill(3); }
        if(hlt_fired) { 
        if (l1double_fired == 1) hL1[2]->Fill(1);
        if (l1triple_fired == 1) hL1[2]->Fill(2);
        if (l1double_fired == 1 || l1triple_fired == 1) hL1[2]->Fill(3); }

        //Filling number of trigger objects histograms before and after removing duplicates
        std::vector<double> trigPt = *MuonPt_HLT;
        hNtrigObj->Fill(trigPt.size());
        hNtrigObJ_L1->Fill(trigPt.size(), n_l1fired);
        hNtrigObJ_HLT->Fill(trigPt.size(), n_hltfired);
        //sort and remove duplicates
        sort( trigPt.begin(), trigPt.end() );
        trigPt.erase( unique( trigPt.begin(), trigPt.end() ), trigPt.end() );
        hNtrigObj_skimmed->Fill(trigPt.size());

        //cout<<"l1triple_fired "<<l1triple_fired<<" l1double_fired "<<l1double_fired<<endl;
        pileupFactor = 1;
       // if(!(datasetName.Contains("2018"))) {
       //     hPileUp_BC->Fill(nPileUpInt);
       //     //cout << "N. pileUpInt = " << nPileUpInt << endl;
       //     cout<<"nPileUpInt "<<nPileUpInt<<endl;
       //     if(nPileUpInt<80) pileupFactor = pileup_weight.at(nPileUpInt);
       //     cout << "PileUpFactor : " << pileupFactor << endl << endl;
       //     hNPrVert_BC->Fill(PVCollection_Size, pileupFactor);
       // }
        //Loop over the TRIPLETS
        for (int j=0; j<TripletVtx_Chi2->size(); j++){
            //Matching between index of single mu of the triplet (mu#_Ind) & that of 'MUONID' (mu#)
            MatchIndex("ID", j, mu_Ind, mu);

            if(isVerbose){
                cout<<"triplet candidate "<<j<<"\n | Pt | Eta | Phi"<<endl;
                cout<<" | "<<Mu1_Pt->at(j)<<" | "<<Mu1_Eta->at(j)<<" | "<<Mu1_Phi->at(j)<<endl;
                cout<<" | "<<Mu2_Pt->at(j)<<" | "<<Mu2_Eta->at(j)<<" | "<<Mu2_Phi->at(j)<<endl;
                cout<<" | "<<Mu3_Pt->at(j)<<" | "<<Mu3_Eta->at(j)<<" | "<<Mu3_Phi->at(j)<<endl;
                cout<<" | invariant mass: "<<Triplet_Mass->at(j)<<endl;
                cout<<" | muon IDs \n | isPF | isGlobal | isTrackerMuon "<<endl;
                cout<<" | "<<Muon_isPF->at(mu[0])<<" | "<<Muon_isGlobal->at(mu[0])<<" | "<<Muon_isTrackerMuon->at(mu[0])<<endl;
                cout<<" | "<<Muon_isPF->at(mu[1])<<" | "<<Muon_isGlobal->at(mu[1])<<" | "<<Muon_isTrackerMuon->at(mu[1])<<endl;
                cout<<" | "<<Muon_isPF->at(mu[2])<<" | "<<Muon_isGlobal->at(mu[2])<<" | "<<Muon_isTrackerMuon->at(mu[2])<<endl;
                if(isMC && !(datasetName.Contains("2016"))){
                    cout<<" | simInfo:\n | mu1_pdgID | mu2_pdgId | mu3_pdgId "<<endl;
                    cout<<" | "<<Muon_simPdgId->at(mu[0])<<" | "<<Muon_simPdgId->at(mu[0])<<" | "<<Muon_simPdgId->at(mu[0])<<endl;
                    cout<<" | mu1_MotherpdgID | mu2_MotherpdgId | mu3_MotherpdgId "<<endl;
                    cout<<" | "<<Muon_simMotherPdgId->at(mu[0])<<" | "<<Muon_simMotherPdgId->at(mu[0])<<" | "<<Muon_simMotherPdgId->at(mu[0])<<endl;
                }
            }

            bool isSignal = false;
            if(isMC && !(datasetName.Contains("2016"))){
                if(std::abs(Muon_simPdgId->at(mu[0])) == 13 &&
                   std::abs(Muon_simPdgId->at(mu[1])) == 13 &&
                   std::abs(Muon_simPdgId->at(mu[2])) == 13 &&
                   std::abs(Muon_simMotherPdgId->at(mu[0])) == 15 &&
                   std::abs(Muon_simMotherPdgId->at(mu[1])) == 15 &&
                   std::abs(Muon_simMotherPdgId->at(mu[2])) == 15 ) isSignal = true;
            }
            //trigIndex contains the 3 indeces for the trigger objects with minimum deltaR
            std::vector< std::size_t > trigIndex_deltaR = trigMatchDeltaR(j, Muon_HLT, false);

            //trigIndex contains the 3 indeces for the trigger objects with minimum deltaP/P
            std::vector< std::size_t > trigIndex_deltaP = trigMatchDeltaP(j, Muon_HLT, false);

            // Fill histograms before selections
            FillHistoBC("MC", j, hMass_tripl_BC, hChi2Vertex, hMassvsChi2, hMass_quad_BC, hMass_quad_Zero_BC, hMass_di_Zero_BC, hMass_di_Zero2_BC, hPtRes_BC, hPtRes_BC_mu, hPtResBarrel_BC, hPtResBarrel_BC_mu, hPtResEndcap_BC, hPtResEndcap_BC_mu, IdsummaryDaughter_Gen, IdsummaryMother_Gen, Idsummary2D_Gen);

            //CUT 1 : L1 fired
            if(!l1double_fired && !l1triple_fired) continue;
              Ncut=1; cut[Ncut]++; cutevt[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cutevtSignal[Ncut]++;}
              if (isVerbose) cout<<jentry<<" evt "<<evt<<" triplet "<<j<<"   passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 2 : HLT fired
            if(!hlt_fired || (MuonPt_HLT->size()==0) ) continue;//in 2018 ntuples no filter on HLT --> many events have 0 trigger objects
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cutevtSignal[Ncut]++;}
              if (isVerbose) cout<<jentry<<" evt "<<evt<<" triplet "<<j<<"   passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 3: Chi2 > 0
            if( TripletVtx_Chi2->at(j) < 0 ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<" passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 4: SV displaced from the beamspot by at least 2 std dev
            if( FlightDistBS_SV_Significance->at(j) < 2 ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"  passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 5 : all mu are PF
            if( (Muon_isPF->at(mu[0]) == 0) || (Muon_isPF->at(mu[1]) == 0) || (Muon_isPF->at(mu[2]) == 0) ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"   passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 6 // Check that mu1 is glb & pt>2 & |eta|<Etamax
            if( !((Muon_isGlobal->at(mu[0]) == 1) && (MuonPt->at(mu[0]) > ptmin) && (abs(MuonEta->at(mu[0])) < EtaMax) )) continue;
                     // Check that mu2 is glb & pt>2 & |eta|<Etamax
            if( !((Muon_isGlobal->at(mu[1]) == 1) && (MuonPt->at(mu[1]) > ptmin) && (abs(MuonEta->at(mu[1])) < EtaMax) )) continue;
                     // Check that mu3 is glb & pt>2 & |eta|<Etamax
            if( !((Muon_isGlobal->at(mu[2]) == 1) && (MuonPt->at(mu[2]) > ptmin) && (abs(MuonEta->at(mu[2])) < EtaMax) )) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"    passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 7: check condition on trimuon mass
            if( !(Triplet_Mass->at(j) > 1.62 && Triplet_Mass->at(j) < 2) ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"     passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 8: Loop on PAIRS of muons of the triplet & check DeltaR
            if( !(isPairDeltaRGood(j, DeltaRmax)) ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"      passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 9: Check |Delta Z|
            Float_t vz1 = Muon_vz->at(mu[0]);
            Float_t vz2 = Muon_vz->at(mu[1]);
            Float_t vz3 = Muon_vz->at(mu[2]);
            if( !(isPairDeltaZGood(vz1, vz2, vz3, DeltaZmax)) ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"       passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 10: VETO on Phi(1020) mass
            dimu = Compute_DimuonMass(mu_Ind, mu);
            if(isPhi(dimu)) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"        passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 11: VETO on Omega(782) mass
            //if(isOmega(dimu)) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"         passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            Float_t dR1 = dR( Mu1_Eta->at(j), Muon_HLT[trigIndex_deltaR[0]][1], Mu1_Phi->at(j), Muon_HLT[trigIndex_deltaR[0]][2]);
            Float_t dR2 = dR( Mu2_Eta->at(j), Muon_HLT[trigIndex_deltaR[1]][1], Mu2_Phi->at(j), Muon_HLT[trigIndex_deltaR[1]][2]);
            Float_t dR3 = dR( Mu3_Eta->at(j), Muon_HLT[trigIndex_deltaR[2]][1], Mu3_Phi->at(j), Muon_HLT[trigIndex_deltaR[2]][2]);

            Float_t dP1 = std::abs(Muon_HLT[trigIndex_deltaP[0]][0] - MuonPt->at(mu[0]))/MuonPt->at(mu[0]);
            Float_t dP2 = std::abs(Muon_HLT[trigIndex_deltaP[1]][0] - MuonPt->at(mu[1]))/MuonPt->at(mu[1]);
            Float_t dP3 = std::abs(Muon_HLT[trigIndex_deltaP[2]][0] - MuonPt->at(mu[2]))/MuonPt->at(mu[2]);

            // CUT 12: Mu1 Trigger Matching
            if( !(dR1<0.03) ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"          passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 13: Mu2 Trigger Matching
            if( !(dR2<0.03) ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"           passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 14: Mu3 Trigger Matching
            if( !(dR3<0.03) ) continue;
              Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
              if(isSignal) {cutSignal[Ncut]++; cuttriplSignal[Ncut]++;}
              if (isVerbose) cout<<j<<"            passed cut "<<Ncut<<endl;
              FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
 
            ntripl++; triplIndex.push_back(j);

        }// end loop on triplets

        // N. events that passed each selection
        for (int k=1; k<NCUTS; k++){
            if(cuttripl[k] > 0) cutevt[k]++;
            if(cuttriplSignal[k] > 0) cutevtSignal[k]++;
        }
        // Histo N. triplets passed for each event
        hNtripl->Fill(ntripl);
        if(ntripl > 0 && isTrigger_forAna) {
            //Best triplet selected based on smaller Chi2
            ind = BestTripletFinder(triplIndex);
            //RiMatching between index of single mu of the triplet (mu#_Ind) & that of  'MUONID' (mu#) & Ricomputing the 3 possible dimuon masses
            MatchIndex("ID", ind, mu_Ind, mu);

            bool isSB_tripletMass = 0;
            if( (Triplet_Mass->at(ind) >= 1.62 && Triplet_Mass->at(ind) <= 1.75) || (Triplet_Mass->at(ind) >= 1.80 && Triplet_Mass->at(ind) <= 2.00) ) isSB_tripletMass = 1;

            //Filling number of trigger objects histograms before and after removing duplicates
            std::vector<double> trigPt = *MuonPt_HLT;
            hNtrigObj_AC->Fill(trigPt.size());
            hNtrigObJ_AC_L1->Fill(trigPt.size(), n_l1fired);
             //sort and remove duplicates
            sort( trigPt.begin(), trigPt.end() );
            trigPt.erase( unique( trigPt.begin(), trigPt.end() ), trigPt.end() );
            hNtrigObj_AC_skimmed->Fill(trigPt.size());

            // Resolution & final histograms
            //if(isMC) MatchIndex("Gen", ind, mu_Ind, muGen);
            dimu = Compute_DimuonMass(mu_Ind, mu);

            if( isSB_tripletMass ){
               //CUT 15 : final plot sideband
               cutevt[NCUTS-3]++; cut[NCUTS-3]++;
               FillHistoStepByStep(isMC, ind, mu_Ind, mu, NCUTS-3, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
            }
            if( !isSB_tripletMass ){
               //CUT 16 : final plot peak
               cutevt[NCUTS-2]++; cut[NCUTS-2]++;
               FillHistoStepByStep(isMC, ind, mu_Ind, mu, NCUTS-2, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
            }
            //CUT 17 : final plot full invariant mass range
            cutevt[NCUTS-1]++; cut[NCUTS-1]++;
            FillHistoStepByStep(isMC, ind, mu_Ind, mu, NCUTS-1, l1double_fired, l1triple_fired, hL1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hMass_pair, hDeltaR_pair, hDeltaZ_pair, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, hChi2_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            if( isMC && datasetName.Contains("2017") ) FillHistoResoPt_AC(muGen, hPtRes_AC, hPtRes_AC_mu, hPtResBarrel_AC, hPtResBarrel_AC_mu, hPtResEndcap_AC, hPtResEndcap_AC_mu);
            if( isMC && datasetName.Contains("Ds")) FillHistoResoTriplMass(mu_Ind, mu, hMassTriRes, hMassTriResBarrel, hMassTriResEndcap);

            double tripReso = ResoTriplMass(mu_Ind, mu);

            run_n = run; lumi_n = lumi; evt_n = evt;

            if(tripReso < catA){
                hTripMassA->Fill(Triplet_Mass->at(ind), pileupFactor);
                hEtaA->Fill(std::max(std::max(abs(MuonEta->at(mu[0])),abs(MuonEta->at(mu[1]))) , abs(MuonEta->at(mu[2])) ), pileupFactor);                        
                TreeFin_Fill(treeA, ind, mu_Ind, mu, lumi_n, run_n, evt_n, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, Vx1, Vx2, Vx3, Vy1, Vy2, Vy3, Vz1, Vz2, Vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
            }
            else if (tripReso >= catA && tripReso <= catB){
                hTripMassB->Fill(Triplet_Mass->at(ind), pileupFactor);
                hEtaB->Fill(std::max(std::max(abs(MuonEta->at(mu[0])),abs(MuonEta->at(mu[1]))) , abs(MuonEta->at(mu[2])) ), pileupFactor);                        
                TreeFin_Fill(treeB, ind, mu_Ind, mu, lumi_n, run_n, evt_n, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, Vx1, Vx2, Vx3, Vy1, Vy2, Vy3, Vz1, Vz2, Vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
            }
            else if(tripReso > catB){
                hTripMassC->Fill(Triplet_Mass->at(ind), pileupFactor);
                hEtaC->Fill(std::max(std::max(abs(MuonEta->at(mu[0])),abs(MuonEta->at(mu[1]))) , abs(MuonEta->at(mu[2])) ), pileupFactor);
                TreeFin_Fill(treeC, ind, mu_Ind, mu, lumi_n, run_n, evt_n, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, Vx1, Vx2, Vx3, Vy1, Vy2, Vy3, Vz1, Vz2, Vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
            }
           //
            for(int k=0; k<NMU; k++){
                hIsolation_03->Fill(Muon_emEt03->at(mu[k]), pileupFactor);
                hIsolation_05->Fill(Muon_emEt05->at(mu[k]), pileupFactor);
            }
            FillHistoResoTriplMass(mu_Ind, mu, hMassTriRes, hMassTriResBarrel, hMassTriResEndcap);
            FillHistoAC(ind, mu, hSegmComp, hfv_d3Dsig, hChi2Track, hNMatchedStat, hFlightDist, hFlightDist_Signif, hFlightDistvsP, hPtErrOverPt, hmassdi, dimu, hmassQuad, hmassQuad_Zero);
            hPileUp_AC->Fill(nPileUpInt);
            hNPrVert_AC->Fill(PVCollection_Size, pileupFactor);
         
        }
        
        if (ientry < 0) break;
    }//end loop on events
    //Print general info
    cout << endl;
    cout << "TOTAL N. EVENTS -> " << cutevt[0] << endl << endl;
    cout << "TOTAL N. TRIPLETS -> " << cut[0] << endl << endl;
    cout << "Triplets survived: " << cutevt[13] << " || Good: " << NgoodTripl << " , Bad: " << NbadTripl << endl;
    //Histo of cuts Efficiency
    TCanvas *canvEvt = new TCanvas("CutEfficiency_Nevents", "CutEfficiency_Nevents", 0, 0, 1200, 1000);
    Draw_CutEffCanvas(canvEvt, hCutEffEvt, cutevt, listCut);
    TCanvas *canv = new TCanvas("CutEfficiency_Ntriplets", "CutEfficiency_Ntriplets", 0, 0, 1200, 1000);
    Draw_CutEffCanvas(canv, hCutEff, cut, listCut);
    //Histo of cuts Purity
    if(isMC){
        TCanvas *canvEvtSignal = new TCanvas("CutEfficiency_Nevents_isSignal", "CutEfficiency_Nevents_isSignal", 0, 0, 1200, 1000);
        Draw_CutEffCanvas(canvEvtSignal, hCutEffEvtSignal, cutevtSignal, listCut);
        TCanvas *canvSignal = new TCanvas("CutEfficiency_Ntriplets_isSignal", "CutEfficiency_Ntriplets_isSignal", 0, 0, 1200, 1000);
        Draw_CutEffCanvas(canvSignal, hCutEffSignal, cutSignal, listCut);
    }
    //PdgId histos
    dirStepByStep->cd(); dirPdgId->cd();
    Draw_PdgIdCanvas_StepByStep(PdgIdCanvas_cut, hPdgId_cut, IdsummaryDaughter, PdgIdMotherCanvas_cut, hMotherPdgId_cut, IdsummaryMother, PdgIdCanvas2D_cut, hPdgId2D_cut, Idsummary2D, pId);
    dirStepByStep->cd(); fout->cd();
    //PdgId histos Gen BC
    dirBeforeCuts->cd();
    Draw_PdgIdCanvasGen(PdgIdCanvas_Gen, hPdgId_Gen, IdsummaryDaughter_Gen, PdgIdMotherCanvas_Gen, hMotherPdgId_Gen, IdsummaryMother_Gen, PdgIdCanvas2D_Gen, hPdgId2D_Gen, Idsummary2D_Gen, pId);
    fout->cd();
    //Write and close the file
    fout->Write();
    fout->Close();
}

// #########################################
