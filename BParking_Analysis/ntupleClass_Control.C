#define ntupleClass_Control_cxx
#define NCUTS 7
#define NMU_C 2
#define NTOT 3
#define mumass 0.1056583715 // Muon mass in GeV
#define PhiMass 1.019461 // Phi mass in GeV
#define OmegaMass 0.78265 // Omega mass in GeV
#define PiMass 0.1396 // Pi mass in GeV
#define Ds_Mass 1.968 // Ds mass in GeV
#define ptmin 2.0

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include "ntupleClass_Control.h"
#include "Utilities_Control.C"

using namespace std;

void Fill_CutName_Control(TString listCut[NCUTS]);

void ntupleClass_Control::LoopControl(TString type, TString datasetName){
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntries();
    // Variables definition & init
    int mu_Ind[NMU] = {0}, mu[NMU] = {0}, cutevt_onlyParking[NCUTS] = {0};
    bool isMC = false, miniAOD = false;
    // Variables for the final tree
    double run_n = 0, lumi_n = 0, evt_n = 0, nHLT = 0, deltaR_max = 0, deltaZ_max = 0, Pmu3 = 0, cLP = 0, tKink = 0, segmComp = 0, tripletMass = 0, tripletMassReso = 0, fv_nC = 0, fv_dphi3D = 0, fv_d3D = 0, fv_d3Dsig = 0, d0 = 0, d0sig = 0, mindca_iso = 0, trkRel = 0, Pmu1 = 0, Ptmu1 = 0, etamu1 = 0, Pmu2 = 0, Ptmu2 = 0, etamu2 = 0, Ptmu3 = 0, etamu3 = 0, P_trip = 0, Pt_trip = 0, eta_trip = 0, nStationsMu1 = 0, nStationsMu2 = 0, nStationsMu3 = 0, Iso03Mu1 = 0, Iso03Mu2 = 0, Iso03Mu3 = 0, Iso05Mu1 = 0, Iso05Mu2 = 0, Iso05Mu3 = 0, nMatchesMu1 = 0, nMatchesMu2 = 0, nMatchesMu3 = 0, timeAtIpInOutMu_sig1 = 0, timeAtIpInOutMu_sig2 = 0, timeAtIpInOutMu_sig3 = 0, cQ_uS = 0, cQ_tK, cQ_gK = 0, cQ_tRChi2 = 0, cQ_sRChi2 = 0, cQ_Chi2LP = 0, cQ_Chi2LM = 0, cQ_lD = 0, cQ_gDEP = 0, cQ_tM = 0, cQ_gTP = 0, calEn_emMu1 = 0, calEn_emMu2 = 0, calEn_emMu3 = 0, calEn_hadMu1 = 0, calEn_hadMu2 = 0, calEn_hadMu3 = 0, caloComp = 0, fliDistPVSV_Chi2 = 0, isGlb1 = 0, isTracker1 = 0, isLoose1 = 0,  isSoft1 = 0, isPF1 = 0, isRPC1 = 0, isSA1 = 0, isCalo1 = 0, isGlb2 = 0, isTracker2 = 0, isLoose2 = 0,  isSoft2 = 0, isPF2 = 0, isRPC2 = 0, isSA2 = 0, isCalo2 = 0, isGlb3 = 0, isTracker3 = 0, isLoose3 = 0, isSoft3 = 0, isPF3 = 0, isRPC3 = 0, isSA3 = 0, isCalo3 = 0, vx1 = 0, vx2 = 0, vx3 = 0, vy1 = 0, vy2 = 0, vy3 = 0, vz1 = 0, vz2 = 0, vz3 = 0, Refvx1 = 0, Refvx2 = 0, Refvx3 = 0, Refvy1 = 0, Refvy2 = 0, Refvy3 = 0, Refvz1 = 0, Refvz2 = 0, Refvz3 = 0, SVx = 0, SVy = 0, SVz = 0, had03 = 0, had05 = 0, nJets03 = 0, nJets05 = 0, nTracks03 = 0, nTracks05 = 0, sumPt03 = 0, sumPt05 = 0, hadVeto03 = 0, hadVeto05 = 0, emVeto03 = 0, emVeto05 = 0, trVeto03 = 0, trVeto05 = 0, EnMu1 = 0, EnMu2 = 0, EnMu3 = 0, ChargeMu1 = 0, ChargeMu2 = 0, ChargeMu3 = 0, isQValid1 = 0, isTValid1 = 0, isIsoValid1 = 0, GLnormChi2_mu1 = 0, GL_nValidMuHits1 = 0, trkLayersWMeas1 = 0, nValidPixelHits1 = 0, outerTrk_P_1 = 0, outerTrk_Eta_1 = 0, outerTrk_normChi2_1 = 0, outerTrk_muStValidHits_1 = 0, innerTrk_P_1 = 0, innerTrk_Eta_1 = 0, innerTrk_normChi2_1 = 0, QInnerOuter_1 = 0, cQ_uS_1 = 0, cQ_tK_1 = 0, cQ_gK_1 = 0, cQ_tRChi2_1 = 0, cQ_sRChi2_1 = 0, cQ_Chi2LP_1 = 0, cQ_Chi2LM_1 = 0, cQ_lD_1 = 0, cQ_gDEP_1 = 0, cQ_tM_1 = 0, cQ_gTP_1 = 0, segmComp_1 = 0, caloComp_1 = 0, isQValid2 = 0, isTValid2 = 0, isIsoValid2 = 0, GLnormChi2_mu2 = 0, GL_nValidMuHits2 = 0, trkLayersWMeas2 = 0, nValidPixelHits2 = 0, outerTrk_P_2 = 0, outerTrk_Eta_2 = 0, outerTrk_normChi2_2 = 0, outerTrk_muStValidHits_2 = 0, innerTrk_P_2 = 0, innerTrk_Eta_2 = 0, innerTrk_normChi2_2 = 0, QInnerOuter_2 = 0, cQ_uS_2 = 0, cQ_tK_2 = 0, cQ_gK_2 = 0, cQ_tRChi2_2 = 0, cQ_sRChi2_2 = 0, cQ_Chi2LP_2 = 0, cQ_Chi2LM_2 = 0, cQ_lD_2 = 0, cQ_gDEP_2 = 0, cQ_tM_2 = 0, cQ_gTP_2 = 0, segmComp_2 = 0, caloComp_2 = 0, isQValid3 = 0, isTValid3 = 0, isIsoValid3 = 0, GLnormChi2_mu3 = 0, GL_nValidMuHits3 = 0, trkLayersWMeas3 = 0, nValidPixelHits3 = 0, outerTrk_P_3 = 0, outerTrk_Eta_3 = 0, outerTrk_normChi2_3 = 0, outerTrk_muStValidHits_3 = 0, innerTrk_P_3 = 0, innerTrk_Eta_3 = 0, innerTrk_normChi2_3 = 0, QInnerOuter_3 = 0, cQ_uS_3 = 0, cQ_tK_3 = 0, cQ_gK_3 = 0, cQ_tRChi2_3 = 0, cQ_sRChi2_3 = 0, cQ_Chi2LP_3 = 0, cQ_Chi2LM_3 = 0, cQ_lD_3 = 0, cQ_gDEP_3 = 0, cQ_tM_3 = 0, cQ_gTP_3 = 0, segmComp_3 = 0, caloComp_3 = 0, trk_dZ = 0, trk_dXY = 0;
    
    TString listCutEffEvt[NCUTS];
    Fill_CutName_Control(listCutEffEvt);
    if( datasetName.Contains("DsPhiPi") ) isMC = true;
    if( datasetName.Contains("mini") ) miniAOD = true;
    // Creation of the output file
    TString root_fileName = fileName;
    TFile *fout = new TFile(root_fileName, "RECREATE");
    fout->cd();
    TTree *tree = new TTree("FinalTree","FinalTree");
    TreeFin_Init(tree, lumi_n, run_n, evt_n, pileupFactor, nHLT, deltaR_max, deltaZ_max, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu_sig1, timeAtIpInOutMu_sig2, timeAtIpInOutMu_sig3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LP, cQ_Chi2LM, cQ_lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb1, isTracker1, isLoose1,  isSoft1, isPF1, isRPC1, isSA1, isCalo1, isGlb2, isTracker2, isLoose2,  isSoft2, isPF2, isRPC2, isSA2, isCalo2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05, EnMu1, EnMu2, EnMu3, ChargeMu1, ChargeMu2, ChargeMu3, isQValid1, isTValid1, isIsoValid1, GLnormChi2_mu1, GL_nValidMuHits1, trkLayersWMeas1, nValidPixelHits1, outerTrk_P_1, outerTrk_Eta_1, outerTrk_normChi2_1, outerTrk_muStValidHits_1, innerTrk_P_1, innerTrk_Eta_1, innerTrk_normChi2_1, QInnerOuter_1, cQ_uS_1, cQ_tK_1, cQ_gK_1, cQ_tRChi2_1, cQ_sRChi2_1, cQ_Chi2LP_1, cQ_Chi2LM_1, cQ_lD_1, cQ_gDEP_1, cQ_tM_1, cQ_gTP_1, segmComp_1, caloComp_1, isQValid2, isTValid2, isIsoValid2, GLnormChi2_mu2, GL_nValidMuHits2, trkLayersWMeas2, nValidPixelHits2, outerTrk_P_2, outerTrk_Eta_2, outerTrk_normChi2_2, outerTrk_muStValidHits_2, innerTrk_P_2, innerTrk_Eta_2, innerTrk_normChi2_2, QInnerOuter_2, cQ_uS_2, cQ_tK_2, cQ_gK_2, cQ_tRChi2_2, cQ_sRChi2_2, cQ_Chi2LP_2, cQ_Chi2LM_2, cQ_lD_2, cQ_gDEP_2, cQ_tM_2, cQ_gTP_2, segmComp_2, caloComp_2, isQValid3, isTValid3, isIsoValid3, GLnormChi2_mu3, GL_nValidMuHits3, trkLayersWMeas3, nValidPixelHits3, outerTrk_P_3, outerTrk_Eta_3, outerTrk_normChi2_3, outerTrk_muStValidHits_3, innerTrk_P_3, innerTrk_Eta_3, innerTrk_normChi2_3, QInnerOuter_3, cQ_uS_3, cQ_tK_3, cQ_gK_3, cQ_tRChi2_3, cQ_sRChi2_3, cQ_Chi2LP_3, cQ_Chi2LM_3, cQ_lD_3, cQ_gDEP_3, cQ_tM_3, cQ_gTP_3, segmComp_3, caloComp_3, trk_dZ, trk_dXY);
    
    //########### Directory Trigger studies
    TDirectory *dirTrigger_studies = fout->mkdir("Trigger_studies");
    dirTrigger_studies->cd();
    TH1D *hMulti = new TH1D("Multi_evt", "Multi_evt", 20, -0.5, 19.5);
    TH1D *hHLT_deepStudy[NCUTS];
    for(int i=0; i<NCUTS; i++){
        TString h_name = "HLT_cut"; h_name+= i;
        hHLT_deepStudy[i] = new TH1D(h_name, h_name, 8, -1.5, 6.5); }
    fout->cd();
    
    //########### Total histogram
    TH1I *hCutEffEvt_onlyParking = new TH1I("CutEff_onlyParking", "CutEff_onlyParking", NCUTS, 0.5, (NCUTS+0.5));
    
    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        vector<int> ind_GoodTriplets;
//        cout << "Evt n. " << jentry+1 << endl;
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetTree()->GetEntry(ientry);
        //CUT 0 : Before cuts - skip event if no good triplets
        if(NGoodTriplets->at(0) < 1) continue;
        cutevt_onlyParking[0]++;
        //CUT 1 : Check HLT and L1 decision
        int hlt_old = 0, hlt_parking = 0, hlt_Mu_fired[6] = {0}, goodTriplet = 0;
        bool good_muID = false, good_Chi2 = false, good_DimuMass = false; vector<int> ind_tripl;
        double rand = gRandom->Uniform();
        for(int h=0; h<Trigger_hltname->size(); h++) {
            TString hltName = Trigger_hltname->at(h);
            // std HLT
            if( hltName.Contains("HLT_DoubleMu3_TkMu_DsTau3Mu_v") && Trigger_hltdecision->at(h) == 1) { hlt_old = 1; }
            // Parking HLT
            if( (hltName.Contains("HLT_Mu7_IP4") || hltName.Contains("HLT_Mu8_IP") || hltName.Contains("HLT_Mu9_IP") || hltName.Contains("HLT_Mu12_IP")) && Trigger_hltdecision->at(h) == 1){
                if(isMC == true){
                    if(rand <= 0.167 && hltName.Contains("HLT_Mu7_IP4")) { hlt_parking = 1; hlt_Mu_fired[0] = 1; hHLT_deepStudy[1]->Fill(0); break; }
                    if(rand > 0.167 && rand <= 0.186 && hltName.Contains("HLT_Mu8_IP3")) { hlt_parking = 1; hlt_Mu_fired[1] = 1; hHLT_deepStudy[1]->Fill(1); break; }
                    if(rand > 0.186 && rand <= 0.266 && hltName.Contains("HLT_Mu8_IP5")) { hlt_parking = 1; hlt_Mu_fired[2] = 1; hHLT_deepStudy[1]->Fill(2); break; }
                    if(rand > 0.266 && rand <= 0.521 && hltName.Contains("HLT_Mu9_IP5")) { hlt_parking = 1; hlt_Mu_fired[3] = 1; hHLT_deepStudy[1]->Fill(3); break; }
                    if(rand > 0.521 && rand <= 0.807 && hltName.Contains("HLT_Mu9_IP6")) { hlt_parking = 1; hlt_Mu_fired[4] = 1; hHLT_deepStudy[1]->Fill(4); break; }
                    if(rand > 0.807 && rand <= 1 && hltName.Contains("HLT_Mu12_IP6")) { hlt_parking = 1; hlt_Mu_fired[5] = 1; hHLT_deepStudy[1]->Fill(5); break; }
                }
                else{
                    if(hltName.Contains("HLT_Mu7_IP4")) { hlt_parking = 1; hlt_Mu_fired[0] = 1; hHLT_deepStudy[1]->Fill(0); break; }
                    if(hltName.Contains("HLT_Mu8_IP3")) { hlt_parking = 1; hlt_Mu_fired[1] = 1; hHLT_deepStudy[1]->Fill(1); break; }
                    if(hltName.Contains("HLT_Mu8_IP5")) { hlt_parking = 1; hlt_Mu_fired[2] = 1; hHLT_deepStudy[1]->Fill(2); break; }
                    if(hltName.Contains("HLT_Mu9_IP5")) { hlt_parking = 1; hlt_Mu_fired[3] = 1; hHLT_deepStudy[1]->Fill(3); break; }
                    if(hltName.Contains("HLT_Mu9_IP6")) { hlt_parking = 1; hlt_Mu_fired[4] = 1; hHLT_deepStudy[1]->Fill(4); break; }
                    if(hltName.Contains("HLT_Mu12_IP6")) { hlt_parking = 1; hlt_Mu_fired[5] = 1; hHLT_deepStudy[1]->Fill(5); break; }
                }
            }
        }
        
        // Keep only events triggered exclusively by B-Parking triggers
        if(!(hlt_parking == 1 && hlt_old == 0)) continue;
        else cutevt_onlyParking[1]++;
        
        //Loop over the TRIPLETS
        for (int j=0; j<TripletVtx2_Chi2->size(); j++){
            double pt[NMU] = {0}, eta[NMU] = {0}, phi[NMU] = {0};
            MatchIndex("ID", j, mu_Ind, mu);
            Fill_MuonAndTrackVariables(mu_Ind, pt, eta, phi);
            // CUT 2 : 2 Glb & PF mu w/ pt1>7 & pt2>2 + 1 trk w/ pt3>2
            if(Muon_isGlobal->at(mu[0]) == 1 && Muon_isPF->at(mu[0]) == 1 && Muon_isGlobal->at(mu[1]) == 1 && Muon_isPF->at(mu[1]) == 1 && pt[0]>7 && pt[1]>2 && pt[2]>2){
                good_muID = true;
                // CUT 3 : Chi2 of the triplet in [0-15]
                if(TripletVtx2_Chi2->at(j) > 0 && TripletVtx2_Chi2->at(j) < 15){
                    good_Chi2 = true;
                    // CUT 4 : Dimuon mass in [1-1.04] GeV
                    double dimu1_2 = DimuonMass(MuonCharge->at(mu[0]), MuonCharge->at(mu[1]), pt[0], pt[1], eta[0], eta[1], phi[0], phi[1], MuonEnergy->at(mu[0]), MuonEnergy->at(mu[1]));
                    if(dimu1_2 >= 1. && dimu1_2 <= 1.04){
                        good_DimuMass = true;
                        ind_tripl.push_back(j);
                    }
                }
            }
        } // end loop triplets
        
        for(int k=0; k<6; k++){
            if(hlt_Mu_fired[k] == 1){
                if(good_muID == true) hHLT_deepStudy[2]->Fill(k);
                if(good_Chi2 == true) hHLT_deepStudy[3]->Fill(k);
                if(good_DimuMass == true) hHLT_deepStudy[4]->Fill(k);
            }
        }
        if(good_muID == true) cutevt_onlyParking[2]++;
        if(good_Chi2 == true) cutevt_onlyParking[3]++;
        if(good_DimuMass == true) cutevt_onlyParking[4]++;
        if(!(good_muID == true && good_Chi2 == true && good_DimuMass == true)) continue;
        
        //Arbitration
        hMulti->Fill(ind_tripl.size());
        int ind_rif = -999;
        if(ind_tripl.size() > 1) {
            double Chi2_rif = TripletVtx2_Chi2->at(ind_tripl.at(0));
            ind_rif = 0;
            for(int d=1; d<ind_tripl.size(); d++){
                if(TripletVtx2_Chi2->at(ind_tripl.at(d)) < Chi2_rif) {
                    Chi2_rif = TripletVtx2_Chi2->at(ind_tripl.at(d));
                    ind_rif = ind_tripl.at(d); }
            }
        }
        else ind_rif = ind_tripl.at(0);
        // CUT 5 : TripletMass in [1.62-2.20] GeV
        if(Triplet2_Mass->at(ind_rif)<2.20 && Triplet2_Mass->at(ind_rif)>1.62) {
            cutevt_onlyParking[5]++;
            for(int k=0; k<6; k++){ if(hlt_Mu_fired[k] == 1) hHLT_deepStudy[5]->Fill(k); }
        }
        else continue;
        for(int k=0; k<NMU; k++) mu_Ind[k] = ind_rif;
        MatchIndex("ID", ind_rif, mu_Ind, mu);
            
        // CUT 6 : TriggerMatching
        // AOD Trigger Matching
        bool parking_trig = false;
        if(!miniAOD){
            if(hlt_parking == 1 && (
                                    ((Mu1_dRtriggerMatch_Mu7->at(ind_rif)<0.03 || Mu2_dRtriggerMatch_Mu7->at(ind_rif)<0.03) && hlt_Mu_fired[0] == 1) ||
                                    ((Mu1_dRtriggerMatch_Mu8->at(ind_rif)<0.03 || Mu2_dRtriggerMatch_Mu8->at(ind_rif)<0.03) && hlt_Mu_fired[1] == 1) ||
                                    (Mu1_dRtriggerMatch_Mu8_IP5->at(ind_rif)<0.03 && hlt_Mu_fired[2] == 1) ||
                                    (Mu1_dRtriggerMatch_Mu9_IP5->at(ind_rif)<0.03 && hlt_Mu_fired[3] == 1) ||
                                    (Mu1_dRtriggerMatch_Mu9_IP6->at(ind_rif)<0.03 && hlt_Mu_fired[4] == 1) ||
                                    (Mu1_dRtriggerMatch_Mu12_IP6->at(ind_rif)<0.03 && hlt_Mu_fired[5] == 1) )){
                parking_trig = true;
                for(int k=0; k<6; k++){
                    if(hlt_Mu_fired[k] == 1) { hHLT_deepStudy[6]->Fill(k); nHLT = k;}
                }
            }
        }
            
        // Fill tree w/ evt exclusively triggered by Bparking
        if(parking_trig == true && miniAOD == false){
            cutevt_onlyParking[6]++;
            run_n = run; lumi_n = lumi; evt_n = evt;
            TreeFin_Fill(tree, ind_rif, mu_Ind, mu, lumi_n, run_n, evt_n, pileupFactor, nHLT, deltaR_max, deltaZ_max, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu_sig1, timeAtIpInOutMu_sig2, timeAtIpInOutMu_sig3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LP, cQ_Chi2LM, cQ_lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb1, isTracker1, isLoose1,  isSoft1, isPF1, isRPC1, isSA1, isCalo1, isGlb2, isTracker2, isLoose2,  isSoft2, isPF2, isRPC2, isSA2, isCalo2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05, EnMu1, EnMu2, EnMu3, ChargeMu1, ChargeMu2, ChargeMu3, isQValid1, isTValid1, isIsoValid1, GLnormChi2_mu1, GL_nValidMuHits1, trkLayersWMeas1, nValidPixelHits1, outerTrk_P_1, outerTrk_Eta_1, outerTrk_normChi2_1, outerTrk_muStValidHits_1, innerTrk_P_1, innerTrk_Eta_1, innerTrk_normChi2_1, QInnerOuter_1, cQ_uS_1, cQ_tK_1, cQ_gK_1, cQ_tRChi2_1, cQ_sRChi2_1, cQ_Chi2LP_1, cQ_Chi2LM_1, cQ_lD_1, cQ_gDEP_1, cQ_tM_1, cQ_gTP_1, segmComp_1, caloComp_1, isQValid2, isTValid2, isIsoValid2, GLnormChi2_mu2, GL_nValidMuHits2, trkLayersWMeas2, nValidPixelHits2, outerTrk_P_2, outerTrk_Eta_2, outerTrk_normChi2_2, outerTrk_muStValidHits_2, innerTrk_P_2, innerTrk_Eta_2, innerTrk_normChi2_2, QInnerOuter_2, cQ_uS_2, cQ_tK_2, cQ_gK_2, cQ_tRChi2_2, cQ_sRChi2_2, cQ_Chi2LP_2, cQ_Chi2LM_2, cQ_lD_2, cQ_gDEP_2, cQ_tM_2, cQ_gTP_2, segmComp_2, caloComp_2, isQValid3, isTValid3, isIsoValid3, GLnormChi2_mu3, GL_nValidMuHits3, trkLayersWMeas3, nValidPixelHits3, outerTrk_P_3, outerTrk_Eta_3, outerTrk_normChi2_3, outerTrk_muStValidHits_3, innerTrk_P_3, innerTrk_Eta_3, innerTrk_normChi2_3, QInnerOuter_3, cQ_uS_3, cQ_tK_3, cQ_gK_3, cQ_tRChi2_3, cQ_sRChi2_3, cQ_Chi2LP_3, cQ_Chi2LM_3, cQ_lD_3, cQ_gDEP_3, cQ_tM_3, cQ_gTP_3, segmComp_3, caloComp_3, trk_dZ, trk_dXY);
        }
    } // end loop on events
            
    //Histo of cuts Efficiency
    for(int k=0; k<NCUTS; k++){
        hCutEffEvt_onlyParking->Fill(k+1, cutevt_onlyParking[k]);
        hCutEffEvt_onlyParking->GetXaxis()->SetBinLabel(k+1, listCutEffEvt[k]);
    }
    //Write and close the file
    fout->Write();
    fout->Close();
    
}

// ###############################################################

void Fill_CutName_Control(TString listCut[NCUTS]){
    listCut[0] = "BeforeCuts";
    listCut[1] = "HLT";
    listCut[2] = "GlbPFMu";
    listCut[3] = "Chi2[0-15]";
    listCut[4] = "DiMu_Mass";
    listCut[5] = "Tripl_Mass";
    listCut[6] = "TrigMatching";
}
