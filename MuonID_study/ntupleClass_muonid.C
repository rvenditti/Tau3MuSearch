#define ntupleClass_muonid_cxx
#define NCUTS 19
#define NPARTICLES 560
#define NMU 3
#define mumass 0.1056583715
#define PhiMass 1.019461 // Phi mass in GeV
#define OmegaMass 0.78265 // Omega mass in GeV
#define ptmin 2.0
#define catA 0.007
#define catB 0.0105

#include "ntupleClass_muonid.h"
#include "Utilities.C"
#include <stdio.h>
#include <iostream>

int Idsummary2D[NCUTS][NPARTICLES][NPARTICLES] = {0};
int Idsummary2D_Gen[NPARTICLES][NPARTICLES] = {0};


void ntupleClass_muonid::Loop(TString type, TString datasetName){
   
    bool isMC = false;
    if(strcmp(type, "MC") == 0 ) isMC = true;

    bool isVerbose = false;
    
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntries();
    // Variables definition

    // Variables for the final tree
    double run_n, lumi_n, evt_n;
    double puFactor, genP_PdgId, genP_Pt, genP_Eta, genP_Phi, genP_MotherPdgId, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_simPdgId, mu_simMotherPdgId, mu_simFlavour, mu_simType, mu_simBX, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPCMuon, mu_isStandAloneMuon, mu_isTrackerMuon, mu_isCaloMuon, mu_isQualityValid, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05;

    // Creation of output file & final tree
    TString root_fileName = fileName;
    TFile *fout = new TFile(root_fileName, "RECREATE");
    fout->cd();
    TTree *tree_final = new TTree("FinalTree","FinalTree");

    TreeFin_Init(tree_final, run_n, lumi_n, evt_n, puFactor, genP_PdgId, genP_Pt, genP_Eta, genP_Phi, genP_MotherPdgId, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_simPdgId, mu_simMotherPdgId, mu_simFlavour, mu_simType, mu_simBX, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPCMuon, mu_isStandAloneMuon, mu_isTrackerMuon, mu_isCaloMuon, mu_isQualityValid, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05);

    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        //cout << "Event n. " << jentry << endl;
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetTree()->GetEntry(ientry);
        if(!MuonPt->size()>0) continue;

        if(isVerbose) cout<<"=================================\nevt "<<evt<<" run "<<run<<" lumi "<<lumi<<endl;
        run_n = run; evt_n = evt; lumi_n = lumi;

        //Check HLT and L1 decision
        bool hlt_fired = 0;
        for(int h=0; h<Trigger_hltname->size(); h++) {
            TString hltName = Trigger_hltname->at(h);
            //HLT 2018
            if( (datasetName.Contains("2018") != std::string::npos) && strncmp(hltName, "HLT_DoubleMu3_TkMu_DsTau3Mu_v", 29) == 0 && Trigger_hltdecision->at(h) == 1) {
                hlt_fired = 1;
            }
        }

        std::vector<Int_t> genInd_list;
        //Loop on reco muons
        for(int j=0; j<MuonPt->size(); j++){
             if(isVerbose){cout<<"recoMu in event: "<<MuonPt->size()<<endl;}
             Int_t genInd = genParticleMatch(j, isVerbose);
             genInd_list.push_back(genInd); 
             puFactor = 1;       
             //Fill MiniTree
             TreeFin_Fill(tree_final, genInd, j, run_n, lumi_n, evt_n, puFactor, genP_PdgId, genP_Pt, genP_Eta, genP_Phi, genP_MotherPdgId, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_simPdgId, mu_simMotherPdgId, mu_simFlavour, mu_simType, mu_simBX, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPCMuon, mu_isStandAloneMuon, mu_isTrackerMuon, mu_isCaloMuon, mu_isQualityValid, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05);
        }
        
        if (ientry < 0) break;
    }//end loop on events
    cout<<"Done"<<endl;
    //Write and close the file
    fout->Write();
    fout->Close();
}

// #########################################
