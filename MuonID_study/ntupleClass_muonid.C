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
    bool isT3M = false;

    bool doWeights = true;

    if(datasetName.Contains("2018Ds") || datasetName.Contains("2018B0") || datasetName.Contains("2018Bp")) isT3M = true;
 
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntries();
    // Variables definition
    int ind = 0, mu_Ind[NMU] = {0}, mu[NMU] = {0};

    // Variables for the final tree
    double run_n, lumi_n, evt_n;

    double ptetaWeight, genP_PdgId, genP_Pt, genP_Eta, genP_Phi, genP_MotherPdgId, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_simPdgId, mu_simMotherPdgId, mu_simFlavour, mu_simType, mu_simBX, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPC, mu_isStandAlone, mu_isTracker, mu_isCalo, mu_isQualityValid, mu_SoftMVA, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_Numberofvalidtrackerhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_validFraction, mu_innerTrack_highPurity, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_IP3D_BS, mu_IP2D_BS, mu_IP3D_PV, mu_IP2D_PV, mu_validMuonHitComb, mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05;

    // Creation of output file & final tree
    TString root_fileName = fileName;
    TFile *fout = new TFile(root_fileName, "RECREATE");
    fout->cd();
    TTree *tree_final = new TTree("FinalTree","FinalTree");

    TreeFin_Init(tree_final, run_n, lumi_n, evt_n, ptetaWeight, genP_PdgId, genP_Pt, genP_Eta, genP_Phi, genP_MotherPdgId, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_simPdgId, mu_simMotherPdgId, mu_simFlavour, mu_simType, mu_simBX, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPC, mu_isStandAlone, mu_isTracker, mu_isCalo, mu_isQualityValid, mu_SoftMVA, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_Numberofvalidtrackerhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_validFraction, mu_innerTrack_highPurity, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_IP3D_BS, mu_IP2D_BS, mu_IP3D_PV, mu_IP2D_PV, mu_validMuonHitComb, mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05);

    //Open files fot pt-eta reweighting
    TFile *f_weights_barrel = TFile::Open("/lustrehome/fsimone/MuonID_study/Phase_space_reweighing_tools/PT_reweighting_muonid_barrel.root", "read");
    TH2D *h_weights_barrel = (TH2D*)f_weights_barrel->Get("hweight_pteta");
    TFile *f_weights_endcap = TFile::Open("/lustrehome/fsimone/MuonID_study/Phase_space_reweighing_tools/PT_reweighting_muonid_endcap.root", "read");
    TH2D *h_weights_endcap = (TH2D*)f_weights_endcap->Get("hweight_pteta");
    Int_t bin_pt, bin_eta;

    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        //cout << "Event n. " << jentry << endl;
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetTree()->GetEntry(ientry);
        if(!MuonPt->size()>0) continue;

        std::vector<Int_t> triplIndex;

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

        double mu3pt=0; double mu3eta=0; double mu3phi=0;

        if(isT3M){
            //require HLT
            if(!hlt_fired) continue;
            //Loop over the TRIPLETS
            int ntripl = 0;

            if(isVerbose) cout<<"Triplets in the event "<<TripletVtx_Chi2->size()<<endl;
            for (int j=0; j<TripletVtx_Chi2->size(); j++){
                //Matching between index of single mu of the triplet (mu#_Ind) & that of 'MUONID' (mu#)
                MatchIndex("ID", j, mu_Ind, mu);

                if(isVerbose){
                    cout<<"\ntriplet candidate "<<j<<"\n | Pt | Eta | Phi"<<endl;
                    cout<<" | "<<Mu1_Pt->at(j)<<" | "<<Mu1_Eta->at(j)<<" | "<<Mu1_Phi->at(j)<<endl;
                    cout<<" | "<<Mu2_Pt->at(j)<<" | "<<Mu2_Eta->at(j)<<" | "<<Mu2_Phi->at(j)<<endl;
                    cout<<" | "<<Mu3_Pt->at(j)<<" | "<<Mu3_Eta->at(j)<<" | "<<Mu3_Phi->at(j)<<endl;
                    cout<<" | invariant mass: "<<Triplet_Mass->at(j)<<endl;
                    cout<<" | vertex chi2: "<<TripletVtx_Chi2->at(j)<<endl;
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
                if(std::abs(Muon_simPdgId->at(mu[0])) == 13 &&
                   std::abs(Muon_simPdgId->at(mu[1])) == 13 &&
                   std::abs(Muon_simPdgId->at(mu[2])) == 13 &&
                   std::abs(Muon_simMotherPdgId->at(mu[0])) == 15 &&
                   std::abs(Muon_simMotherPdgId->at(mu[1])) == 15 &&
                   std::abs(Muon_simMotherPdgId->at(mu[2])) == 15 ) isSignal = true;
                if(!isSignal) continue; 
                
                ntripl++; triplIndex.push_back(j);
            }//end triplet loop

            if(ntripl==0) continue;
            //Best triplet selected based on smaller Chi2
            ind = BestTripletFinder(triplIndex);
            //RiMatching between index of single muons vs triplet
            MatchIndex("ID", ind, mu_Ind, mu);
            //saving info for mu3
            mu3pt = Mu3_Pt->at(ind);
            mu3eta = Mu3_Eta->at(ind);
            mu3phi = Mu3_Phi->at(ind);
            if(isVerbose) cout<<"mu3pt "<<mu3pt<<" mu3eta "<<mu3eta<<" mu3phi "<<mu3phi<<endl;
        } //isT3M

        std::vector<Int_t> genInd_list;
        //Loop on reco muons
        for(int j=0; j<MuonPt->size(); j++){
             if(isVerbose){cout<<"recoMu in event: "<<MuonPt->size()<<endl;}
             Int_t genInd = genParticleMatch(j, isVerbose);
             genInd_list.push_back(genInd);

             if(!doWeights || isT3M) ptetaWeight = 1;
             else{
                 //compute weights
                 if(abs(MuonEta->at(j))<1.2) {
                     bin_pt =  h_weights_barrel->GetXaxis()->FindBin(MuonPt->at(j)); 
                     bin_eta = h_weights_barrel->GetYaxis()->FindBin(MuonEta->at(j)); 
                     ptetaWeight = h_weights_barrel->GetBinContent(bin_pt, bin_eta);
                 }
                 else {
                     bin_pt =  h_weights_endcap->GetXaxis()->FindBin(MuonPt->at(j)); 
                     bin_eta = h_weights_endcap->GetYaxis()->FindBin(MuonEta->at(j)); 
                     ptetaWeight = h_weights_endcap->GetBinContent(bin_pt, bin_eta);
                 }
             }
             if(isT3M)
                 if(!(MuonPt->at(j)==mu3pt && MuonEta->at(j)==mu3eta && MuonPhi->at(j)==mu3phi)) continue;     

             //Fill MiniTree
             TreeFin_Fill(tree_final, genInd, j, run_n, lumi_n, evt_n, ptetaWeight, genP_PdgId, genP_Pt, genP_Eta, genP_Phi, genP_MotherPdgId, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_simPdgId, mu_simMotherPdgId, mu_simFlavour, mu_simType, mu_simBX, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPC, mu_isStandAlone, mu_isTracker, mu_isCalo, mu_isQualityValid, mu_SoftMVA, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_Numberofvalidtrackerhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_validFraction, mu_innerTrack_highPurity, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_IP3D_BS, mu_IP2D_BS, mu_IP3D_PV, mu_IP2D_PV, mu_validMuonHitComb, mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05);
        }
        
        if (ientry < 0) break;
    }//end loop on events
    cout<<"Done"<<endl;
    //Write and close the file
    fout->Write();
    fout->Close();
}

// #########################################
