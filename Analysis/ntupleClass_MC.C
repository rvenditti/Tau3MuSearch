#define ntupleClass_MC_cxx
#include "ntupleClass_MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



// #########################################  MC DsTau3mu ANALYSIS NEW CUTFLOW

// Cuts: (over triplets)
// * cut[1] -> Chi2 triplet vertex (in 0 - 15)
// * cut[2] -> There are (2 mu glb w/ pt>ptmin (=2) e 1 mu Tracker con pt>ptminTrack (=0.5)) & |eta|<2.4
// * cut[3] -> Triplet mass (in 1.62 - 2 GeV)
// * cut[4] -> The 3 possible pairs of mu of the triplet have proper DeltaR (<0.8)
// * cut[5] -> The 3 possible pairs of mu of the triplet have proper |DeltaZ| (<0.5)
// * cut[6] -> Cut on the dimuon mass w.r.t. Phi(1020) per pairs of mu of the triplet w/ opposite sign
// * cut[7] -> Cut on the dimuon mass w.r.t. Omega(782) per pairs of mu of the triplet w/ opposite sign

void ntupleClass_MC::LoopMC_New(){
    
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    // Variables definition
    int ntripl, trInd = 0, ind = 0, NgoodTripl = 0, NbadTripl = 0, cut[8] = {0}, cutevt[8] = {0}, triplIndex[1000] = {0};
    int IdsummaryDaughter_BC[38] = {0}, IdsummaryMother_BC[38] = {0}, IdsummaryDaughter_AC[38] = {0}, IdsummaryMother_AC[38] = {0};
    float ptmin = 2.0, ptminTrack = 0.5, DeltaRmax = 0.8, DeltaZmax = 0.5, DeltaZ1 = 0, DeltaZ2 = 0, DeltaZ3 = 0, ntripltot = 0;
    double massmin = 1.62, massmax = 2.00, sigmaPhi = 0.011, sigmaOmega = 0.0085, mumass = 0.1056583715; // Mass values in GeV
    double TripletVtx_Chi2max = 15, EtaMax = 2.4;
    TString pId[38], listCut[8];
    //Variables inizialization
    cutevt[0] = nentries;
    particleName(pId);
    cutName_New(listCut);
    // Creation of the file with analysed data
    TFile *fout = new TFile("AnalysedTree_MC_DsTau3mu.root", "RECREATE");
    fout->cd();
    // Creation of histograms for variables BEFORE cuts
    TDirectory *dirBeforeCuts = fout->mkdir("BeforeCuts");
    dirBeforeCuts->cd();
    TH1D *hPtcomp_mu1_BC = new TH1D("PtResolution_Mu1_BC", "PtResolution_Mu1_BC", 160, -2., 2.);
    TH1D *hPtcomp_mu2_BC = new TH1D("PtResolution_Mu2_BC", "PtResolution_Mu2_BC", 160, -2., 2.);
    TH1D *hPtcomp_mu3_BC = new TH1D("PtResolution_Mu3_BC", "PtResolution_Mu3_BC", 160, -2., 2.);
    TH1D *hPtcomp_BC = new TH1D("PtResolution_BC", "PtResolution_BC", 160, -2., 2.);
    TH1D *hPtcompBarrel_mu1_BC = new TH1D("PtResolutionBarrel_Mu1_BC", "PtResolutionBarrel_Mu1_BC", 160, -2., 2.);
    TH1D *hPtcompBarrel_mu2_BC = new TH1D("PtResolutionBarrel_Mu2_BC", "PtResolutionBarrel_Mu2_BC", 160, -2., 2.);
    TH1D *hPtcompBarrel_mu3_BC = new TH1D("PtResolutionBarrel_Mu3_BC", "PtResolutionBarrel_Mu3_BC", 160, -2., 2.);
    TH1D *hPtcompBarrel_BC = new TH1D("PtResolutionBarrel_BC", "PtResolutionBarrel_BC", 160, -2., 2.);
    TH1D *hPtcompEndcap_mu1_BC = new TH1D("PtResolutionEndcap_Mu1_BC", "PtResolutionEndcap_Mu1_BC", 160, -2., 2.);
    TH1D *hPtcompEndcap_mu2_BC = new TH1D("PtResolutionEndcap_Mu2_BC", "PtResolutionEndcap_Mu2_BC", 160, -2., 2.);
    TH1D *hPtcompEndcap_mu3_BC = new TH1D("PtResolutionEndcap_Mu3_BC", "PtResolutionEndcap_Mu3_BC", 160, -2., 2.);
    TH1D *hPtcompEndcap_BC = new TH1D("PtResolutionEndcap_BC", "PtResolutionEndcap_BC", 160, -2., 2.);
    TH1D *hmasstriBeforeCuts = new TH1D("TripletMass_BC", "TripletMass_BC", 400, -0.05, 19.95); // binning 50 MeV
    TH1D *hChi2Vertex = new TH1D("Chi2Vertex_BC", "Chi2Vertex_BC", 40, -1.5, 18.5); // binning 0.5
    TH2D *hMassvsChi2 = new TH2D("MassvsChi2_BC", "MassvsChi2_BC", 200, 1., 2.5, 200, -5., 20.);
    TH1F *hmassQuad_Other = new TH1F("QuadMuonMass_BC", "QuadMuonMass_BC", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Other_Zero = new TH1F("QuadMuonMass_Zero_BC", "QuadMuonMass_Zero_BC", 400, -0.05, 79.95); // binning 200 MeV
    TH1D *hmassDi_Other_Zero = new TH1D("DiMuon_2glbMu_Zero_BC", "DiMuon_2glbMu_Zero_BC", 500, -0.05, 24.95); // binning 50 MeV
    TH1D *hmassDi_Other2_Zero = new TH1D("DiMuon_Other", "DiMuon_Other", 500, -0.05, 24.95); // binning 50 MeV
    TH1I *hPdgId_BC = new TH1I("DaughterPdgId_BC", "DaughterPdgId_BC", 38, -0.5, 37.5);
    TH1I *hMotherPdgId_BC = new TH1I("MotherPdgId_BC", "MotherPdgId_BC", 38, -0.5, 37.5);
    hPtcomp_mu1_BC->Sumw2();
    hPtcomp_mu2_BC->Sumw2();
    hPtcomp_mu3_BC->Sumw2();
    hPtcomp_BC->Sumw2();
    hPtcompBarrel_mu1_BC->Sumw2();
    hPtcompBarrel_mu2_BC->Sumw2();
    hPtcompBarrel_mu3_BC->Sumw2();
    hPtcompBarrel_BC->Sumw2();
    hPtcompEndcap_mu1_BC->Sumw2();
    hPtcompEndcap_mu2_BC->Sumw2();
    hPtcompEndcap_mu3_BC->Sumw2();
    hPtcompEndcap_BC->Sumw2();
    hmasstriBeforeCuts->Sumw2();
    hChi2Vertex->Sumw2();
    hmassQuad_Other->Sumw2();
    hmassQuad_Other_Zero->Sumw2();
    hmassDi_Other_Zero->Sumw2();
    hmassDi_Other2_Zero->Sumw2();
    fout->cd();
    // Creation of INTERMEDIATE histograms
    TDirectory *dirIntermediate = fout->mkdir("Intermediate");
    dirIntermediate->cd();
    TH1D *hmassdi_int = new TH1D("DimuonMass_int", "DimuonMass_int", 500, -0.05, 24.95); // binning 50 MeV
    TH1F *hmassQuad_int = new TH1F("QuadMuonMass_int", "QuadMuonMass_int", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Zero_int = new TH1F("QuadMuonMass_Zero_int", "QuadMuonMass_Zero_int", 400, -0.05, 79.95); // binning 200 MeV
    hmassdi_int->Sumw2();
    hmassQuad_int->Sumw2();
    hmassQuad_Zero_int->Sumw2();
    fout->cd();
    // Creation of histograms for variables AFTER cuts
    TDirectory *dirAfterCuts = fout->mkdir("AfterCuts");
    dirAfterCuts->cd();
    TH1I *hNtripl = new TH1I("Ntripl_AC", "Ntripl_AC", 5, -0.5, 4.5);
    TH1D *hPt = new TH1D("Pt_AC", "Pt_AC", 100, -0.05, 24.95); // binning 250 MeV
    TH1D *hPt_tripl = new TH1D("Pt_tripl_AC", "Pt_tripl_AC", 160, -0.05, 39.95); // binning 250 MeV
    TH1D *hPt_mu1 = new TH1D("Pt_mu1_AC", "Pt_mu1_AC", 100, -0.05, 24.95); // binning 250 MeV
    TH1D *hPt_mu2 = new TH1D("Pt_mu2_AC", "Pt_mu2_AC", 100, -0.05, 24.95); // binning 250 MeV
    TH1D *hPt_mu3 = new TH1D("Pt_mu3_AC", "Pt_mu3_AC", 100, -0.05, 24.95); // binning 250 MeV
    TH1D *hPtcomp_mu1 = new TH1D("PtResolution_Mu1_AC", "PtResolution_Mu1_AC", 160, -2., 2.);
    TH1D *hPtcomp_mu2 = new TH1D("PtResolution_Mu2_AC", "PtResolution_Mu2_AC", 160, -2., 2.);
    TH1D *hPtcomp_mu3 = new TH1D("PtResolution_Mu3_AC", "PtResolution_Mu3_AC", 160, -2., 2.);
    TH1D *hPtcomp = new TH1D("PtResolution_AC", "PtResolution_AC", 160, -2., 2.);
    TH1D *hPtcompBarrel_mu1 = new TH1D("PtResolutionBarrel_Mu1_AC", "PtResolutionBarrel_Mu1_AC", 160, -2., 2.);
    TH1D *hPtcompBarrel_mu2 = new TH1D("PtResolutionBarrel_Mu2_AC", "PtResolutionBarrel_Mu2_AC", 160, -2., 2.);
    TH1D *hPtcompBarrel_mu3 = new TH1D("PtResolutionBarrel_Mu3_AC", "PtResolutionBarrel_Mu3_AC", 160, -2., 2.);
    TH1D *hPtcompBarrel = new TH1D("PtResolutionBarrel_AC", "PtResolutionBarrel_AC", 160, -2., 2.);
    TH1D *hPtcompEndcap_mu1 = new TH1D("PtResolutionEndcap_Mu1_AC", "PtResolutionEndcap_Mu1_AC", 160, -2., 2.);
    TH1D *hPtcompEndcap_mu2 = new TH1D("PtResolutionEndcap_Mu2_AC", "PtResolutionEndcap_Mu2_AC", 160, -2., 2.);
    TH1D *hPtcompEndcap_mu3 = new TH1D("PtResolutionEndcap_Mu3_AC", "PtResolutionEndcap_Mu3_AC", 160, -2., 2.);
    TH1D *hPtcompEndcap = new TH1D("PtResolutionEndcap_AC", "PtResolutionEndcap_AC", 160, -2., 2.);
    TH1D *hEta = new TH1D("Eta_AC", "Eta_AC", 100, -2.5, 2.5); // binning 0.05
    TH1D *hEta_tripl = new TH1D("Eta_tripl_AC", "Eta_tripl_AC", 100, -2.5, 2.5); // binning 0.05
    TH1D *hEta_mu1 = new TH1D("Eta_mu1_AC", "Eta_mu1_AC", 100, -2.5, 2.5); // binning 0.05
    TH1D *hEta_mu2 = new TH1D("Eta_mu2_AC", "Eta_mu2_AC", 100, -2.5, 2.5); // binning 0.05
    TH1D *hEta_mu3 = new TH1D("Eta_mu3_AC", "Eta_mu3_AC", 100, -2.5, 2.5); // binning 0.05
    TH1D *hPhi = new TH1D("Phi_AC", "Phi_AC", 140, -3.5, 3.5); // binning 0.05
    TH1D *hPhi_tripl = new TH1D("Phi_tripl_AC", "Phi_tripl_AC", 140, -3.5, 3.5); // binning 0.05
    TH1F *hChi2Track = new TH1F("Chi2Track", "Chi2Track", 27, -0.3, 5.1); //binning di 0.2
    TH1D *hmasstri = new TH1D("TripletMass_AC", "TripletMass_AC", 42, 1.60, 2.02); // binning 10 MeV
    TH1D *hMassTriResolution = new TH1D("ResolutionTripletMass", "ResolutionTripletMass", 12000, -2., 2.);
    TH1D *hMassTriResolutionBarrel = new TH1D("ResolutionTripletMass_Barrel", "ResolutionTripletMass_Barrel", 12000, -2., 2.);
    TH1D *hMassTriResolutionEndcap = new TH1D("ResolutionTripletMass_Endcap", "ResolutionTripletMass_Endcap", 12000, -2., 2.);
    TH1D *hmassdi = new TH1D("DimuonMass_AC", "DimuonMass_AC", 60, -0.05, 2.95); //binning 50 MeV
    TH1F *hmassQuad = new TH1F("QuadMuonMass_AC", "QuadMuonMass_AC", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Zero = new TH1F("QuadMuonMass_Zero", "QuadMuonMass_Zero", 400, -0.05, 79.95); // binning 200 MeV
    TH1D *hNMatchedStat = new TH1D("NofMatchedStations", "NofMatchedStations", 6, -0.5, 5.5);
    TH1D *hFlightDist = new TH1D("FlightDist_AC", "FlightDist_AC", 90, 0., 3.);
    TH1D *hFlightDist_Signif = new TH1D("FlightDist_Signif_AC", "FlightDist_Signif_AC", 100, 0., 100.);
    TH2D *hFlightDistvsP = new TH2D("FlightDistvsP_AC", "FlightDistvsP_AC", 20, 0., 0.3, 20, 0., 45.);
    TH1I *hPdgId_AC = new TH1I("DaughterPdgId_AC", "DaughterPdgId_AC", 38, -0.5, 37.5);
    TH1I *hMotherPdgId_AC = new TH1I("MotherPdgId_AC", "MotherPdgId_AC", 38, -0.5, 37.5);
    hPt->Sumw2();
    hPt_tripl->Sumw2();
    hPt_mu1->Sumw2();
    hPt_mu2->Sumw2();
    hPt_mu3->Sumw2();
    hPtcomp_mu1->Sumw2();
    hPtcomp_mu2->Sumw2();
    hPtcomp_mu3->Sumw2();
    hPtcomp->Sumw2();
    hPtcompBarrel_mu1->Sumw2();
    hPtcompBarrel_mu2->Sumw2();
    hPtcompBarrel_mu3->Sumw2();
    hPtcompBarrel->Sumw2();
    hPtcompEndcap_mu1->Sumw2();
    hPtcompEndcap_mu2->Sumw2();
    hPtcompEndcap_mu3->Sumw2();
    hPtcompEndcap->Sumw2();
    hEta->Sumw2();
    hEta_tripl->Sumw2();
    hEta_mu1->Sumw2();
    hEta_mu2->Sumw2();
    hEta_mu3->Sumw2();
    hPhi->Sumw2();
    hPhi_tripl->Sumw2();
    hChi2Track->Sumw2();
    hmasstri->Sumw2();
    hMassTriResolution->Sumw2();
    hMassTriResolutionBarrel->Sumw2();
    hMassTriResolutionEndcap->Sumw2();
    hmassdi->Sumw2();
    hmassQuad->Sumw2();
    hmassQuad_Zero->Sumw2();
    hNMatchedStat->Sumw2();
    hFlightDist->Sumw2();
    hFlightDist_Signif->Sumw2();
    fout->cd();
    // Creation of total histograms
    TH1I *hCutEff = new TH1I("CutEff_Ntriplets", "CutEff_Ntriplets", 8, 0.5, 8.5);
    TH1I *hCutEffEvt = new TH1I("CutEff_NEvents", "CutEff_NEvents", 8, 0.5, 8.5);
    
    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        cout << "Events n." << jentry+1 << endl;
        ntripl = 0, trInd = 0; int cutevt2[8] = {0};
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetEntry(ientry);
        //Check sul numero di tracce nel vertice primario
        if(PV_NTracks > 3){
            //Loop over the TRIPLETS
            for (int j=0; j<TripletVtx_Chi2->size(); j++){
                ntripltot++;
                // BEFORE cuts
                // 4 muons inv mass
                if(MuonPt->size() > 3) QuadMuMass(hmassQuad_Other, hmassQuad_Other_Zero);
                // 2 GLB muons inv mass
                if(MuonPt->size() > 1) DiMuMass(hmassDi_Other_Zero, hmassDi_Other2_Zero, ptmin);
                //Matching tra gli indici dei singoli muoni dei tripletto (mu#_Ind) e quelli 'generali' (mu#)
                int mu1_Ind = Mu1_TripletIndex->at(j);
                int mu2_Ind = Mu2_TripletIndex->at(j);
                int mu3_Ind = Mu3_TripletIndex->at(j);
                if (mu1_Ind != j || mu2_Ind != j || mu3_Ind != j) cout << "Indici diversi ?!?" << endl;
                int mu1 = MuonFinder(Mu1_Pt->at(mu1_Ind), Mu1_Eta->at(mu1_Ind), Mu1_Phi->at(mu1_Ind));
                int mu2 = MuonFinder(Mu2_Pt->at(mu2_Ind), Mu2_Eta->at(mu2_Ind), Mu2_Phi->at(mu2_Ind));
                int mu3 = MuonFinder(Mu3_Pt->at(mu3_Ind), Mu3_Eta->at(mu3_Ind), Mu3_Phi->at(mu3_Ind));
                //PdgId histograms
                particleId(Muon_PdgId->at(mu1), IdsummaryDaughter_BC);
                particleId(Muon_PdgId->at(mu2), IdsummaryDaughter_BC);
                particleId(Muon_PdgId->at(mu3), IdsummaryDaughter_BC);
                particleId(Muon_MotherPdgId->at(mu1), IdsummaryMother_BC);
                particleId(Muon_MotherPdgId->at(mu2), IdsummaryMother_BC);
                particleId(Muon_MotherPdgId->at(mu3), IdsummaryMother_BC);
                //Triplet histograms
                hmasstriBeforeCuts->Fill(Triplet_Mass->at(j));
                hChi2Vertex->Fill(TripletVtx_Chi2->at(j));
                hMassvsChi2->Fill(Triplet_Mass->at(j), TripletVtx_Chi2->at(j));
                //Pt resolution histograms
                for (int i=0; i<GenMatchMu1_SimPt->size(); i++){
                    hPtcomp_mu1_BC->Fill((GenMatchMu1_SimPt->at(i)-GenMatchMu1_Pt->at(i))/GenMatchMu1_SimPt->at(i));
                    hPtcomp_BC->Fill((GenMatchMu1_SimPt->at(i)-GenMatchMu1_Pt->at(i))/GenMatchMu1_SimPt->at(i));
                    if(abs(GenMatchMu1_Eta->at(i)) < 1.4){
                        hPtcompBarrel_mu1_BC->Fill((GenMatchMu1_SimPt->at(i)-GenMatchMu1_Pt->at(i))/GenMatchMu1_SimPt->at(i));
                        hPtcompBarrel_BC->Fill((GenMatchMu1_SimPt->at(i)-GenMatchMu1_Pt->at(i))/GenMatchMu1_SimPt->at(i));
                    }
                    else {
                        hPtcompEndcap_mu1_BC->Fill((GenMatchMu1_SimPt->at(i)-GenMatchMu1_Pt->at(i))/GenMatchMu1_SimPt->at(i));
                        hPtcompEndcap_BC->Fill((GenMatchMu1_SimPt->at(i)-GenMatchMu1_Pt->at(i))/GenMatchMu1_SimPt->at(i));
                    }
                }
                for (int i=0; i<GenMatchMu2_SimPt->size(); i++){
                    hPtcomp_mu2_BC->Fill((GenMatchMu2_SimPt->at(i)-GenMatchMu2_Pt->at(i))/GenMatchMu2_SimPt->at(i));
                    hPtcomp_BC->Fill((GenMatchMu2_SimPt->at(i)-GenMatchMu2_Pt->at(i))/GenMatchMu2_SimPt->at(i));
                    if(abs(GenMatchMu2_Eta->at(i)) < 1.4){
                        hPtcompBarrel_mu2_BC->Fill((GenMatchMu2_SimPt->at(i)-GenMatchMu2_Pt->at(i))/GenMatchMu2_SimPt->at(i));
                        hPtcompBarrel_BC->Fill((GenMatchMu2_SimPt->at(i)-GenMatchMu2_Pt->at(i))/GenMatchMu2_SimPt->at(i));
                    }
                    else {
                        hPtcompEndcap_mu2_BC->Fill((GenMatchMu2_SimPt->at(i)-GenMatchMu2_Pt->at(i))/GenMatchMu2_SimPt->at(i));
                        hPtcompEndcap_BC->Fill((GenMatchMu2_SimPt->at(i)-GenMatchMu2_Pt->at(i))/GenMatchMu2_SimPt->at(i));
                    }
                }
                for (int i=0; i<GenMatchMu3_SimPt->size(); i++){
                    hPtcomp_mu3_BC->Fill((GenMatchMu3_SimPt->at(i)-GenMatchMu3_Pt->at(i))/GenMatchMu3_SimPt->at(i));
                    hPtcomp_BC->Fill((GenMatchMu3_SimPt->at(i)-GenMatchMu3_Pt->at(i))/GenMatchMu3_SimPt->at(i));
                    if(abs(GenMatchMu3_Eta->at(i)) < 1.4){
                        hPtcompBarrel_mu3_BC->Fill((GenMatchMu3_SimPt->at(i)-GenMatchMu3_Pt->at(i))/GenMatchMu3_SimPt->at(i));
                        hPtcompBarrel_BC->Fill((GenMatchMu3_SimPt->at(i)-GenMatchMu3_Pt->at(i))/GenMatchMu3_SimPt->at(i));
                    }
                    else {
                        hPtcompEndcap_mu3_BC->Fill((GenMatchMu3_SimPt->at(i)-GenMatchMu3_Pt->at(i))/GenMatchMu3_SimPt->at(i));
                        hPtcompEndcap_BC->Fill((GenMatchMu3_SimPt->at(i)-GenMatchMu3_Pt->at(i))/GenMatchMu3_SimPt->at(i));
                    }
                }
                //check condition on * Chi2 vertex ( 0 < Chi2 < 15)
                if (TripletVtx_Chi2->at(j) < TripletVtx_Chi2max && TripletVtx_Chi2->at(j) > 0 ){
                    cut[1]++; cutevt2[1]++;
                    // Controllo che il muone 1 sia glb e abbia pt>ptmax & |eta|<Etamax
                    if((Muon_isGlobal->at(mu1) == 1) && (MuonPt->at(mu1) > ptmin) && abs(Mu1_Eta->at(j)) < EtaMax){
                        DeltaZ1 = Muon_vz->at(mu1);
                        // Controllo che il muone 2 sia glb e abbia pt>ptmax & |eta|<Etamax
                        if((Muon_isGlobal->at(mu2) == 1) && (MuonPt->at(mu2) > ptmin) && abs(Mu2_Eta->at(j)) < EtaMax){
                            DeltaZ2 = Muon_vz->at(mu2);
                            //Controllo che il muone 3 sia tracker e abbia pt>0.5 & |eta|<Etamax
                            if((Muon_isTrackerMuon->at(mu3) == 1) && (MuonPt->at(mu3) > ptminTrack) && abs(Mu3_Eta->at(j)) < EtaMax){
                                DeltaZ3 = Muon_vz->at(mu3);
                                cut[2]++; cutevt2[2]++;
                                //Riempio intermediate histo
                                // Calcolo le 3 possibili dimuon mass
                                double dimu1_2 = DimuonMass(MuonCharge->at(mu1), MuonCharge->at(mu2), Mu1_Pt->at(mu1_Ind), Mu2_Pt->at(mu2_Ind), Mu1_Eta->at(mu1_Ind), Mu2_Eta->at(mu2_Ind), Mu1_Phi->at(mu1_Ind), Mu2_Phi->at(mu2_Ind), mumass);
                                double dimu2_3 = DimuonMass(MuonCharge->at(mu2), MuonCharge->at(mu3), Mu2_Pt->at(mu2_Ind), Mu3_Pt->at(mu3_Ind), Mu2_Eta->at(mu2_Ind), Mu3_Eta->at(mu3_Ind), Mu2_Phi->at(mu2_Ind), Mu3_Phi->at(mu3_Ind), mumass);
                                double dimu1_3 = DimuonMass(MuonCharge->at(mu1), MuonCharge->at(mu3), Mu1_Pt->at(mu1_Ind), Mu3_Pt->at(mu3_Ind), Mu1_Eta->at(mu1_Ind), Mu3_Eta->at(mu3_Ind), Mu1_Phi->at(mu1_Ind), Mu3_Phi->at(mu3_Ind), mumass);
                                DiMuonHisto(hmassdi_int, dimu1_2, dimu2_3, dimu1_3);
                                //Aggiungo al tripletto una 4a traccia e ne faccio la massa inv
                                QuadMuMass2(hmassQuad_int, hmassQuad_Zero_int, mu1, mu2, mu3);
                                //check condition on trimuon mass
                                if(Triplet_Mass->at(j) > massmin && Triplet_Mass->at(j) < massmax){
                                    cut[3]++; cutevt2[3]++;
                                    // Loop sulle COPPIE di muoni del tripletto e check su DeltaR e |Delta Z|
                                    if(isPairDeltaRGood(j, DeltaRmax) == true){
                                        cut[4]++; cutevt2[4]++;
                                        if(isPairDeltaZGood(DeltaZ1, DeltaZ2, DeltaZ3, DeltaZmax) == true){
                                            cut[5]++; cutevt2[5]++;
                                            // Controllo se sono compatibili con la Phi(1020)
                                            if(isPairNotAPhi(dimu1_2, dimu2_3, dimu1_3, sigmaPhi) == true){
                                                cut[6]++; cutevt2[6]++;
                                                // Controllo se sono compatibili con la Omega(782)
                                                if(isPairNotAOmega(dimu1_2, dimu2_3, dimu1_3, sigmaOmega) == true){
                                                    ntripl++;
                                                    triplIndex[trInd] = j;
                                                    trInd++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } // end loop sui tripletti
        }
        // Conto numero di eventi passati dopo ogni taglio
        for (int k=1; k<8; k++){
            if(cutevt2[k] > 0) cutevt[k]++;
        }
        // Riempio isto con il numero di tripletti passati per evento
        hNtripl->Fill(ntripl);
        if(ntripl > 0) {
            cutevt[7]++; cut[7]++;
            ind = BestTripletFinder(triplIndex, ntripl);
            
            double TripletP = MuonP(Triplet_Pt->at(ind), Triplet_Eta->at(ind), Triplet_Phi->at(ind), mumass);
            //RiMatching tra gli indici dei singoli muoni dei tripletto (mu#_Ind) e quelli 'generali' (mu#)
            int mu1_Ind = Mu1_TripletIndex->at(ind);
            int mu2_Ind = Mu2_TripletIndex->at(ind);
            int mu3_Ind = Mu3_TripletIndex->at(ind);
            if (mu1_Ind != ind || mu2_Ind != ind || mu3_Ind != ind) cout << "Indici diversi ?!?" << endl;
            int mu1 = MuonFinder(Mu1_Pt->at(mu1_Ind), Mu1_Eta->at(mu1_Ind), Mu1_Phi->at(mu1_Ind));
            int mu2 = MuonFinder(Mu2_Pt->at(mu2_Ind), Mu2_Eta->at(mu2_Ind), Mu2_Phi->at(mu2_Ind));
            int mu3 = MuonFinder(Mu3_Pt->at(mu3_Ind), Mu3_Eta->at(mu3_Ind), Mu3_Phi->at(mu3_Ind));
            // Ricalcolo le 3 possibili dimuon mass
            double dimu1_2 = DimuonMass(MuonCharge->at(mu1), MuonCharge->at(mu2), Mu1_Pt->at(mu1_Ind), Mu2_Pt->at(mu2_Ind), Mu1_Eta->at(mu1_Ind), Mu2_Eta->at(mu2_Ind), Mu1_Phi->at(mu1_Ind), Mu2_Phi->at(mu2_Ind), mumass);
            double dimu2_3 = DimuonMass(MuonCharge->at(mu2), MuonCharge->at(mu3), Mu2_Pt->at(mu2_Ind), Mu3_Pt->at(mu3_Ind), Mu2_Eta->at(mu2_Ind), Mu3_Eta->at(mu3_Ind), Mu2_Phi->at(mu2_Ind), Mu3_Phi->at(mu3_Ind), mumass);
            double dimu1_3 = DimuonMass(MuonCharge->at(mu1), MuonCharge->at(mu3), Mu1_Pt->at(mu1_Ind), Mu3_Pt->at(mu3_Ind), Mu1_Eta->at(mu1_Ind), Mu3_Eta->at(mu3_Ind), Mu1_Phi->at(mu1_Ind), Mu3_Phi->at(mu3_Ind), mumass);
            // Riempio isto con caratteristiche dei 3 muoni che hanno superato tutti i tagli
            hPt->Fill(Mu1_Pt->at(mu1_Ind));
            hPt->Fill(Mu2_Pt->at(mu2_Ind));
            hPt->Fill(Mu3_Pt->at(mu3_Ind));
            hPt_tripl->Fill(Triplet_Pt->at(ind));
            hPt_mu1->Fill(Mu1_Pt->at(mu1_Ind));
            hPt_mu2->Fill(Mu2_Pt->at(mu2_Ind));
            hPt_mu3->Fill(Mu3_Pt->at(mu3_Ind));
            // Pt Resolution histograms
            int muGen1 = MuonFinderGen1(Mu1_Pt->at(mu1_Ind), Mu1_Eta->at(mu1_Ind), Mu1_Phi->at(mu1_Ind));
            int muGen2 = MuonFinderGen2(Mu2_Pt->at(mu2_Ind), Mu2_Eta->at(mu2_Ind), Mu2_Phi->at(mu2_Ind));
            int muGen3 = MuonFinderGen3(Mu3_Pt->at(mu3_Ind), Mu3_Eta->at(mu3_Ind), Mu3_Phi->at(mu3_Ind));
            if(muGen1 != -999){
                hPtcomp_mu1->Fill((GenMatchMu1_SimPt->at(muGen1)-GenMatchMu1_Pt->at(muGen1))/GenMatchMu1_SimPt->at(muGen1));
                hPtcomp->Fill((GenMatchMu1_SimPt->at(muGen1)-GenMatchMu1_Pt->at(muGen1))/GenMatchMu1_SimPt->at(muGen1));
                if(abs(GenMatchMu1_Eta->at(muGen1)) < 1.4){
                    hPtcompBarrel_mu1->Fill((GenMatchMu1_SimPt->at(muGen1)-GenMatchMu1_Pt->at(muGen1))/GenMatchMu1_SimPt->at(muGen1));
                    hPtcompBarrel->Fill((GenMatchMu1_SimPt->at(muGen1)-GenMatchMu1_Pt->at(muGen1))/GenMatchMu1_SimPt->at(muGen1));
                }
                else {
                    hPtcompEndcap_mu1->Fill((GenMatchMu1_SimPt->at(muGen1)-GenMatchMu1_Pt->at(muGen1))/GenMatchMu1_SimPt->at(muGen1));
                    hPtcompEndcap->Fill((GenMatchMu1_SimPt->at(muGen1)-GenMatchMu1_Pt->at(muGen1))/GenMatchMu1_SimPt->at(muGen1));
                }
            }
            if(muGen2 != -999){
                hPtcomp_mu2->Fill((GenMatchMu2_SimPt->at(muGen2)-GenMatchMu2_Pt->at(muGen2))/GenMatchMu2_SimPt->at(muGen2));
                hPtcomp->Fill((GenMatchMu2_SimPt->at(muGen2)-GenMatchMu2_Pt->at(muGen2))/GenMatchMu2_SimPt->at(muGen2));
                if(abs(GenMatchMu2_Eta->at(muGen2)) < 1.4){
                    hPtcompBarrel_mu2->Fill((GenMatchMu2_SimPt->at(muGen2)-GenMatchMu2_Pt->at(muGen2))/GenMatchMu2_SimPt->at(muGen2));
                    hPtcompBarrel->Fill((GenMatchMu2_SimPt->at(muGen2)-GenMatchMu2_Pt->at(muGen2))/GenMatchMu2_SimPt->at(muGen2));
                }
                else {
                    hPtcompEndcap_mu2->Fill((GenMatchMu2_SimPt->at(muGen2)-GenMatchMu2_Pt->at(muGen2))/GenMatchMu2_SimPt->at(muGen2));
                    hPtcompEndcap->Fill((GenMatchMu2_SimPt->at(muGen2)-GenMatchMu2_Pt->at(muGen2))/GenMatchMu2_SimPt->at(muGen2));
                }
            }
            if(muGen3 != -999){
                hPtcomp_mu3->Fill((GenMatchMu3_SimPt->at(muGen3)-GenMatchMu3_Pt->at(muGen3))/GenMatchMu3_SimPt->at(muGen3));
                hPtcomp->Fill((GenMatchMu3_SimPt->at(muGen3)-GenMatchMu3_Pt->at(muGen3))/GenMatchMu3_SimPt->at(muGen3));
                if(abs(GenMatchMu3_Eta->at(muGen3)) < 1.4){
                    hPtcompBarrel_mu3->Fill((GenMatchMu3_SimPt->at(muGen3)-GenMatchMu3_Pt->at(muGen3))/GenMatchMu3_SimPt->at(muGen3));
                    hPtcompBarrel->Fill((GenMatchMu3_SimPt->at(muGen3)-GenMatchMu3_Pt->at(muGen3))/GenMatchMu3_SimPt->at(muGen3));
                }
                else {
                    hPtcompEndcap_mu3->Fill((GenMatchMu3_SimPt->at(muGen3)-GenMatchMu3_Pt->at(muGen3))/GenMatchMu3_SimPt->at(muGen3));
                    hPtcompEndcap->Fill((GenMatchMu3_SimPt->at(muGen3)-GenMatchMu3_Pt->at(muGen3))/GenMatchMu3_SimPt->at(muGen3));
                }
            }
            if(muGen1 != -999 && muGen2 != -999 && muGen3 != -999) {
                NgoodTripl++;
                double massTriplSim = TriMuonMass(GenMatchMu1_SimPt->at(muGen1), GenMatchMu2_SimPt->at(muGen2), GenMatchMu3_SimPt->at(muGen3), GenMatchMu1_SimEta->at(muGen1), GenMatchMu2_SimEta->at(muGen2), GenMatchMu3_SimEta->at(muGen3), GenMatchMu1_SimPhi->at(muGen1), GenMatchMu2_SimPhi->at(muGen2), GenMatchMu3_SimPhi->at(muGen3), mumass);
                hMassTriResolution->Fill((massTriplSim-Triplet_Mass->at(ind))/massTriplSim);
                if(abs(GenMatchMu1_Eta->at(muGen1)) < 1.4 && abs(GenMatchMu2_Eta->at(muGen2)) < 1.4 && abs(GenMatchMu3_Eta->at(muGen3)) < 1.4) hMassTriResolutionBarrel->Fill((massTriplSim-Triplet_Mass->at(ind))/massTriplSim);
                if(abs(GenMatchMu1_Eta->at(muGen1)) >= 1.4 && abs(GenMatchMu2_Eta->at(muGen2)) >= 1.4 && abs(GenMatchMu3_Eta->at(muGen3)) >= 1.4) hMassTriResolutionEndcap->Fill((massTriplSim-Triplet_Mass->at(ind))/massTriplSim);
            }
            else NbadTripl++;
            
            hEta->Fill(Mu1_Eta->at(mu1_Ind));
            hEta->Fill(Mu2_Eta->at(mu2_Ind));
            hEta->Fill(Mu3_Eta->at(mu3_Ind));
            hEta_tripl->Fill(Triplet_Eta->at(ind));
            hEta_mu1->Fill(Mu1_Eta->at(mu1_Ind));
            hEta_mu2->Fill(Mu2_Eta->at(mu2_Ind));
            hEta_mu3->Fill(Mu3_Eta->at(mu3_Ind));
            hPhi->Fill(Mu1_Phi->at(mu1_Ind));
            hPhi->Fill(Mu2_Phi->at(mu2_Ind));
            hPhi->Fill(Mu3_Phi->at(mu3_Ind));
            hPhi_tripl->Fill(Triplet_Phi->at(ind));
            hChi2Track->Fill(Muon_innerTrack_normalizedChi2->at(mu1));
            hChi2Track->Fill(Muon_innerTrack_normalizedChi2->at(mu2));
            hChi2Track->Fill(Muon_innerTrack_normalizedChi2->at(mu3));
            hNMatchedStat->Fill(Muon_numberOfMatchedStations->at(mu1));
            hNMatchedStat->Fill(Muon_numberOfMatchedStations->at(mu2));
            hNMatchedStat->Fill(Muon_numberOfMatchedStations->at(mu3));
            hmasstri->Fill(Triplet_Mass->at(ind));
            particleId(Muon_PdgId->at(mu1), IdsummaryDaughter_AC);
            particleId(Muon_PdgId->at(mu2), IdsummaryDaughter_AC);
            particleId(Muon_PdgId->at(mu3), IdsummaryDaughter_AC);
            particleId(Muon_MotherPdgId->at(mu1), IdsummaryMother_AC);
            particleId(Muon_MotherPdgId->at(mu2), IdsummaryMother_AC);
            particleId(Muon_MotherPdgId->at(mu3), IdsummaryMother_AC);
            hFlightDist->Fill(FlightDistPVSV->at(ind));
            hFlightDist_Signif->Fill(FlightDistPVSV_Significance->at(ind));
            hFlightDistvsP->Fill(FlightDistPVSV->at(ind), TripletP);
            DiMuonHisto(hmassdi, dimu1_2, dimu2_3, dimu1_3);
            //Aggiungo al tripletto una 4a traccia e ne faccio la massa inv
            QuadMuMass2(hmassQuad, hmassQuad_Zero, mu1, mu2, mu3);
        }
        
        if (ientry < 0) break;
    }//end loop sugli eventi
    cut[0] = ntripltot;
    //Stampo info generali & disegno gli istogrammi
    cout << endl;
    cout << "NUMERO EVENTI -> " << nentries << endl << endl;
    cout << "NUMERO TRIPLETTI -> " << ntripltot << endl << endl;
    cout << "Triplets survived: " << cutevt[7] << " || Good: " << NgoodTripl << " , Bad: " << NbadTripl << endl;
    //Histo AFTER cuts
    gStyle->SetOptStat(1101);
    hNtripl->GetXaxis()->SetTitle("N. triplets per event");
    hPt->GetXaxis()->SetTitle("Pt (GeV)");
    hPt_tripl->GetXaxis()->SetTitle("Pt (GeV)");
    hPt_mu1->GetXaxis()->SetTitle("Pt (GeV)");
    hPt_mu2->GetXaxis()->SetTitle("Pt (GeV)");
    hPt_mu3->GetXaxis()->SetTitle("Pt (GeV)");
    hPtcomp_mu1->GetXaxis()->SetTitle("(Mu1_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcomp_mu2->GetXaxis()->SetTitle("(Mu2_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcomp_mu3->GetXaxis()->SetTitle("(Mu3_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompBarrel_mu1->GetXaxis()->SetTitle("Barrel_(Mu1_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompBarrel_mu2->GetXaxis()->SetTitle("Barrel_(Mu2_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompBarrel_mu3->GetXaxis()->SetTitle("Barrel_(Mu3_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompEndcap_mu1->GetXaxis()->SetTitle("Endcap_(Mu1_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompEndcap_mu2->GetXaxis()->SetTitle("Endcap_(Mu2_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompEndcap_mu3->GetXaxis()->SetTitle("Endcap_(Mu3_Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hmasstri->GetXaxis()->SetTitle("Triplet mass (GeV)");
    hMassTriResolution->GetXaxis()->SetTitle("(Mass_{sim}-Mass_{reco})/Mass_{sim}");
    hMassTriResolutionBarrel->GetXaxis()->SetTitle("Barrel_(Mass_{sim}-Mass_{reco})/Mass_{sim}");
    hMassTriResolutionEndcap->GetXaxis()->SetTitle("Endcap_(Mass_{sim}-Mass_{reco})/Mass_{sim}");
    hmassdi->GetXaxis()->SetTitle("2 mu mass (GeV)");
    hmassQuad->GetXaxis()->SetTitle("4 mu mass (GeV)");
    hmassQuad_Zero->GetXaxis()->SetTitle("4 mu mass (GeV)");
    hNMatchedStat->GetXaxis()->SetTitle("N. of matched muon stations");
    hFlightDist->GetXaxis()->SetTitle("Decay length (cm)");
    hFlightDist_Signif->GetXaxis()->SetTitle("Decay length significance");
    TCanvas *PdgIdCanv_AC = new TCanvas("PdgId_Daughter_AC", "PdgId_Daughter_AC", 0, 0, 1600, 1000);
    PdgIdHisto(PdgIdCanv_AC, hPdgId_AC, IdsummaryDaughter_AC, pId);
    TCanvas *PdgIdMotherCanv_AC = new TCanvas("PdgId_Mother_AC", "PdgId_Mother_AC", 0, 0, 1600, 1000);
    PdgIdHisto(PdgIdMotherCanv_AC, hMotherPdgId_AC, IdsummaryMother_AC, pId);
    hFlightDistvsP->GetXaxis()->SetTitle("FlightDist");
    hFlightDistvsP->GetYaxis()->SetTitle("P (GeV)");
    hFlightDistvsP->DrawCopy("colz");
    //Histo BEFORE cuts
    hPtcomp_mu1_BC->GetXaxis()->SetTitle("Mu1_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcomp_mu2_BC->GetXaxis()->SetTitle("Mu2_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcomp_mu3_BC->GetXaxis()->SetTitle("Mu3_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcomp_BC->GetXaxis()->SetTitle("(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompBarrel_mu1_BC->GetXaxis()->SetTitle("Barrel_Mu1_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompBarrel_mu2_BC->GetXaxis()->SetTitle("Barrel_Mu2_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompBarrel_mu3_BC->GetXaxis()->SetTitle("Barrel_Mu3_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompBarrel_BC->GetXaxis()->SetTitle("Barrel_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompEndcap_mu1_BC->GetXaxis()->SetTitle("Endcap_Mu1_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompEndcap_mu2_BC->GetXaxis()->SetTitle("Endcap_Mu2_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompEndcap_mu3_BC->GetXaxis()->SetTitle("Endcap_Mu3_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hPtcompEndcap_BC->GetXaxis()->SetTitle("Endcap_(Pt_{sim}-Pt_{reco})/Pt_{sim}");
    hmasstriBeforeCuts->GetXaxis()->SetTitle("Triplet Mass (GeV)");
    hMassvsChi2->GetXaxis()->SetTitle("Mass (GeV)");
    hMassvsChi2->GetYaxis()->SetTitle("Chi2triplet");
    hMassvsChi2->DrawCopy("colz");
    hmassQuad_Other->GetXaxis()->SetTitle("4 mu mass (GeV)");
    hmassQuad_Other_Zero->GetXaxis()->SetTitle("4 mu mass (GeV)");
    hmassDi_Other_Zero->GetXaxis()->SetTitle("2 mu mass (GeV)");
    hmassDi_Other2_Zero->GetXaxis()->SetTitle("2 mu mass (GeV)");
    TCanvas *PdgIdCanv_BC = new TCanvas("PdgId_Daughter_BC", "PdgId_Daughter_BC", 0, 0, 1600, 1000);
    PdgIdHisto(PdgIdCanv_BC, hPdgId_BC, IdsummaryDaughter_BC, pId);
    TCanvas *PdgIdMotherCanv_BC = new TCanvas("PdgId_Mother_BC", "PdgId_Mother_BC", 0, 0, 1600, 1000);
    PdgIdHisto(PdgIdMotherCanv_BC, hMotherPdgId_BC, IdsummaryMother_BC, pId);
    //Histo intermediate
    hmassdi_int->GetXaxis()->SetTitle("2 mu mass (GeV)");
    hmassQuad_int->GetXaxis()->SetTitle("4 mu mass (GeV)");
    hmassQuad_Zero_int->GetXaxis()->SetTitle("4 mu mass (GeV)");
    //Histo of cuts Efficiency
    TCanvas *canvEvt = new TCanvas("CutEfficiency_Nevents", "CutEfficiency_Nevents", 0, 0, 1200, 1000);
    CutEffHisto(canvEvt, hCutEffEvt, cutevt, listCut);
    TCanvas *canv = new TCanvas("CutEfficiency_Ntriplets", "CutEfficiency_Ntriplets", 0, 0, 1200, 1000);
    CutEffHisto(canv, hCutEff, cut, listCut);
    //Write and close the file
    fout->Write();
    fout->Close();
}

// #########################################  END SIGNAL MC ANALYSIS  NEW CUTFLOW








