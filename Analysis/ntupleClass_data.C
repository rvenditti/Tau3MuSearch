#define ntupleClass_MC_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


// ######################################### DATA: SIGNAL ANALYSIS CUTFLOW

// Cuts: (over triplets)
// * cut[1] -> Chi2 triplet vertex (in 0 - 15)
// * cut[2] -> There are (2 mu glb w/ pt>ptmin (=2) e 1 mu Tracker con pt>ptminTrack (=0.5)) & |eta|<2.4
// * cut[3] -> Triplet mass (in 1.62 - 2 GeV)
// * cut[4] -> The 3 possible pairs of mu of the triplet have proper DeltaR (<0.8)
// * cut[5] -> The 3 possible pairs of mu of the triplet have proper |DeltaZ| (<0.5)
// * cut[6] -> Cut on the dimuon mass w.r.t. Phi(1020) per pairs of mu of the triplet w/ opposite sign
// * cut[7] -> Cut on the dimuon mass w.r.t. Omega(782) per pairs of mu of the triplet w/ opposite sign
// * cut[8] -> Trigger Matching
// * cut[9] -> Filling HistoStepByStep if invariant mass is in SB and HLT fired and L1Double fired
// * cut[10] -> Filling HistoStepByStep if invariant mass is in SB and HLT fired and L1Double fired
// * cut[11] -> Filling HistoStepByStep if HLT and L1Double OR L1Triple fired, no cut on invariant mass



void ntupleClass_MC::LoopData_New(TString type, TString datasetName){
    
    pileupFactor = 1;
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    // Variables definition
    int ntripl, trInd = 0, ind = 0, mu_Ind[NMU] = {0}, mu[NMU] = {0}, NgoodTripl = 0, NbadTripl = 0, cut[NCUTS] = {0}, cutevt[NCUTS] = {0}, triplIndex[1000] = {0}, Ncut = 0, IdsummaryDaughter[NCUTS][NPARTICLES] = {0}, IdsummaryMother[NCUTS][NPARTICLES] = {0}, IdsummaryDaughter_Gen[NPARTICLES] = {0}, IdsummaryMother_Gen[NPARTICLES] = {0};
    float ptminTrack = 0.5, DeltaRmax = 0.8, DeltaZmax = 0.5, DeltaZ1 = 0, DeltaZ2 = 0, DeltaZ3 = 0, ntripltot = 0;
    double dimu[NMU] = {0}, massmin = 1.62, massmax = 2.00, sigmaPhi = 0.011, sigmaOmega = 0.0085, TripletVtx_Chi2max = 15, EtaMax = 2.4;
    TString listCut[NCUTS], pId[NPARTICLES];
    // Variables for the final tree
    double Pmu3 = 0, cLP = 0, segmComp = 0, fv_nC = 0, fv_dphi3D = 0, fv_d3Dsig = 0, d0 = 0, d0sig = 0, mindca_iso = 0, trkRel = 0, Pmu1 = 0, Ptmu1 = 0, etamu1 = 0, Pmu2 = 0, Ptmu2 = 0, etamu2 = 0, Ptmu3 = 0, etamu3 = 0, P_trip = 0, Pt_trip = 0, eta_trip = 0, nStationsMu1 = 0, nStationsMu2 = 0, nStationsMu3 = 0, Iso03Mu1 = 0, Iso03Mu2 = 0, Iso03Mu3 = 0, Iso05Mu1 = 0, Iso05Mu2 = 0, Iso05Mu3 = 0, nMatchesMu1 = 0, nMatchesMu2 = 0, nMatchesMu3 = 0, timeAtIpInOutMu1 = 0, timeAtIpInOutMu2 = 0, timeAtIpInOutMu3 = 0, cQ_uS = 0, cQ_tK, cQ_gK = 0, cQ_tRChi2 = 0, cQ_sRChi2 = 0, cQ_Chi2LM = 0, cQ_Chi2lD = 0, cQ_gDEP = 0, cQ_tM = 0, cQ_gTP = 0, calEn_emMu1 = 0, calEn_emMu2 = 0, calEn_emMu3 = 0, calEn_hadMu1 = 0, calEn_hadMu2 = 0, calEn_hadMu3 = 0, caloComp = 0, fliDistPVSV_Chi2 = 0, isGlb3 = 0, isTracker3 = 0, isLoose3 = 0, isSoft3 = 0, isPF3 = 0, isRPC3 = 0, isSA3 = 0, isCalo3 = 0, vx1 = 0, vx2 = 0, vx3 = 0, vy1 = 0, vy2 = 0, vy3 = 0, vz1 = 0, vz2 = 0, vz3 = 0, Refvx1 = 0, Refvx2 = 0, Refvx3 = 0, Refvy1 = 0, Refvy2 = 0, Refvy3 = 0, Refvz1 = 0, Refvz2 = 0, Refvz3 = 0, SVx = 0, SVy = 0, SVz = 0, had03 = 0, had05 = 0, nJets03 = 0, nJets05 = 0, nTracks03 = 0, nTracks05 = 0, sumPt03 = 0, sumPt05 = 0, hadVeto03 = 0, hadVeto05 = 0, emVeto03 = 0, emVeto05 = 0, trVeto03 = 0, trVeto05 = 0;
    float tKink = 0;
    //Variables inizialization
    Fill_particleName(pId);
    cutevt[0] = nentries;
    Fill_CutName(listCut);
    // Creation of the output file & final tree
    TString root_fileName = fileName;
    TFile *fout = new TFile(root_fileName, "RECREATE");
    fout->cd();
    TTree *treeA = new TTree("FinalTreeA_Bkg","FinalTreeA_Bkg");
    TTree *treeB = new TTree("FinalTreeB_Bkg","FinalTreeB_Bkg");
    TTree *treeC = new TTree("FinalTreeC_Bkg","FinalTreeC_Bkg");
    TreeFin_Init(treeA, Pmu3, cLP, tKink, segmComp, fv_nC, fv_dphi3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
    TreeFin_Init(treeB, Pmu3, cLP, tKink, segmComp, fv_nC, fv_dphi3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
    TreeFin_Init(treeC, Pmu3, cLP, tKink, segmComp, fv_nC, fv_dphi3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
    // Creation of histograms for variables BEFORE cuts
    TDirectory *dirBeforeCuts = fout->mkdir("BeforeCuts");
    dirBeforeCuts->cd();
    TH1D *hPileUp_BC, *hNPrVert_BC, *hMass_tripl_BC, *hChi2Vertex, *hMass_di_Zero_BC, *hMass_di_Zero2_BC, *hPtRes_BC, *hPtRes_BC_mu[NMU], *hPtResBarrel_BC, *hPtResBarrel_BC_mu[NMU], *hPtResEndcap_BC, *hPtResEndcap_BC_mu[NMU]; TH2D *hMassvsChi2;
    TH1F *hMass_quad_BC, *hMass_quad_Zero_BC;
    TH1I *hPdgId_Gen, *hMotherPdgId_Gen; TH2I *hPdgId2D_Gen;
    hPileUp_BC = new TH1D("hNPileUp", "hNPileUp", 80, -0.5, 79.5);
    hPileUp_BC->Sumw2();
    hNPrVert_BC = new TH1D("hNPrimaryVertices", "hNPrimaryVertices", 100, -0.5, 99.5);
    hNPrVert_BC->Sumw2();
    TCanvas *PdgIdCanvas_Gen, *PdgIdMotherCanvas_Gen, *PdgIdCanvas2D_Gen;
    InitHistoBC(hMass_tripl_BC, hChi2Vertex, hMassvsChi2, hMass_quad_BC, hMass_quad_Zero_BC, hMass_di_Zero_BC, hMass_di_Zero2_BC, hPtRes_BC, hPtRes_BC_mu, hPtResBarrel_BC, hPtResBarrel_BC_mu, hPtResEndcap_BC, hPtResEndcap_BC_mu, hPdgId_Gen, hMotherPdgId_Gen, hPdgId2D_Gen);
    fout->cd();
    // Creation of histograms for variables AFTER cuts
    TDirectory *dirAfterCuts = fout->mkdir("AfterCuts");
    dirAfterCuts->cd();
    TH1I *hNtripl; TH1F *hChi2Track, *hmassQuad, *hmassQuad_Zero;
    TH1D *hPileUp_AC, *hNPrVert_AC, *hTripTriggerMatched, *hMassTriRes, *hMassTriResBarrel, *hMassTriResEndcap, *hmassdi, *hPtRes_AC, *hPtRes_AC_mu[NMU], *hPtResBarrel_AC, *hPtResBarrel_AC_mu[NMU], *hPtResEndcap_AC, *hPtResEndcap_AC_mu[NMU], *hNMatchedStat, *hFlightDist, *hFlightDist_Signif, *hPtErrOverPt, *hPt_tripl_good, *hPt_tripl_fake, *hDeltaX, *hDeltaY, *hDeltaZ, *hDeltaX_fake, *hDeltaY_fake, *hDeltaZ_fake, *hTripMassA, *hTripMassB, *hTripMassC, *hEtaA, *hEtaB, *hEtaC; TH2D *hFlightDistvsP;
    hPileUp_AC = new TH1D("hNPileUp", "hNPileUp", 80, -0.5, 79.5);
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
    InitHistoAC(hNtripl, hChi2Track, hMassTriRes, hMassTriResBarrel, hMassTriResEndcap, hmassdi, hmassQuad, hmassQuad_Zero, hPtRes_AC, hPtRes_AC_mu, hPtResBarrel_AC, hPtResBarrel_AC_mu, hPtResEndcap_AC, hPtResEndcap_AC_mu, hNMatchedStat, hFlightDist, hFlightDist_Signif, hFlightDistvsP, hPtErrOverPt, hPt_tripl_good, hPt_tripl_fake, hDeltaX, hDeltaY, hDeltaZ, hDeltaX_fake, hDeltaY_fake, hDeltaZ_fake);
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
    TDirectory *dirTriplet = dirStepByStep->mkdir("Triplet");  // Triplet variables histo
    dirTriplet->cd();
    TH1D *hPt_tripl[NCUTS], *hEta_tripl[NCUTS], *hPhi_tripl[NCUTS], *hMass_tripl[NCUTS];
    InitHistoStepByStep_Triplet(hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl);
    dirStepByStep->cd();
    TDirectory *dirPdgId = dirStepByStep->mkdir("PdgId"); // PdgId histograms
    dirPdgId->cd();
    TH1I *hPdgId_cut[NCUTS], *hMotherPdgId_cut[NCUTS]; TH2I *hPdgId2D_cut[NCUTS];
    TCanvas *PdgIdCanvas_cut[NCUTS], *PdgIdMotherCanvas_cut[NCUTS], *PdgIdCanvas2D_cut[NCUTS];
    InitHistoStepByStep_PdgId(hPdgId_cut, hMotherPdgId_cut, hPdgId2D_cut);
    dirStepByStep->cd();
    fout->cd();
    // Creation of total histograms
    TH1I *hCutEff = new TH1I("CutEff_Ntriplets", "CutEff_Ntriplets", NCUTS, 0.5, (NCUTS+0.5));
    TH1I *hCutEffEvt = new TH1I("CutEff_NEvents", "CutEff_NEvents", NCUTS, 0.5, (NCUTS+0.5));
    
    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        ntripl = 0, trInd = 0; int cutevt2[NCUTS] = {0};
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetTree()->GetEntry(ientry);
        hPileUp_BC->Fill(nPileUpInt);
        hNPrVert_BC->Fill(PVCollection_Size);
        //Loop over the TRIPLETS
        for (int j=0; j<TripletVtx_Chi2->size(); j++){
            //Check the number of tracks in the primary vertex
            if(RefittedPV_NTracks->at(j) > 1){
                ntripltot++; Ncut = 0;
                // BEFORE cuts
                //Matching between index of single mu of the triplet (mu#_Ind) & that of  'MUONID' (mu#)
                MatchIndex("ID", j, mu_Ind, mu);
                // Fill histograms
                FillHistoBC("data", j, hMass_tripl_BC, hChi2Vertex, hMassvsChi2, hMass_quad_BC, hMass_quad_Zero_BC, hMass_di_Zero_BC, hMass_di_Zero2_BC, hPtRes_BC, hPtRes_BC_mu, hPtResBarrel_BC, hPtResBarrel_BC_mu, hPtResEndcap_BC, hPtResEndcap_BC_mu, IdsummaryDaughter_Gen, IdsummaryMother_Gen, Idsummary2D_Gen);
                FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                //CUT 1 : check condition on * Chi2 vertex ( 0 < Chi2 < 15)
                if (TripletVtx_Chi2->at(j) < TripletVtx_Chi2max && TripletVtx_Chi2->at(j) > 0 ){
                    Ncut++; cut[Ncut]++; cutevt2[Ncut]++;
                    FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                    // CUT 2 :
                    // Check that mu1 is glb & pt>ptmax & |eta|<Etamax
                    if((Muon_isLoose->at(mu[0]) == 1) && (MuonPt->at(mu[0]) > ptmin) && abs(Mu1_Eta->at(mu_Ind[0])) < EtaMax){
                        DeltaZ1 = Muon_vz->at(mu[0]);
                        // Check that mu2 is glb & pt>ptmax & |eta|<Etamax
                        if((Muon_isLoose->at(mu[1]) == 1) && (MuonPt->at(mu[1]) > ptmin) && abs(Mu2_Eta->at(mu_Ind[1])) < EtaMax){
                            DeltaZ2 = Muon_vz->at(mu[1]);
                            // Check that mu3 is tracker & pt>0.5 & |eta|<Etamax
                            if((Muon_isLoose->at(mu[2]) == 1) && (MuonPt->at(mu[2]) > ptminTrack) && abs(Mu3_Eta->at(mu_Ind[2])) < EtaMax){
                            //if((Muon_isGlobal->at(mu[2]) == 1) && (MuonPt->at(mu[2]) > ptminTrack) && abs(Mu3_Eta->at(mu_Ind[2])) < EtaMax){
                                DeltaZ3 = Muon_vz->at(mu[2]);
                                Ncut++; cut[Ncut]++; cutevt2[Ncut]++;
                                FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                                //CUT 3: check condition on trimuon mass
                                if(Triplet_Mass->at(j) > massmin && Triplet_Mass->at(j) < massmax){
                                    Ncut++; cut[Ncut]++; cutevt2[Ncut]++;
                                    FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                                    //CUT 4: Loop on PAIRS of muons of the triplet & check DeltaR
                                    if(isPairDeltaRGood(j, DeltaRmax) == true){
                                        Ncut++; cut[Ncut]++; cutevt2[Ncut]++;
                                        FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                                        //CUT 5 : Check |Delta Z|
                                        if(isPairDeltaZGood(DeltaZ1, DeltaZ2, DeltaZ3, DeltaZmax) == true){
                                            Ncut++; cut[Ncut]++; cutevt2[Ncut]++;
                                            FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                                            //CUT 6: VETO on Phi(1020) mass
                                            Fill_DimuonMass(mu_Ind, mu, dimu);
                                            if(isPairNotAPhi(dimu, sigmaPhi) == true){
                                                Ncut++; cut[Ncut]++; cutevt2[Ncut]++;
                                                FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                                                //CUT 7: VETO on Omega(782) mass
                                                if(isPairNotAOmega(dimu, sigmaOmega) == true){
                                                    Ncut++; cut[Ncut]++; cutevt2[Ncut]++;
                                                    FillHistoStepByStep("data", j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
                                                    //CUT 8: Trigger Matching
                                                    if(Mu1_dRtriggerMatch->at(j)<0.03 && Mu2_dRtriggerMatch->at(j)<0.03 && Mu3_dRtriggerMatch->at(j)<0.03){
                                                       Ncut++; ntripl++; triplIndex[trInd] = j; trInd++;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } // end loop on triplets
        }
        // N. events that passed each selection
        for (int k=1; k<NCUTS; k++){
            if(cutevt2[k] > 0) cutevt[k]++;
        }
        // Histo N. triplets passed for each event
        hNtripl->Fill(ntripl);
        if(ntripl > 0) {
            cutevt[NCUTS-1]++; cut[NCUTS-1]++;
            ind = BestTripletFinder(triplIndex, ntripl);
            //RiMatching between index of single mu of the triplet (mu#_Ind) & that of  'MUONID' (mu#) & Ricomputing the 3 possible dimuon masses
            MatchIndex("ID", ind, mu_Ind, mu);

            //Check HLT and L1 decision
            bool hlt_fired = 0;
            bool l1double_fired = 0;
            bool l1triple_fired = 0;
            cout<<"evt "<<evt<<endl;
            for(int h=0; h<Trigger_hltname->size(); h++) {
               TString hltName = Trigger_hltname->at(h);
               if(strncmp(hltName, "HLT_DoubleMu3_Trk_Tau3mu_v", 26) == 0 && Trigger_hltdecision->at(h) == 1) {
                  hlt_fired = 1;
               }
            }
            cout<<" ========= debug ========= "<<endl;
            for(int k=0; k<Trigger_l1name->size(); k++) {
               TString l1Name = Trigger_l1name->at(k);
               cout<<l1Name<<endl;
               //HLT + L1_DoubleMu
               if(strcmp(l1Name, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") == 0 && Trigger_l1decision->at(k) == 1) l1double_fired = 1;
               //HLT + L1_TripleMu
               if( ( strcmp(l1Name, "L1_TripleMu_5_3_0_DoubleMu_5_3_OS_Mass_Max17") == 0 || strcmp(l1Name, "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9") == 0 ) && Trigger_l1decision->at(k) == 1)  l1triple_fired = 1;
            }

            bool isTrigger_forAna = 0;
            if( hlt_fired == 1 && ( l1double_fired == 1 || l1triple_fired == 1 ) ) isTrigger_forAna = 1;

            bool isSB_tripletMass = 0;
            if( (Triplet_Mass->at(ind) >= 1.65 && Triplet_Mass->at(ind) <= 1.73) || (Triplet_Mass->at(ind) >= 1.82 && Triplet_Mass->at(ind) <= 1.90) ) isSB_tripletMass = 1;
            
           //if(strcmp(type, "data_bkg") == 0){
            if(strcmp(type, "data_bkg") == 0 && isTrigger_forAna && isSB_tripletMass){
                double tripReso = ResoTriplMass(mu_Ind, mu);
                if(tripReso < catA){
                    hTripMassA->Fill(Triplet_Mass->at(ind));
                    hEtaA->Fill(abs(Triplet_Eta->at(ind)), pileupFactor);
                    TreeFin_Fill(treeA, ind, mu_Ind, mu, Pmu3, cLP, tKink, segmComp, fv_nC, fv_dphi3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
                }
                else if (tripReso >= catA && tripReso <= catB){
                    hTripMassB->Fill(Triplet_Mass->at(ind));
                    hEtaB->Fill(abs(Triplet_Eta->at(ind)), pileupFactor);
                    TreeFin_Fill(treeB, ind, mu_Ind, mu, Pmu3, cLP, tKink, segmComp, fv_nC, fv_dphi3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
                }
                else if(tripReso > catB){
                    hTripMassC->Fill(Triplet_Mass->at(ind));
                    hEtaC->Fill(abs(Triplet_Eta->at(ind)), pileupFactor);
                    TreeFin_Fill(treeC, ind, mu_Ind, mu, Pmu3, cLP, tKink, segmComp, fv_nC, fv_dphi3D, fv_d3Dsig, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3, isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
                }
            }
            
            Fill_DimuonMass(mu_Ind, mu, dimu);

            if( isSB_tripletMass && hlt_fired == 1 && l1double_fired == 1  ){
               //plot bkg for yield l1DoubleMu
               FillHistoStepByStep("data", ind, mu_Ind, mu, NCUTS-3, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
            }
            if( isSB_tripletMass && hlt_fired == 1 && l1triple_fired == 1  ){
               //plot bkg for yield l1TripleMu
               FillHistoStepByStep("data", ind, mu_Ind, mu, NCUTS-2, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
            }
            if ( isTrigger_forAna ) {
               FillHistoStepByStep("data", ind, mu_Ind, mu, NCUTS-1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
               for(int k=0; k<NMU; k++){
                   hIsolation_03->Fill(Muon_emEt03->at(mu[k]));
                   hIsolation_05->Fill(Muon_emEt05->at(mu[k]));
               }
               FillHistoResoTriplMass(mu_Ind, mu, hMassTriRes, hMassTriResBarrel, hMassTriResEndcap);
               FillHistoAC(ind, mu, hChi2Track, hNMatchedStat, hFlightDist, hFlightDist_Signif, hFlightDistvsP, hPtErrOverPt, hmassdi, dimu, hmassQuad, hmassQuad_Zero);
               hPileUp_AC->Fill(nPileUpInt);
               hNPrVert_AC->Fill(PVCollection_Size);
            }
            // Filling histogram with triplet mass if hlt is fired AND (l1_DoubleMu || l1_TripleMu) is fireds
            TriggerRequirements(ind, hTripTriggerMatched);
        }
        if (ientry < 0) break;
    }//end loop on events
    cut[0] = ntripltot;
    //Print general info
    cout << endl;
    cout << "N. EVENTS -> " << nentries << endl << endl;
    cout << "N. TRIPLETS -> " << ntripltot << endl << endl;
    //Histo of cuts Efficiency
//    TCanvas *canvEvt = new TCanvas("CutEfficiency_Nevents", "CutEfficiency_Nevents", 0, 0, 1200, 1000);
//    Draw_CutEffCanvas(canvEvt, hCutEffEvt, cutevt, listCut);
//    TCanvas *canv = new TCanvas("CutEfficiency_Ntriplets", "CutEfficiency_Ntriplets", 0, 0, 1200, 1000);
//    Draw_CutEffCanvas(canv, hCutEff, cut, listCut);
    //Write and close the file
    fout->Write();
    fout->Close();
}

// #########################################

