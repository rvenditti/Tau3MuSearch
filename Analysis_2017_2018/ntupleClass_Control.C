#define ntupleClass_Control_cxx
#define NCUTS 16
#define NPARTICLES 560
#define NMU_C 2
#define NTOT 3
#define mumass 0.1056583715 // Muon mass in GeV
#define PhiMass 1.019461 // Phi mass in GeV
#define OmegaMass 0.78265 // Omega mass in GeV
#define ptmin 2.0

#include "ntupleClass_Control.h"
#include "Utilities_Control.C"
#include "PdgId_list_Control.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// #########################################    DsPhiPi ANALYSIS CUTFLOW

// Cuts: (over events)
// * cut[0] -> Before cuts
// * cut[1] -> Event fires L1 seed
// * cut[2] -> Event fires HLT
// Cuts: (over triplets)
// * cut[3] -> all triplets w/ at least 2 track associated with PV
// * cut[4] -> Both mu have to be different from the track
// * cut[5] -> Chi2 triplet vertex (in 0 - 15)
// * cut[6] -> 2 mu w/ opposite charge
// * cut[7] -> Invariant dimuon mass in (1-1.04) GeV
// * cut[8] -> Longitudianl IP of the track w.r.t. the beam spot < 20 cm
// * cut[9] -> Transverse IP of the track w.r.t. the beam spot < 0.3 cm
// * cut[10] -> Null
// * cut[11] -> Mu1 TriggerMatching
// * cut[12] -> Mu2 TriggerMatching
// * cut[13] -> Trk TriggerMatching
// Cuts: (over events containing at least 1 triplet)
// * cut[14] -> cut on invariant mass: sgn
// * cut[15] -> cut on invariant mass: bkg
// * cut[16] -> no cut on invariant mass
//
// N.B.: cut[NCUTS] total number of triplets passing each selection
//       cutevt[NCUTS] total number of events passing each selection
//       cuttripl[NCUTS] number of triplets passing each selection in current event

void ntupleClass_Control::LoopControl(TString type, TString datasetName){

    bool isMC = false;
    if(strcmp(type, "MC") == 0 ) isMC = true;

    bool isVerbose = false;
    bool isMiniAOD = true;

    bool doPUrew = true;
    if(!isMC) doPUrew = false;

     if( isMC && doPUrew ) {
        // Pile-up reweighting
        // data pileup from 2017 era C-D-E
        TFile *fPileUp = new TFile("/lustrehome/fsimone/Analysis/Pile_up_reweighing_tools/PileUp_ReweightingStudy_2017CDE_DsPhiPi.root");
        TH1 *hPileUpRew = (TH1*)fPileUp->Get("PileUp_Reweighting");
        int Nbins = hPileUpRew->GetNbinsX();
        for(int m=0; m<Nbins; m++){
            pileup_weight.push_back(hPileUpRew->GetBinContent(m));
        }
    }
    // End pile-up reweighting

    // opening file containing ScaleFactors
    TFile *_f = new TFile("/lustrehome/fsimone/Analysis/ScaleFactor/Tau3MuMuonSF.root");
    TH2F* SF_h = dynamic_cast<TH2F*> (_f->Get("LooseID_abseta_pt_Bins"));

    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntries();
    // Variables definition
    int ntripl, trInd = 0, ind = 0, mu_Ind[NTOT] = {0}, mu[NTOT] = {0}, muGen[NTOT] = {0}, NgoodTripl = 0, NbadTripl = 0, cut[NCUTS] = {0}, cutevt[NCUTS] = {0}, Ncut = 0, IdsummaryDaughter[NCUTS][NPARTICLES] = {0}, IdsummaryMother[NCUTS][NPARTICLES] = {0}, IdsummaryDaughter_Gen[NPARTICLES] = {0}, IdsummaryMother_Gen[NPARTICLES] = {0};
    float ptminTrack = 0.5, DeltaRmax = 0.8, DeltaZmax = 0.5, DeltaZ1 = 0, DeltaZ2 = 0, DeltaZ3 = 0;
    double dimu[NTOT] = {0};
    double massmin = 1.62, massmax = 2.00, sigma = 0.011, TripletVtx_Chi2max = 15, EtaMax = 2.4;
    TString pId[NPARTICLES], listCut[NCUTS];

    // Variables for the final tree
    double Pmu3 = 0, cLP = 0, segmComp = 0, tripletMass = 0, tripletMassReso = 0, fv_nC = 0, fv_dphi3D = 0, fv_d3D = 0, fv_d3Dsig = 0, bs_sv_d3D = 0, bs_sv_d3Dsig = 0, pv_sv_dxy = 0, pv_sv_dxy_sig = 0, d0 = 0, d0sig = 0, mindca_iso = 0, trkRel = 0, Pmu1 = 0, Ptmu1 = 0, etamu1 = 0, Pmu2 = 0, Ptmu2 = 0, etamu2 = 0, Ptmu3 = 0, etamu3 = 0, P_trip = 0, Pt_trip = 0, eta_trip = 0, nStationsMu1 = 0, nStationsMu2 = 0, nStationsMu3 = 0, Iso03Mu1 = 0, Iso03Mu2 = 0, Iso03Mu3 = 0, Iso05Mu1 = 0, Iso05Mu2 = 0, Iso05Mu3 = 0, nMatchesMu1 = 0, nMatchesMu2 = 0, nMatchesMu3 = 0, timeAtIpInOutMu1 = 0, timeAtIpInOutMu2 = 0, timeAtIpInOutMu3 = 0, cQ_uS = 0, cQ_tK, cQ_gK = 0, cQ_tRChi2 = 0, cQ_sRChi2 = 0, cQ_Chi2LM = 0, cQ_Chi2lD = 0, cQ_gDEP = 0, cQ_tM = 0, cQ_gTP = 0, calEn_emMu1 = 0, calEn_emMu2 = 0, calEn_emMu3 = 0, calEn_hadMu1 = 0, calEn_hadMu2 = 0, calEn_hadMu3 = 0, caloComp = 0, fliDistPVSV_Chi2 = 0, vx1 = 0, vx2 = 0, vx3 = 0, vy1 = 0, vy2 = 0, vy3 = 0, vz1 = 0, vz2 = 0, vz3 = 0, dxy1 = 0, dxy2 = 0, dxy3 = 0, dxyErr1 = 0, dxyErr2 = 0, dxyErr3 = 0, Refvx1 = 0, Refvx2 = 0, Refvx3 = 0, Refvy1 = 0, Refvy2 = 0, Refvy3 = 0, Refvz1 = 0, Refvz2 = 0, Refvz3 = 0, SVx = 0, SVy = 0, SVz = 0, had03 = 0, had05 = 0, nJets03 = 0, nJets05 = 0, nTracks03 = 0, nTracks05 = 0, sumPt03 = 0, sumPt05 = 0, hadVeto03 = 0, hadVeto05 = 0, emVeto03 = 0, emVeto05 = 0, trVeto03 = 0, trVeto05 = 0;
    double isGlb1 = 0, isTracker1 = 0, isLoose1 = 0, isSoft1 = 0, isPF1 = 0, isRPC1 = 0, isSA1 = 0, isCalo1 = 0, isMedium1 = 0, muID1 = 0;
    double isGlb2 = 0, isTracker2 = 0, isLoose2 = 0, isSoft2 = 0, isPF2 = 0, isRPC2 = 0, isSA2 = 0, isCalo2 = 0, isMedium2 = 0, muID2 = 0;
    double isGlb3 = 0, isTracker3 = 0, isLoose3 = 0, isSoft3 = 0, isPF3 = 0, isRPC3 = 0, isSA3 = 0, isCalo3 = 0, isMedium3 = 0, muID3 = 0;
    double tKink = 0;
    double run_n = 0, lumi_n = 0, evt_n = 0;
    double mu_pt = 0, mu_eta = 0, mu_phi = 0, mu_energy = 0, mu_charge = 0, mu_isGlobal = 0, mu_isSoft = 0, mu_isLoose = 0, mu_isTight = 0, mu_isPF = 0, mu_isRPC = 0, mu_isStandAlone = 0, mu_isTracker = 0, mu_isCalo = 0, mu_isQualityValid = 0, mu_SoftMVA = 0, mu_isTimeValid = 0, mu_isIsolationValid = 0, mu_numberOfMatchedStations = 0, mu_numberOfMatches = 0, mu_timeAtIpInOut = 0, mu_timeAtIpInOutErr = 0, mu_GLnormChi2 = 0, mu_GLhitPattern_numberOfValidMuonHits = 0, mu_trackerLayersWithMeasurement = 0, mu_Numberofvalidpixelhits = 0, mu_Numberofvalidtrackerhits = 0, mu_outerTrack_p = 0, mu_outerTrack_eta = 0, mu_outerTrack_phi = 0, mu_outerTrack_normalizedChi2 = 0, mu_outerTrack_muonStationsWithValidHits = 0, mu_innerTrack_p = 0, mu_innerTrack_eta = 0, mu_innerTrack_phi = 0, mu_innerTrack_validFraction = 0, mu_innerTrack_highPurity = 0, mu_innerTrack_normalizedChi2 = 0, mu_QInnerOuter = 0, mu_combinedQuality_updatedSta = 0, mu_combinedQuality_trkKink = 0, mu_combinedQuality_glbKink = 0, mu_combinedQuality_trkRelChi2 = 0, mu_combinedQuality_staRelChi2 = 0, mu_combinedQuality_chi2LocalPosition = 0, mu_combinedQuality_chi2LocalMomentum = 0, mu_combinedQuality_localDistance = 0, mu_combinedQuality_globalDeltaEtaPhi = 0, mu_combinedQuality_tightMatch = 0, mu_combinedQuality_glbTrackProbability = 0, mu_IP3D_BS = 0, mu_IP2D_BS = 0, mu_IP3D_PV = 0, mu_IP2D_PV = 0, mu_validMuonHitComb,  mu_calEnergy_em = 0, mu_calEnergy_emS9 = 0, mu_calEnergy_emS25 = 0, mu_calEnergy_had = 0, mu_calEnergy_hadS9 = 0, mu_segmentCompatibility = 0, mu_caloCompatibility = 0, mu_ptErrOverPt = 0, mu_BestTrackPt = 0, mu_BestTrackPtErr = 0, mu_BestTrackEta = 0, mu_BestTrackEtaErr = 0, mu_BestTrackPhi = 0, mu_BestTrackPhiErr = 0, mu_emEt03 = 0, mu_hadEt03 = 0, mu_nJets03 = 0, mu_nTracks03 = 0, mu_sumPt03 = 0, mu_hadVetoEt03 = 0, mu_emVetoEt03 = 0, mu_trackerVetoPt03 = 0, mu_emEt05 = 0, mu_hadEt05 = 0, mu_nJets05 = 0, mu_nTracks05 = 0, mu_sumPt05 = 0, mu_hadVetoEt05 = 0, mu_emVetoEt05 = 0, mu_trackerVetoPt05 = 0;

    //Variables inizialization
    Fill_particleName(pId);
    Fill_CutName(listCut);
    // Creation of the output file
    TString root_fileName = fileName;
    TFile *fout = new TFile(root_fileName, "RECREATE");
    fout->cd();
    TTree *tree = new TTree("FinalTree_Control","FinalTree_Control");
    TreeFin_Init(tree, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, bs_sv_d3Dsig, bs_sv_d3D, pv_sv_dxy_sig, pv_sv_dxy, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
    // Creation of histograms for variables BEFORE cuts
    TDirectory *dirBeforeCuts = fout->mkdir("BeforeCuts");
    dirBeforeCuts->cd();
    TH1D *hPileUp_BC, *hNPrVert_BC, *hMass_tripl_BC, *hChi2Vertex, *hMass_di_Zero_BC, *hMass_di_Zero2_BC, *hPtRes_BC, *hPtRes_BC_mu[NMU_C], *hPtResBarrel_BC, *hPtResBarrel_BC_mu[NMU_C], *hPtResEndcap_BC, *hPtResEndcap_BC_mu[NMU_C]; TH2D *hMassvsChi2;
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
    TH1D *hPileUp_AC, *hNPrVert_AC, *hTripTriggerMatched, *hMassTriRes, *hMassTriResBarrel, *hMassTriResEndcap, *hmassdi, *hPtRes_AC, *hPtRes_AC_mu[NMU_C], *hPtResBarrel_AC, *hPtResBarrel_AC_mu[NMU_C], *hPtResEndcap_AC, *hPtResEndcap_AC_mu[NMU_C], *hNMatchedStat, *hFlightDist, *hFlightDist_Signif, *hPtErrOverPt, *hPt_tripl_good, *hPt_tripl_fake, *hDeltaX, *hDeltaY, *hDeltaZ, *hDeltaX_fake, *hDeltaY_fake, *hDeltaZ_fake, *hChi2VertexNorm, *hSegmComp, *hDeltaR, *hTrIPSign;
    TH1F *hMinSegmComp;
    TH2D *hFlightDistvsP;
    hPileUp_AC = new TH1D("hNPileUp", "hNPileUp", 80, -0.5, 79.5);
    hPileUp_AC->Sumw2();
    hNPrVert_AC = new TH1D("hNPrimaryVertices", "hNPrimaryVertices", 100, -0.5, 99.5);
    hNPrVert_AC->Sumw2();
    hTripTriggerMatched = new TH1D("hTriplMassTrigMatched", "hTriplMassTrigMatched", 42, 1.60, 2.02); // binning 10 MeV
    hTripTriggerMatched->Sumw2();
    TH1D *hPropDecayL_AC = new TH1D("hProperDecayLength", "hProperDecayLength", 100, -0.5, 99.5);
    hPropDecayL_AC->Sumw2();
    InitHistoAC("", hNtripl, hMinSegmComp, hChi2Track, hMassTriRes, hMassTriResBarrel, hMassTriResEndcap, hmassdi, hmassQuad, hmassQuad_Zero, hPtRes_AC, hPtRes_AC_mu, hPtResBarrel_AC, hPtResBarrel_AC_mu, hPtResEndcap_AC, hPtResEndcap_AC_mu, hNMatchedStat, hFlightDist, hFlightDist_Signif, hFlightDistvsP, hPtErrOverPt, hPt_tripl_good, hPt_tripl_fake, hDeltaX, hDeltaY, hDeltaZ, hDeltaX_fake, hDeltaY_fake, hDeltaZ_fake, hChi2VertexNorm, hSegmComp, hDeltaR, hTrIPSign);
    fout->cd();
    // Creation of histo StepByStep
    TDirectory *dirStepByStep = fout->mkdir("StepByStep");
    dirStepByStep->cd();
    TDirectory *dirSingleMu = dirStepByStep->mkdir("SingleMu"); // Single mu variables histo
    dirSingleMu->cd();
    TH1D *hPt[NCUTS], *hPt_mu[NCUTS][NMU_C], *hEta[NCUTS], *hEta_mu[NCUTS][NMU_C], *hPhi[NCUTS], *hVx[NCUTS], *hVy[NCUTS], *hVz[NCUTS], *hPt_Tr[NCUTS], *hEta_Tr[NCUTS];
    InitHistoStepByStep_SingleMu(hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr);
    dirStepByStep->cd();
    TDirectory *dirTriplet = dirStepByStep->mkdir("Triplet");  // Triplet variables histo
    dirTriplet->cd();
    TH1D *hPt_tripl[NCUTS], *hEta_tripl[NCUTS], *hPhi_tripl[NCUTS], *hMass_tripl[NCUTS];
    InitHistoStepByStep_Triplet(hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl);
    dirStepByStep->cd();
    TDirectory *dirPdgId = dirStepByStep->mkdir("PdgId"); // PdgId histo
    dirPdgId->cd();
    TH1I *hPdgId_cut[NCUTS], *hMotherPdgId_cut[NCUTS]; TH2I *hPdgId2D_cut[NCUTS];
    TCanvas *PdgIdCanvas_cut[NCUTS], *PdgIdMotherCanvas_cut[NCUTS], *PdgIdCanvas2D_cut[NCUTS];
    InitHistoStepByStep_PdgId(hPdgId_cut, hMotherPdgId_cut, hPdgId2D_cut);
    dirStepByStep->cd();
    fout->cd();
    // Creation of total histograms
    TH1I *hCutEff = new TH1I("CutEff_Ntriplets", "CutEff_Ntriplets", NCUTS , 0.5, (NCUTS+0.5));
    TH1I *hCutEffEvt = new TH1I("CutEff_NEvents", "CutEff_NEvents", NCUTS , 0.5, (NCUTS+0.5));
   
    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        std::vector<Int_t> triplIndex;
        ntripl = 0; trInd = 0; int cuttripl[NCUTS] = {0};
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetTree()->GetEntry(ientry);

        //assigning factor for pileup re-weighting
        pileupFactor = 1;
        if(doPUrew){
            if(nPileUpInt < 80) pileupFactor = pileup_weight.at(nPileUpInt);
            else pileupFactor = pileup_weight.at(79);
            if(isVerbose){
                cout << "nPileUpInt " << nPileUpInt << endl;
                cout << "PileUpFactor : " << pileupFactor << endl << endl;
            }
        }
        hPileUp_BC->Fill(nPileUpInt);
        hNPrVert_BC->Fill(PVCollection_Size, pileupFactor);

        //Skip event if no good triplets
        if(NGoodTriplets->at(0) == 0) continue;
        if(isVerbose) cout<<"=================================\nevt "<<evt<<" run "<<run<<" lumi "<<lumi<<endl;
   
        //Different branches for 2017 triggerObject
        std::vector< std::array<double, 3> > Muon_HLT2017;
        if( isMiniAOD && datasetName.Contains("2017") && (MuonPt_HLT2017->size()>0) ) {
            //Define vector of arrays (Pt, Eta, Phi) for trigger objects in the event 
            for(std::size_t k=0; k<MuonPt_HLT2017->size(); k++){
                std::array<double, 3> temp;
                temp[0] = MuonPt_HLT2017->at(k); temp[1] = MuonEta_HLT2017->at(k); temp[2] = MuonPhi_HLT2017->at(k);
                Muon_HLT2017.push_back(temp);
            }
            if(isVerbose){
                cout<<"\n"<<MuonPt_HLT2017->size()<<" 2017 TriggerObjects in the event:\n | Pt | Eta | Phi "<<endl;
                for( auto const& mu: Muon_HLT2017 ){
                    for( auto const& var: mu ){
                         cout<<" | "<<var;
                    } cout<<"\n";
                } cout<<"\n";
            }
            //Sort the vector and removes the duplicates
            std::sort(Muon_HLT2017.begin(), Muon_HLT2017.end());
            Muon_HLT2017.erase(std::unique(Muon_HLT2017.begin(), Muon_HLT2017.end()), Muon_HLT2017.end());
            if(isVerbose){
                cout<<"After sorting and removing duplicates "<<endl;
                for( auto const& mu: Muon_HLT2017 ){
                    for( auto const& var: mu ){
                         cout<<" | "<<var;
                    } cout<<"\n";
                } cout<<"\n";
            }
        }
 
        //Check HLT and L1 decision
        bool hlt_fired = 0;
        bool l1double_fired = 0;
        for(int h=0; h<Trigger_hltname->size(); h++) {
            TString hltName = Trigger_hltname->at(h);
            //HLT 2017
            if( (datasetName.Contains("2017") != std::string::npos) && strncmp(hltName, "HLT_DoubleMu3_Trk_Tau3mu_v", 26) == 0 && Trigger_hltdecision->at(h) == 1) {
                hlt_fired = 1;
            }
        }
        for(int k=0; k<Trigger_l1name->size(); k++) {
            TString l1Name = Trigger_l1name->at(k);
            //2017 L1_DoubleMu
            if( (
                  strcmp(l1Name, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4") == 0 ||//2017 + 2018
                  strcmp(l1Name, "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4") == 0 ||//2017 + 2018 Backup
                  strcmp(l1Name, "L1_DoubleMu4_SQ_OS_dR_Max1p2") == 0 ||     //2017 + 2018
                  strcmp(l1Name, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2") == 0 )    //2017 + 2018 Backup
                  && Trigger_l1decision->at(k) == 1){
                l1double_fired = 1;
                if(isVerbose) cout<<"L1 fired: "<<l1Name<<endl;
            //    if ( hlt_fired == 1 ) h_l1prescaleDouble->Fill(Trigger_l1prescale->at(k));
            }
        }
        if(isVerbose) cout<<"l1double_fired "<<l1double_fired<<endl;
        bool isTrigger_forAna = 0;
        if( hlt_fired == 1 &&  l1double_fired == 1 ) isTrigger_forAna = 1;


        //Loop over the TRIPLETS
        if(isVerbose) cout<<"Triplets in the event "<<TripletVtx2_Chi2->size()<<endl;
        for (int j=0; j<TripletVtx2_Chi2->size(); j++){
            //Matching between index of the 2 mu (mu#_Ind) & that of 'MUONID' (mu#)
            //and that of the track of the triplet
            MatchIndex("ID", j, mu_Ind, mu);

            if(isVerbose){
                cout<<"\ntriplet candidate "<<j<<"\n | index| Pt | Eta | Phi"<<endl;
                cout<<" | "<<mu[0]<<" | "<<Mu01_Pt->at(j)<<" | "<<Mu01_Eta->at(j)<<" | "<<Mu01_Phi->at(j)<<endl;
                cout<<" | "<<mu[1]<<" | "<<Mu02_Pt->at(j)<<" | "<<Mu02_Eta->at(j)<<" | "<<Mu02_Phi->at(j)<<endl;
                cout<<" | "<<mu[2]<<" | "<<Tr_Pt->at(j)<<" | "<<Tr_Eta->at(j)<<" | "<<Tr_Phi->at(j)<<endl;
                cout<<" | invariant mass: "<<Triplet2_Mass->at(j)<<endl;
                cout<<" | vertex chi2: "<<TripletVtx2_Chi2->at(j)<<endl;
                cout<<" | muon IDs \n | isPF | isGlobal | isTrackerMuon "<<endl;
                cout<<" | "<<Muon_isPF->at(mu[0])<<" | "<<Muon_isGlobal->at(mu[0])<<" | "<<Muon_isTrackerMuon->at(mu[0])<<endl;
                cout<<" | "<<Muon_isPF->at(mu[1])<<" | "<<Muon_isGlobal->at(mu[1])<<" | "<<Muon_isTrackerMuon->at(mu[1])<<endl;
                //cout<<" | "<<Muon_isPF->at(mu[2])<<" | "<<Muon_isGlobal->at(mu[2])<<" | "<<Muon_isTrackerMuon->at(mu[2])<<endl;
                if(isMC && !(datasetName.Contains("2016"))){
                    cout<<" | simInfo:\n | mu1_pdgID | mu2_pdgId | trk_pdgId "<<endl;
                    cout<<" | "<<Muon_PdgId->at(mu[0])<<" | "<<Muon_PdgId->at(mu[1])<<" | "<<Track_pdgId->at(mu[2])<<endl;
                    cout<<" | mu1_MotherpdgID | mu2_MotherpdgId  "<<endl;
                    cout<<" | "<<Muon_MotherPdgId->at(mu[0])<<" | "<<Muon_MotherPdgId->at(mu[1])<<endl;
                }
            }

            bool isSignal = false;
            if(isMC){
                if(std::abs(Muon_PdgId->at(mu[0])) == 13 &&
                   std::abs(Muon_PdgId->at(mu[1])) == 13 &&
                   std::abs(Track_pdgId->at(mu[2])) == 211 &&
                   std::abs(Muon_MotherPdgId->at(mu[0])) == 333 &&
                   std::abs(Muon_MotherPdgId->at(mu[1])) == 333 ) isSignal = true;
            }

            //CUT 0 : Before cuts
            Ncut = 0; cut[Ncut]++; cuttripl[Ncut]++;

            // Fill histograms before selections
            FillHistoBC("MC", j, hMass_tripl_BC, hChi2Vertex, hMassvsChi2, hMass_quad_BC, hMass_quad_Zero_BC, hMass_di_Zero_BC, hMass_di_Zero2_BC, hPtRes_BC, hPtRes_BC_mu, hPtResBarrel_BC, hPtResBarrel_BC_mu, hPtResEndcap_BC, hPtResEndcap_BC_mu, IdsummaryDaughter_Gen, IdsummaryMother_Gen, Idsummary2D_Gen);

            //CUT 1 : L1 fired
            if(!l1double_fired) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<" passed L1 cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //CUT 2 : HLT fired
            if(!hlt_fired) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<" passed HLT cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 3: Chi2 > 0
            if( TripletVtx2_Chi2->at(j) < 0 ) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<" passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 4: all triplets w/ at least 2 track associated with PV
            if(RefittedPV2_NTracks->at(j) < 2) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"   passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 5: check both mu are different from the track
            if( !(DuplicateFinder(Mu01_Eta->at(mu_Ind[0]), Mu02_Eta->at(mu_Ind[1]), Tr_Eta->at(mu_Ind[2]), Mu01_Phi->at(mu_Ind[0]), Mu02_Phi->at(mu_Ind[1]), Tr_Phi->at(mu_Ind[2]), Mu01_Pt->at(mu_Ind[0]), Mu02_Pt->at(mu_Ind[1]), Tr_Pt->at(mu_Ind[2])) ) ) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"     passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 6: check both mu are global
            if( !(Muon_isGlobal->at(mu[0])) || !(Muon_isGlobal->at(mu[1])) ) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"       passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 7: cut on vertex chi2 (Chi2 < 15)
            if( TripletVtx2_Chi2->at(j) >= 15 ) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"         passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 8: muons have opposite charge and mass within 1..1.04 GeV
            double dimass = DimuonMass(MuonCharge->at(mu[0]), MuonCharge->at(mu[1]), Mu01_Pt->at(mu_Ind[0]), Mu02_Pt->at(mu_Ind[1]), Mu01_Eta->at(mu_Ind[0]), Mu02_Eta->at(mu_Ind[1]), Mu01_Phi->at(mu_Ind[0]), Mu02_Phi->at(mu_Ind[1]));
            if(dimass < 1 || dimass > 1.04) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"           passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 9: condition on track Impact Parameter (dz<20 and dxy<0.3)
            if(Track_dz->at(mu[2]) >= 20 || Track_dxy->at(mu[2]) >= 0.3) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"             passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            //trigIndex contains the 3 indeces for the trigger objects with minimum deltaR
            std::vector< std::size_t > trigIndex_deltaR;
            if(isMiniAOD && datasetName.Contains("2017") && MuonPt_HLT2017->size()>0) trigIndex_deltaR = trigMatchDeltaR(j, Muon_HLT2017, false);

            Float_t dR1 = 999;
            Float_t dR2 = 999;
            Float_t dR3 = 999;

            if(isMiniAOD && datasetName.Contains("2017") && MuonPt_HLT2017->size()>0) {
                dR1 = dR( Mu01_Eta->at(j), Muon_HLT2017[trigIndex_deltaR[0]][1], Mu01_Phi->at(j), Muon_HLT2017[trigIndex_deltaR[0]][2]);
                dR2 = dR( Mu02_Eta->at(j), Muon_HLT2017[trigIndex_deltaR[1]][1], Mu02_Phi->at(j), Muon_HLT2017[trigIndex_deltaR[1]][2]);
                dR3 = dR( Tr_Eta->at(j), Muon_HLT2017[trigIndex_deltaR[2]][1], Tr_Phi->at(j), Muon_HLT2017[trigIndex_deltaR[2]][2]);
            }

            Float_t dP1 = 999;
            Float_t dP2 = 999;
            Float_t dP3 = 999;

            if(isMiniAOD && datasetName.Contains("2017") && MuonPt_HLT2017->size()>0) {
                dP1 = std::abs(Muon_HLT2017[trigIndex_deltaR[0]][0] - MuonPt->at(mu[0]))/MuonPt->at(mu[0]);
                dP2 = std::abs(Muon_HLT2017[trigIndex_deltaR[1]][0] - MuonPt->at(mu[1]))/MuonPt->at(mu[1]);
                dP3 = std::abs(Muon_HLT2017[trigIndex_deltaR[2]][0] - Track_pt->at(mu[2]))/Track_pt->at(mu[2]);
            }

            if (isVerbose) {
                cout<<j<<" dR1="<<dR1<<" deltaP/P1="<<dP1<<endl;
                cout<<j<<" dR2="<<dR2<<" deltaP/P2="<<dP2<<endl;
                cout<<j<<" dR3="<<dR3<<" deltaP/P3="<<dP3<<endl;
            }

            // CUT 10: Mu01 dR and dP Trigger Matching
            if( !(dR1<0.03 && dP1<0.1) ) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"               passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 11: Mu02 dR and dP Trigger Matching
            if( !(dR2<0.03 && dP2<0.1) ) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"                 passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            // CUT 12: Track dR and dP Trigger Matching
            if( !(dR3<0.03 && dP3<0.1) ) continue;
            Ncut++; cut[Ncut]++; cuttripl[Ncut]++;
            if (isVerbose) cout<<j<<"                   passed cut "<<Ncut<<endl;
            FillHistoStepByStep(isMC, j, mu_Ind, mu, Ncut, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            ntripl++; triplIndex.push_back(j);
        } // end loop on triplets
        // N. events that passed each selection
        for (int k=0; k<NCUTS; k++){
            if(cuttripl[k] > 0) cutevt[k]++;
        }

        if(ntripl > 0 && isTrigger_forAna) {
            // Histo N. triplets passed for each event
            hNtripl->Fill(ntripl);
            //Best triplet selected based on smaller Chi2    
            ind = BestTripletFinder(triplIndex);
	    double dimass = DimuonMass(MuonCharge->at(mu[0]), MuonCharge->at(mu[1]), Mu01_Pt->at(mu_Ind[0]), Mu02_Pt->at(mu_Ind[1]), Mu01_Eta->at(mu_Ind[0]), Mu02_Eta->at(mu_Ind[1]), Mu01_Phi->at(mu_Ind[0]), Mu02_Phi->at(mu_Ind[1]));
            //RiMatching between index of single mu of the triplet (mu#_Ind) & that of  'MUONID' (mu#) & Ricomputing the 3 possible dimuon masses
            MatchIndex("ID", ind, mu_Ind, mu);
            if(isVerbose) cout<<"Triplets passing presel = "<<ntripl<<" Index of best triplet is "<<ind<<endl;

            //muon scale factors                
            sfMuon scaleFactor;
            scaleFactor = GetMuonSF(SF_h, Mu01_Pt->at(mu_Ind[0]), Mu01_Eta->at(mu_Ind[0]));
            if(isVerbose) cout<<"scaleFactor "<<scaleFactor.value<<" +- "<<scaleFactor.error<<endl;
            scaleFactor = GetMuonSF(SF_h, Mu02_Pt->at(mu_Ind[1]), Mu02_Eta->at(mu_Ind[1]));
            if(isVerbose) cout<<"scaleFactor "<<scaleFactor.value<<" +- "<<scaleFactor.error<<endl;
            scaleFactor = GetMuonSF(SF_h, Tr_Pt->at(mu_Ind[2]), Tr_Eta->at(mu_Ind[2]));
            if(isVerbose) cout<<"scaleFactor "<<scaleFactor.value<<" +- "<<scaleFactor.error<<endl;

            bool isSB_tripletMass = 0;
            if( Triplet2_Mass->at(ind) >= 1.7 && Triplet2_Mass->at(ind) <= 1.8 ) isSB_tripletMass = 1;
            bool isDs_tripletMass = 0;
            if( Triplet2_Mass->at(ind) >= 1.93 && Triplet2_Mass->at(ind) <= 2.01 ) isDs_tripletMass = 1;
    
            TreeFin_Fill(tree, ind, mu_Ind, mu, pileupFactor, Pmu3, cLP, tKink, segmComp, tripletMass, tripletMassReso, fv_nC, fv_dphi3D, fv_d3D, fv_d3Dsig, bs_sv_d3Dsig, bs_sv_d3D, pv_sv_dxy_sig, pv_sv_dxy, d0, d0sig, mindca_iso, trkRel, Pmu1, Ptmu1, etamu1, Pmu2, Ptmu2, etamu2, Ptmu3, etamu3, P_trip, Pt_trip, eta_trip, nStationsMu1, nStationsMu2, nStationsMu3, Iso03Mu1, Iso03Mu2, Iso03Mu3, Iso05Mu1, Iso05Mu2, Iso05Mu3, nMatchesMu1, nMatchesMu2, nMatchesMu3, timeAtIpInOutMu1, timeAtIpInOutMu2, timeAtIpInOutMu3, cQ_uS, cQ_tK, cQ_gK, cQ_tRChi2, cQ_sRChi2, cQ_Chi2LM, cQ_Chi2lD, cQ_gDEP, cQ_tM, cQ_gTP, calEn_emMu1, calEn_emMu2, calEn_emMu3, calEn_hadMu1, calEn_hadMu2, calEn_hadMu3, caloComp, fliDistPVSV_Chi2, isGlb3, isTracker3, isLoose3,  isSoft3, isPF3, isRPC3, isSA3, isCalo3, vx1, vx2, vx3, vy1, vy2, vy3, vz1, vz2, vz3, Refvx1, Refvx2, Refvx3, Refvy1, Refvy2, Refvy3, Refvz1, Refvz2, Refvz3, SVx, SVy, SVz, had03, had05, nJets03, nJets05, nTracks03, nTracks05, sumPt03, sumPt05, hadVeto03, hadVeto05, emVeto03, emVeto05, trVeto03, trVeto05);
            //Muon MiniTree
            //TreeMuon_Fill(tree_mu1, mu[0], run_n, lumi_n, evt_n, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPC, mu_isStandAlone, mu_isTracker, mu_isCalo, mu_isQualityValid, mu_SoftMVA, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_Numberofvalidtrackerhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_validFraction, mu_innerTrack_highPurity, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_IP3D_BS, mu_IP2D_BS, mu_IP3D_PV, mu_IP2D_PV, mu_validMuonHitComb,  mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05);
            //TreeMuon_Fill(tree_mu2, mu[1], run_n, lumi_n, evt_n, mu_pt, mu_eta, mu_phi, mu_energy, mu_charge, mu_isGlobal, mu_isSoft, mu_isLoose, mu_isTight, mu_isPF, mu_isRPC, mu_isStandAlone, mu_isTracker, mu_isCalo, mu_isQualityValid, mu_SoftMVA, mu_isTimeValid, mu_isIsolationValid, mu_numberOfMatchedStations, mu_numberOfMatches, mu_timeAtIpInOut, mu_timeAtIpInOutErr, mu_GLnormChi2, mu_GLhitPattern_numberOfValidMuonHits, mu_trackerLayersWithMeasurement, mu_Numberofvalidpixelhits, mu_Numberofvalidtrackerhits, mu_outerTrack_p, mu_outerTrack_eta, mu_outerTrack_phi, mu_outerTrack_normalizedChi2, mu_outerTrack_muonStationsWithValidHits, mu_innerTrack_p, mu_innerTrack_eta, mu_innerTrack_phi, mu_innerTrack_validFraction, mu_innerTrack_highPurity, mu_innerTrack_normalizedChi2, mu_QInnerOuter, mu_combinedQuality_updatedSta, mu_combinedQuality_trkKink, mu_combinedQuality_glbKink, mu_combinedQuality_trkRelChi2, mu_combinedQuality_staRelChi2, mu_combinedQuality_chi2LocalPosition, mu_combinedQuality_chi2LocalMomentum, mu_combinedQuality_localDistance, mu_combinedQuality_globalDeltaEtaPhi, mu_combinedQuality_tightMatch, mu_combinedQuality_glbTrackProbability, mu_IP3D_BS, mu_IP2D_BS, mu_IP3D_PV, mu_IP2D_PV, mu_validMuonHitComb,  mu_calEnergy_em, mu_calEnergy_emS9, mu_calEnergy_emS25, mu_calEnergy_had, mu_calEnergy_hadS9, mu_segmentCompatibility, mu_caloCompatibility, mu_ptErrOverPt, mu_BestTrackPt, mu_BestTrackPtErr, mu_BestTrackEta, mu_BestTrackEtaErr, mu_BestTrackPhi, mu_BestTrackPhiErr, mu_emEt03, mu_hadEt03, mu_nJets03, mu_nTracks03, mu_sumPt03, mu_hadVetoEt03, mu_emVetoEt03, mu_trackerVetoPt03, mu_emEt05, mu_hadEt05, mu_nJets05, mu_nTracks05, mu_sumPt05, mu_hadVetoEt05, mu_emVetoEt05, mu_trackerVetoPt05);              

            Compute_DimuonMass(mu_Ind, mu, dimu);

            if(isDs_tripletMass) {
                //CUT 13 : final plot Phi mass region
                cutevt[NCUTS-3]++; cut[NCUTS-3]++;
                FillHistoStepByStep(isMC, ind, mu_Ind, mu, NCUTS-3, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
            }
            else if(isSB_tripletMass){
                //CUT 14 : final plot sideband
                cutevt[NCUTS-2]++; cut[NCUTS-2]++;
                FillHistoStepByStep(isMC, ind, mu_Ind, mu, NCUTS-2, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);
            }
            //CUT 15 : final plot full invariant mass range
            cutevt[NCUTS-1]++; cut[NCUTS-1]++;
            FillHistoStepByStep(isMC, ind, mu_Ind, mu, NCUTS-1, hPt, hPt_mu, hEta, hEta_mu, hPhi, hVx, hVy, hVz, hPt_Tr, hEta_Tr, hPt_tripl, hEta_tripl, hPhi_tripl, hMass_tripl, IdsummaryDaughter, IdsummaryMother, Idsummary2D);

            FillHistoAC(ind, mu, hMinSegmComp, hChi2Track, hNMatchedStat, hFlightDist, hFlightDist_Signif, hFlightDistvsP, hPtErrOverPt, hmassdi, dimu, hmassQuad, hmassQuad_Zero, hChi2VertexNorm, hSegmComp, hDeltaR, hTrIPSign);
            hPileUp_AC->Fill(nPileUpInt);
            hNPrVert_AC->Fill(PVCollection_Size, pileupFactor);
            hPropDecayL_AC->Fill((Triplet2_Mass->at(ind)*FlightDistPVSV2->at(ind))/MuonP(Triplet2_Pt->at(ind), Triplet2_Eta->at(ind), Triplet2_Phi->at(ind)), pileupFactor);
            // Trigger requirements
            //TriggerRequirements(ind, hTripTriggerMatched);
        }
        if (ientry < 0) break;
    }//end loop on events
    //Print general info
    cout << endl;
    cout << "TOTAL N. EVENTS -> " << cutevt[0] << endl;
    cout << "TOTAL N. TRIPLETS -> " << cut[0] << endl;
    cout << "Triplets survived: " << cutevt[NCUTS-1] << " || Good: " << NgoodTripl << " , Bad: " << NbadTripl << endl;
    //Histo of cuts Efficiency
    TCanvas *canvEvt = new TCanvas("CutEfficiency_Nevents", "CutEfficiency_Nevents", 0, 0, 1200, 1000);
    Draw_CutEffCanvas(canvEvt, hCutEffEvt, cutevt, listCut);
    TCanvas *canv = new TCanvas("CutEfficiency_Ntriplets", "CutEfficiency_Ntriplets", 0, 0, 1200, 1000);
    Draw_CutEffCanvas(canv, hCutEff, cut, listCut);
    ////PdgId histos
    //dirStepByStep->cd(); dirPdgId->cd();
    //Draw_PdgIdCanvas_StepByStep(PdgIdCanvas_cut, hPdgId_cut, IdsummaryDaughter, PdgIdMotherCanvas_cut, hMotherPdgId_cut, IdsummaryMother, PdgIdCanvas2D_cut, hPdgId2D_cut, Idsummary2D, pId);
    //dirStepByStep->cd(); fout->cd();
    ////PdgId histos Gen BC
    //dirBeforeCuts->cd();
    //Draw_PdgIdCanvasGen(PdgIdCanvas_Gen, hPdgId_Gen, IdsummaryDaughter_Gen, PdgIdMotherCanvas_Gen, hMotherPdgId_Gen, IdsummaryMother_Gen, PdgIdCanvas2D_Gen, hPdgId2D_Gen, Idsummary2D_Gen, pId);
    //fout->cd();
    //Write and close the file
    fout->Write();
    fout->Close();
    
}

// #########################################  END DsPhiPi ANALYSIS CUTFLOW
