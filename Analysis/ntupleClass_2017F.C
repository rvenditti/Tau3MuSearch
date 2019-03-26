#define ntupleClass_2017F_cxx
#include "ntupleClass_2017F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// #########################################    SIGNAL ANALYSIS NEW CUTFLOW

// Cuts: (over triplets)
// * cut[1] -> Chi2 triplet vertex (in 0 - 15)
// * cut[2] -> There are (2 mu glb w/ pt>ptmin (=2) e 1 mu Tracker con pt>ptminTrack (=0.5)) & |eta|<2.4
// * cut[3] -> Triplet mass (in 1.62 - 2 GeV)
// * cut[4] -> The 3 possible pairs of mu of the triplet have proper DeltaR (<0.8)
// * cut[5] -> The 3 possible pairs of mu of the triplet have proper |DeltaZ| (<0.5)
// * cut[6] -> Cut on the dimuon mass w.r.t. Phi(1020) per pairs of mu of the triplet w/ opposite sign

void ntupleClass_2017F::Loop_New(){
    
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    // Variables definition
    int ntripl, cut[7] = {0}, cutevt[7] = {0};
    float ptmin = 2.0, ptminTrack = 0.5, DeltaRmax = 0.8, DeltaZmax = 0.5, DeltaZ1 = 0, DeltaZ2 = 0, DeltaZ3 = 0, ntripltot = 0;
    double massmin = 1.62, massmax = 2.00, sigma = 0.011, mumass = 0.1056583715; // Mass values in GeV
    double TripletVtx_Chi2max = 15, EtaMax = 2.4;
    TString listCut[7];
    //Variables inizialization
    cutevt[0] = nentries;
    cutName_New(listCut);
    // Creation of the file with analysed data
    TFile *fout = new TFile("AnalysedTree_Data2017F_final.root", "RECREATE");
    fout->cd();
    // Creation of histograms for variables BEFORE cuts
    TDirectory *dirBeforeCuts = fout->mkdir("BeforeCuts");
    dirBeforeCuts->cd();
    TH1D *hmasstriBeforeCuts = new TH1D("TripletMass_BC", "TripletMass_BC", 1000, -0.05, 49.95); // binning 50 MeV
    TH1D *hChi2Vertex = new TH1D("Chi2Vertex_BC", "Chi2Vertex_BC", 40, -1.5, 18.5); // binning 0.5
    TH2D *hMassvsChi2 = new TH2D("MassvsChi2_BC", "MassvsChi2_BC", 200, 1., 2.5, 200, -5., 20.);
    TH1F *hmassQuad_Other = new TH1F("QuadMuonMass_BC", "QuadMuonMass_BC", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Other_Zero = new TH1F("QuadMuonMass_Zero_BC", "QuadMuonMass_Zero_BC", 400, -0.05, 79.95); // binning 200 MeV
    TH1D *hmassDi_Other_Zero = new TH1D("DiMuon_2glbMu_Zero_BC", "DiMuon_2glbMu_Zero_BC", 500, -0.05, 24.95); // binning 50 MeV
    TH1D *hmassDi_Other2_Zero = new TH1D("DiMuon_Other", "DiMuon_Other", 500, -0.05, 24.95); // binning 50 MeV
    fout->cd();
    // Creation of INTERMEDIATE histograms
    TDirectory *dirIntermediate = fout->mkdir("Intermediate");
    dirIntermediate->cd();
    TH1D *hmassdi_int = new TH1D("DimuonMass_int", "DimuonMass_int", 500, -0.05, 24.95); // binning 50 MeV
    TH1F *hmassQuad_int = new TH1F("QuadMuonMass_int", "QuadMuonMass_int", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Zero_int = new TH1F("QuadMuonMass_Zero_int", "QuadMuonMass_Zero_int", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Due_int = new TH1F("QuadMuonMass_Due_int", "QuadMuonMass_Due_int", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_MenoDue_int = new TH1F("QuadMuonMass_MenoDue_int", "QuadMuonMass_MenoDue_int", 400, -0.05, 79.95); // binning 200 MeV
    fout->cd();
    // Creation of histograms for variables AFTER cuts
    TDirectory *dirAfterCuts = fout->mkdir("AfterCuts");
    dirAfterCuts->cd();
    TH1I *hNtripl = new TH1I("Ntripl_AC", "Ntripl_AC", 5, -0.5, 4.5);
    TH1D *hPt = new TH1D("Pt_AC", "Pt_AC", 500, -0.05, 24.95); // binning 50 MeV
    TH1D *hEta = new TH1D("Eta_AC", "Eta_AC", 100, -2.5, 2.5); // binning 0.05
    TH1D *hPhi = new TH1D("Phi_AC", "Phi_AC", 140, -3.5, 3.5); // binning 0.05
    TH1D *hIsolation_03 = new TH1D("Isolation03_AC", "Isolation03_AC", 30, -0.5, 29.5); // binning di 1
    TH1D *hIsolation_05 = new TH1D("Isolation05_AC", "Isolation05_AC", 30, -0.5, 29.5); // binning di 1
    TH1F *hChi2Track = new TH1F("Chi2Track", "Chi2Track", 27, -0.3, 5.1); //binning di 0.2
    TH1D *hmasstri = new TH1D("TripletMass_AC", "TripletMass_AC", 42, 1.60, 2.02); // binning 10 MeV
    TH1D *hmassdi = new TH1D("DimuonMass_AC", "DimuonMass_AC", 500, -0.05, 24.95); //binning 50 MeV
    TH1F *hmassQuad = new TH1F("QuadMuonMass_AC", "QuadMuonMass_AC", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Zero = new TH1F("QuadMuonMass_Zero", "QuadMuonMass_Zero", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_Due = new TH1F("QuadMuonMass_Due", "QuadMuonMass_Due", 400, -0.05, 79.95); // binning 200 MeV
    TH1F *hmassQuad_MenoDue = new TH1F("QuadMuonMass_MenoDue", "QuadMuonMass_MenoDue", 400, -0.05, 79.95); // binning 200 MeV
    TH1D *hNMatchedStat = new TH1D("NofMatchedStations", "NofMatchedStations", 6, -0.5, 5.5);
    TH1D *hFlightDist = new TH1D("FlightDist_AC", "FlightDist_AC", 90, 0., 3.);
    TH1D *hFlightDist_Signif = new TH1D("FlightDist_Signif_AC", "FlightDist_Signif_AC", 100, 0., 100.);
    TH2D *hFlightDistvsP = new TH2D("FlightDistvsP_AC", "FlightDistvsP_AC", 20, 0., 0.3, 20, 0., 45.);
    fout->cd();
    // Creation of total histograms
    TH1I *hCutEff = new TH1I("CutEff_Ntriplets", "CutEff_Ntriplets", 7, 0.5, 7.5);
    TH1I *hCutEffEvt = new TH1I("CutEff_NEvents", "CutEff_NEvents", 7, 0.5, 7.5);
    
    //Loop over the events
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        if (jentry == 2500000) cout << "We are at 2.5 M" << endl;
        ntripl = 0; int cutevt2[7] = {0};
        Long64_t ientry = fChain->LoadTree(jentry);
        fChain->GetEntry(ientry);
        //Loop over the TRIPLETS
        for (int j=0; j<TripletVtx_Chi2->size(); j++){
            //Check sul numero di tracce nel vertice primario
            if(RefittedPV_NTracks->at(j) > 1){
                ntripltot++;
                
//                // Esperimento ... (4 muons inv mass)
//                if(MuonPt->size() > 3){
//                    for (int k=0; k<(MuonPt->size()-3); k++){
//                        for (int l=0; l <(MuonPt->size()-2); l++){
//                            for (int m=0; m<(MuonPt->size()-1); m++){
//                                for (int g=0; g<(MuonPt->size()); g++){
//                                    float quadMuMass = QuadMuonMass(MuonPt->at(k), MuonPt->at(l), MuonPt->at(m), MuonPt->at(g), MuonEta->at(k), MuonEta->at(l), MuonEta->at(m), MuonEta->at(g), MuonPhi->at(k), MuonPhi->at(l), MuonPhi->at(m), MuonPhi->at(g), mumass);
//                                    hmassQuad_Other->Fill(quadMuMass);
//                                    if((MuonCharge->at(k) + MuonCharge->at(l) + MuonCharge->at(m) + MuonCharge->at(g)) == 0)
//                                        hmassQuad_Other_Zero->Fill(quadMuMass);
//                                }
//                            }
//                        }
//                    }
//                }
                
                // 2 GLB muons inv mass
                if(MuonPt->size() > 1){
                    for (int k=0; k<(MuonPt->size()-1); k++){
                        for (int l=0; l<(MuonPt->size()); l++){
                            if(Muon_isGlobal->at(k) == 1 && Muon_isGlobal->at(l) == 1 && MuonPt->at(k) > ptmin && MuonPt->at(l) > ptmin){
                                double dimass = DimuonMass(MuonCharge->at(k), MuonCharge->at(l), MuonPt->at(k), MuonPt->at(l), MuonEta->at(k), MuonEta->at(l), MuonPhi->at(k), MuonPhi->at(l), mumass);
                                if(dimass != 0) hmassDi_Other_Zero->Fill(dimass);
                                if(MuonPt->at(k) >= 5 && MuonPt->at(l) >= 5){
                                    if(dimass != 0) hmassDi_Other2_Zero->Fill(dimass);
                                }
                            }
                        }
                    }
                }
                
                
                //Matching tra gli indici dei singoli muoni dei tripletto (mu#_Ind) e quelli 'generali' (mu#)
                int mu1_Ind = Mu1_TripletIndex->at(j);
                int mu2_Ind = Mu2_TripletIndex->at(j);
                int mu3_Ind = Mu3_TripletIndex->at(j);
                if (mu1_Ind != j || mu2_Ind != j || mu3_Ind != j) cout << "Indici diversi ?!?" << endl;
                int mu1 = MuonFinder(Mu1_Pt->at(mu1_Ind), Mu1_Eta->at(mu1_Ind), Mu1_Phi->at(mu1_Ind));
                int mu2 = MuonFinder(Mu2_Pt->at(mu2_Ind), Mu2_Eta->at(mu2_Ind), Mu2_Phi->at(mu2_Ind));
                int mu3 = MuonFinder(Mu3_Pt->at(mu3_Ind), Mu3_Eta->at(mu3_Ind), Mu3_Phi->at(mu3_Ind));
                // BEFORE cuts
                //Riempio alcuni isto del tripletto
                hmasstriBeforeCuts->Fill(Triplet_Mass->at(j));
                hChi2Vertex->Fill(TripletVtx_Chi2->at(j));
                hMassvsChi2->Fill(Triplet_Mass->at(j), TripletVtx_Chi2->at(j));
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
                                for(int k=0; k<MuonPt->size(); k++){
                                    if(k != mu1 && k!= mu2 && k!= mu3){
                                        if(Muon_isTrackerMuon->at(k) == 1 && MuonPt->at(k) > ptminTrack){
                                            float quadMass = QuadMuonMass(MuonPt->at(mu1), MuonPt->at(mu2), MuonPt->at(mu3), MuonPt->at(k), MuonEta->at(mu1), MuonEta->at(mu2), MuonEta->at(mu3), MuonEta->at(k), MuonPhi->at(mu1), MuonPhi->at(mu2), MuonPhi->at(mu3), MuonPhi->at(k), mumass);
                                            hmassQuad_int->Fill(quadMass);
                                            if((MuonCharge->at(mu1)+MuonCharge->at(mu2)+MuonCharge->at(mu3)+MuonCharge->at(k)) == 0){
                                                hmassQuad_Zero_int->Fill(quadMass);
                                            }
                                            if((MuonCharge->at(mu1)+MuonCharge->at(mu2)+MuonCharge->at(mu3)+MuonCharge->at(k)) == 2)
                                                hmassQuad_Due_int->Fill(quadMass);
                                            if((MuonCharge->at(mu1)+MuonCharge->at(mu2)+MuonCharge->at(mu3)+MuonCharge->at(k)) == -2)
                                                hmassQuad_MenoDue_int->Fill(quadMass);
                                        }
                                    }
                                }
                                //check condition on trimuon mass
                                if(Triplet_Mass->at(j) > massmin && Triplet_Mass->at(j) < massmax){
                                    cut[3]++; cutevt2[3]++;
                                    // Loop sulle COPPIE di muoni del tripletto e check su DeltaR e |Delta Z|
                                    if(isPairDeltaRGood(j, DeltaRmax) == true){
                                        cut[4]++; cutevt2[4]++;
                                        if(isPairDeltaZGood(DeltaZ1, DeltaZ2, DeltaZ3, DeltaZmax) == true){
                                            cut[5]++; cutevt2[5]++;
                                            // Controllo se sono compatibili con la Phi(1020)
                                            if(isPairNotAPhi(dimu1_2, dimu2_3, dimu1_3, sigma) == true){
                                                cut[6]++; ntripl++;
                                                double TripletP = MuonP(Triplet_Pt->at(j), Triplet_Eta->at(j), Triplet_Phi->at(j), mumass);
                                                // Riempio isto con caratteristiche dei 3 muoni che hanno superato tutti i tagli
                                                hPt->Fill(Mu1_Pt->at(mu1_Ind));
                                                hPt->Fill(Mu2_Pt->at(mu2_Ind));
                                                hPt->Fill(Mu3_Pt->at(mu3_Ind));
                                                hEta->Fill(Mu1_Eta->at(mu1_Ind));
                                                hEta->Fill(Mu2_Eta->at(mu2_Ind));
                                                hEta->Fill(Mu3_Eta->at(mu3_Ind));
                                                hPhi->Fill(Mu1_Phi->at(mu1_Ind));
                                                hPhi->Fill(Mu2_Phi->at(mu2_Ind));
                                                hPhi->Fill(Mu3_Phi->at(mu3_Ind));
                                                hIsolation_03->Fill(Muon_emEt03->at(mu1));
                                                hIsolation_03->Fill(Muon_emEt03->at(mu2));
                                                hIsolation_03->Fill(Muon_emEt03->at(mu3));
                                                hIsolation_05->Fill(Muon_emEt03->at(mu1));
                                                hIsolation_05->Fill(Muon_emEt03->at(mu2));
                                                hIsolation_05->Fill(Muon_emEt03->at(mu3));
                                                hChi2Track->Fill(Muon_innerTrack_normalizedChi2->at(mu1));
                                                hChi2Track->Fill(Muon_innerTrack_normalizedChi2->at(mu2));
                                                hChi2Track->Fill(Muon_innerTrack_normalizedChi2->at(mu3));
                                                hNMatchedStat->Fill(Muon_numberOfMatchedStations->at(mu1));
                                                hNMatchedStat->Fill(Muon_numberOfMatchedStations->at(mu2));
                                                hNMatchedStat->Fill(Muon_numberOfMatchedStations->at(mu3));
                                                hmasstri->Fill(Triplet_Mass->at(j));
                                                hFlightDist->Fill(FlightDistPVSV->at(j));
                                                hFlightDist_Signif->Fill(FlightDistPVSV_Significance->at(j));
                                                hFlightDistvsP->Fill(FlightDistPVSV->at(j), TripletP);
                                                DiMuonHisto(hmassdi, dimu1_2, dimu2_3, dimu1_3);
                                                //Aggiungo al tripletto una 4a traccia e ne faccio la massa inv
                                                for(int k=0; k<MuonPt->size(); k++){
                                                    if(k != mu1 && k!= mu2 && k!= mu3){
                                                        if(Muon_isTrackerMuon->at(k) == 1 && MuonPt->at(k) > ptminTrack){
                                                            float quadMass = QuadMuonMass(MuonPt->at(mu1), MuonPt->at(mu2), MuonPt->at(mu3), MuonPt->at(k), MuonEta->at(mu1), MuonEta->at(mu2), MuonEta->at(mu3), MuonEta->at(k), MuonPhi->at(mu1), MuonPhi->at(mu2), MuonPhi->at(mu3), MuonPhi->at(k), mumass);
                                                            hmassQuad->Fill(quadMass);
                                                            if((MuonCharge->at(mu1)+MuonCharge->at(mu2)+MuonCharge->at(mu3)+MuonCharge->at(k)) == 0){
                                                                hmassQuad_Zero->Fill(quadMass);
                                                            }
                                                            if((MuonCharge->at(mu1)+MuonCharge->at(mu2)+MuonCharge->at(mu3)+MuonCharge->at(k)) == 2)
                                                                hmassQuad_Due->Fill(quadMass);
                                                            if((MuonCharge->at(mu1)+MuonCharge->at(mu2)+MuonCharge->at(mu3)+MuonCharge->at(k)) == -2)
                                                                hmassQuad_MenoDue->Fill(quadMass);
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
                }
            }
        } // end loop sui tripletti
        // Conto numero di eventi passati dopo ogni taglio
        for (int k=1; k<7; k++){
            if(cutevt2[k] > 0) cutevt[k]++;
        }
        // Riempio isto con il numero di tripletti passati per evento
        hNtripl->Fill(ntripl);
        if(ntripl > 0) cutevt[6]++;
        if (ientry < 0) break;
    }//end loop sugli eventi
    cut[0] = ntripltot;
    //Stampo info generali & disegno gli istogrammi
    cout << endl;
    cout << "NUMERO EVENTI -> " << nentries << endl << endl;
    cout << "NUMERO TRIPLETTI -> " << ntripltot << endl << endl;
    //Histo AFTER cuts
    hNtripl->DrawCopy();
    hPt->DrawCopy();
    hEta->DrawCopy();
    hPhi->DrawCopy();
    hChi2Track->DrawCopy();
    hNMatchedStat->DrawCopy();
    hNMatchedStat->DrawCopy("HIST TEXT0");
    hmasstri->DrawCopy();
    hmassdi->DrawCopy();
    hmassQuad->DrawCopy();
    hmassQuad_Zero->DrawCopy();
    hmassQuad_Due->DrawCopy();
    hmassQuad_MenoDue->DrawCopy();
    hFlightDist->DrawCopy();
    hFlightDist_Signif->DrawCopy();
    hFlightDistvsP->GetXaxis()->SetTitle("FlightDist");
    hFlightDistvsP->GetYaxis()->SetTitle("P (GeV)");
    hFlightDistvsP->DrawCopy("colz");
    //Histo BEFORE cuts
    hmasstriBeforeCuts->DrawCopy();
    hChi2Vertex->DrawCopy();
    hMassvsChi2->GetXaxis()->SetTitle("Mass (GeV)");
    hMassvsChi2->GetYaxis()->SetTitle("Chi2triplet");
    hMassvsChi2->DrawCopy("colz");
    hmassQuad_Other->DrawCopy();
    hmassQuad_Other_Zero->DrawCopy();
    hmassDi_Other_Zero->DrawCopy();
    //Histo intermediate
    hmassdi_int->DrawCopy();
    hmassQuad_int->DrawCopy();
    hmassQuad_Zero_int->DrawCopy();
    hmassQuad_Due_int->DrawCopy();
    hmassQuad_MenoDue_int->DrawCopy();
    //Histo of cuts Efficiency
    TCanvas *canvEvt = new TCanvas("CutEfficiency_Nevents", "CutEfficiency_Nevents", 0, 0, 1200, 1000);
    CutEffHisto(canvEvt, hCutEffEvt, cutevt, listCut);
    TCanvas *canv = new TCanvas("CutEfficiency_Ntriplets", "CutEfficiency_Ntriplets", 0, 0, 1200, 1000);
    CutEffHisto(canv, hCutEff, cut, listCut);
    //Write and close the file
    fout->Write();
    fout->Close();
}

// #########################################  END SIGNAL ANALYSIS  NEW CUTFLOW
