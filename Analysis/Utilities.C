
#define ntupleClass_MC_cxx
#define NCUTS 8
#define NPARTICLES 560
#define NMU 3
#define mumass 0.1056583715 // Muon mass in GeV
#define PhiMass 1.019461 // Phi mass in GeV
#define OmegaMass 0.78265 // Omega mass in GeV
#define ptmin 2.0

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


Int_t ntupleClass_MC::BestTripletFinder(Int_t triplIndex[1000], Int_t n){
    // Given the index of all the triplets of an event that passed all the cuts, it returns the index of the one with the smallest Chi2 of the vertex
    int index = 0; double bestChi2 = 15;
    if (n>1000) cout << "There are too many good triplets!" << endl;
    for(int i=0; i<n; i++){
        if(TripletVtx_Chi2->at(triplIndex[i]) < bestChi2){
            bestChi2 = TripletVtx_Chi2->at(triplIndex[i]);
            index = triplIndex[i];
        }
    }
    return index;
}

Double_t ntupleClass_MC::DimuonMass(Double_t charge1, Double_t charge2, Double_t pt1, Double_t pt2, Double_t eta1, Double_t eta2, Double_t phi1, Double_t phi2){
    // Given the characteristics of 2 muons, if their charge is opposite the function returns their invariant mass, otherwise it returns 0
    double inv = 0;
    if(charge1 + charge2 != 0)  return inv;
    else {
        TLorentzVector mu1, mu2, mutot;
        mu1.SetPtEtaPhiM(pt1, eta1, phi1, mumass);
        mu2.SetPtEtaPhiM(pt2, eta2, phi2, mumass);
        mutot = mu1 + mu2;
        return mutot.M();
    }
}

void ntupleClass_MC::Draw_CutEffCanvas(TCanvas *canv, TH1I *hist, Int_t cut[NCUTS], TString listCut[NCUTS]){
    // This function writes on the canvas the histo of the cuts efficiency
    for(int k=0; k<NCUTS; k++){
        hist->Fill(k+1, cut[k]);
        hist->GetXaxis()->SetBinLabel(k+1, listCut[k]);
    }
    canv->SetLogy();
    hist->DrawCopy("HIST TEXT0");
    canv->Write();
    canv->Close();
}

void ntupleClass_MC::Draw_PdgIdCanvas(TCanvas *canv, TH1I *hist, Int_t Idsummary[NPARTICLES], TString pIdList[NPARTICLES]){
    // This function writes on the canvas the histo w/ the particle names
    int l=0;
    for(int k=0; k<NPARTICLES; k++){
        if(Idsummary[k]>0){
            hist->Fill(l, Idsummary[k]);
            hist->GetXaxis()->SetBinLabel(l+1, pIdList[k]);
            l++;
        }
    }
    canv->SetLogy();
    hist->SetMarkerSize(0.7);
    hist->DrawCopy();
    hist->DrawCopy("HIST TEXT0");
    canv->Write();
    canv->Close();
}

void ntupleClass_MC::Draw_PdgIdCanvas_2D(TCanvas *canv, TH2I *hist, Int_t Idsummary2D[NPARTICLES][NPARTICLES], TString pIdList[NPARTICLES]){
    // This function writes on the canvas the histo 2D w/ the particle names
    int l=0, m=0; int fillX = 0, fillY = 0;
    for(int k=0; k<NPARTICLES; k++){ // loop on x
        fillX = 0; m=0;
        for(int i=0; i<NPARTICLES; i++){ // loop on y (check if there is at least 1 element(y) [for x=k] that is !=0)
            if(Idsummary2D[k][i] > 0)
                fillX = 1;
        }
        if(fillX == 1){
            hist->GetXaxis()->SetBinLabel(l+1, pIdList[k]); //put name on x
            for(int i=0; i<NPARTICLES; i++){ // loop on y
                fillY = 0;
                for(int k1=0; k1<NPARTICLES; k1++){ // loop on x (in order to check if there is at least a case in which that y is !=0)
                    if(Idsummary2D[k1][i] > 0)
                        fillY = 1;
                }
                if (fillY==1){
                    if(Idsummary2D[k][i] > 0){
                        hist->GetYaxis()->SetBinLabel(m+1, pIdList[i]);
                        hist->Fill(l, m, Idsummary2D[k][i]);
                    }
                    m++;
                }
            }
            l++;
        }
    }
    gStyle->SetPalette(kBlackBody);
    hist->SetMarkerSize(0.7);
    hist->GetXaxis()->SetTitleOffset(0.8);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->DrawCopy("colz TEXT0");
    canv->Write();
    canv->Close();
}

void ntupleClass_MC::Draw_PdgIdCanvas_StepByStep(TCanvas *PdgIdCanvas_cut[NCUTS], TH1I *hPdgId_cut[NCUTS], Int_t IdsummaryDaughter[NCUTS][NPARTICLES], TCanvas *PdgIdMotherCanvas_cut[NCUTS], TH1I *hMotherPdgId_cut[NCUTS], Int_t IdsummaryMother[NCUTS][NPARTICLES], TCanvas *PdgIdCanvas2D_cut[NCUTS], TH2I *hPdgId2D_cut[NCUTS], Int_t Idsummary2D[NCUTS][NPARTICLES][NPARTICLES], TString pId[NPARTICLES]){
    // This function draws the Pdg histograms
    for (int i=0; i<NCUTS; i++){
        TString canvName = "PdgId_Daughter_cut"; canvName += i;
        TString canvNameMother = "PdgId_Mother_cut"; canvNameMother += i;
        TString canvName2D = "PdgId2D_cut"; canvName2D += i;
        PdgIdCanvas_cut[i] = new TCanvas(canvName, canvName, 0, 0, 1600, 1000);
        Draw_PdgIdCanvas(PdgIdCanvas_cut[i], hPdgId_cut[i], IdsummaryDaughter[i], pId);
        PdgIdMotherCanvas_cut[i] = new TCanvas(canvNameMother, canvNameMother, 0, 0, 1600, 1000);
        Draw_PdgIdCanvas(PdgIdMotherCanvas_cut[i], hMotherPdgId_cut[i], IdsummaryMother[i], pId);
        PdgIdCanvas2D_cut[i] = new TCanvas(canvName2D, canvName2D, 0, 0, 1600, 1000);
        Draw_PdgIdCanvas_2D(PdgIdCanvas2D_cut[i], hPdgId2D_cut[i], Idsummary2D[i], pId);
    }
}

void ntupleClass_MC::Draw_PdgIdCanvasGen(TCanvas *PdgIdCanvas_Gen, TH1I *hPdgId_Gen, Int_t IdsummaryDaughter_Gen[NPARTICLES], TCanvas *PdgIdMotherCanvas_Gen, TH1I *hMotherPdgId_Gen, Int_t IdsummaryMother_Gen[NPARTICLES], TCanvas *PdgIdCanvas2D_Gen, TH2I *hPdgId2D_Gen, Int_t Idsummary2D_Gen[NPARTICLES][NPARTICLES], TString pId[NPARTICLES]){
    // This function draws the Pdg histograms for Gen particles
    PdgIdCanvas_Gen = new TCanvas("PdgId_Daughter_Gen", "PdgId_Daughter_Gen", 0, 0, 1600, 1000);
    Draw_PdgIdCanvas(PdgIdCanvas_Gen, hPdgId_Gen, IdsummaryDaughter_Gen, pId);
    PdgIdMotherCanvas_Gen = new TCanvas("PdgId_Mother_Gen", "PdgId_Mother_Gen", 0, 0, 1600, 1000);
    Draw_PdgIdCanvas(PdgIdMotherCanvas_Gen, hMotherPdgId_Gen, IdsummaryMother_Gen, pId);
    PdgIdCanvas2D_Gen = new TCanvas("PdgId2D_Gen", "PdgId2D_Gen", 0, 0, 1600, 1000);
    Draw_PdgIdCanvas_2D(PdgIdCanvas2D_Gen, hPdgId2D_Gen, Idsummary2D_Gen, pId);
}

void ntupleClass_MC::Fill_CutName(TString listCut[NCUTS]){
    // Init a vector of strings w/ the names of the cuts
    listCut[0] = "BeforeCuts";
    listCut[1] = "#chi^{2} triplet";
    listCut[2] = "2glb+1trk";
    listCut[3] = "MassTriplet";
    listCut[4] = "DeltaR";
    listCut[5] = "DeltaZ";
    listCut[6] = "VETO #Phi mass";
    listCut[7] = "VETO #omega mass";
}

void ntupleClass_MC::Fill_DimuonMass(Int_t mu_Ind[NMU], Int_t mu[NMU], Double_t dimu[NMU]){
    // Fills the vector w/ the 3 possible dimuon masses of the muons of the triplet
    double pt[NMU] = {0}, eta[NMU] = {0}, phi[NMU] = {0};
    Fill_MuonVariables(mu_Ind, pt, eta, phi);
    dimu[0] = DimuonMass(MuonCharge->at(mu[0]), MuonCharge->at(mu[1]), pt[0], pt[1], eta[0], eta[1], phi[0], phi[1]); // dimuon mass 1-2
    dimu[1] = DimuonMass(MuonCharge->at(mu[1]), MuonCharge->at(mu[2]), pt[1], pt[2], eta[1], eta[2], phi[1], phi[2]); // dimuon mass 2-3
    dimu[2] = DimuonMass(MuonCharge->at(mu[0]), MuonCharge->at(mu[2]), pt[0], pt[2], eta[0], eta[2], phi[0], phi[2]); // dimuon mass 1-3
}

void ntupleClass_MC::Fill_MuonVariables(Int_t mu_Ind[NMU], Double_t pt[NMU], Double_t eta[NMU], Double_t phi[NMU]){
    // Fills vectors w/ the variables of the muons of the triplet
    pt[0] = Mu1_Pt->at(mu_Ind[0]);
    pt[1] = Mu2_Pt->at(mu_Ind[1]);
    pt[2] = Mu3_Pt->at(mu_Ind[2]);
    eta[0] = Mu1_Eta->at(mu_Ind[0]);
    eta[1] = Mu2_Eta->at(mu_Ind[1]);
    eta[2] = Mu3_Eta->at(mu_Ind[2]);
    phi[0] = Mu1_Phi->at(mu_Ind[0]);
    phi[1] = Mu2_Phi->at(mu_Ind[1]);
    phi[2] = Mu3_Phi->at(mu_Ind[2]);
}

void ntupleClass_MC::Fill_MuonVariablesGen(Int_t muGen[NMU], Double_t ptGEN[NMU], Double_t etaGEN[NMU], Double_t phiGEN[NMU]){
    // Fills vectors w/ the variables GEN of the muons of the triplet
    ptGEN[0] = GenMatchMu1_Pt->at(muGen[0]);
    ptGEN[1] = GenMatchMu2_Pt->at(muGen[1]);
    ptGEN[2] = GenMatchMu3_Pt->at(muGen[2]);
    etaGEN[0] = GenMatchMu1_Eta->at(muGen[0]);
    etaGEN[1] = GenMatchMu2_Eta->at(muGen[1]);
    etaGEN[2] = GenMatchMu3_Eta->at(muGen[2]);
    phiGEN[0] = GenMatchMu1_Phi->at(muGen[0]);
    phiGEN[1] = GenMatchMu2_Phi->at(muGen[1]);
    phiGEN[2] = GenMatchMu3_Phi->at(muGen[2]);
}

void ntupleClass_MC::Fill_MuonVariablesGen_Sim(Int_t muGen[NMU], Double_t ptSimGEN[NMU], Double_t etaSimGEN[NMU], Double_t phiSimGEN[NMU]){
    // Fills vectors w/ the variables GEN of the muons of the triplet @gen level
    ptSimGEN[0] = GenMatchMu1_SimPt->at(muGen[0]);
    ptSimGEN[1] = GenMatchMu2_SimPt->at(muGen[1]);
    ptSimGEN[2] = GenMatchMu3_SimPt->at(muGen[2]);
    etaSimGEN[0] = GenMatchMu1_SimEta->at(muGen[0]);
    etaSimGEN[1] = GenMatchMu2_SimEta->at(muGen[1]);
    etaSimGEN[2] = GenMatchMu3_SimEta->at(muGen[2]);
    phiSimGEN[0] = GenMatchMu1_SimPhi->at(muGen[0]);
    phiSimGEN[1] = GenMatchMu2_SimPhi->at(muGen[1]);
    phiSimGEN[2] = GenMatchMu3_SimPhi->at(muGen[2]);
}

void ntupleClass_MC::Fill_ParticleIdSummary(Int_t mu[NMU], Int_t IdsummaryDaughter[NPARTICLES], Int_t IdsummaryMother[NPARTICLES], Int_t IdSummary2D[NPARTICLES][NPARTICLES]){
    // This function fills the IdSummary vector w/ the all the particles (daughter and mother) present
    for(int i=0; i<NMU; i++){
        Fill_particleId(Muon_PdgId->at(mu[i]), IdsummaryDaughter);
        Fill_particleId(Muon_MotherPdgId->at(mu[i]), IdsummaryMother);
        Fill_particleId_2D(Muon_PdgId->at(mu[i]), Muon_MotherPdgId->at(mu[i]), IdSummary2D);
    }
}

void ntupleClass_MC::FillHistoAC(Int_t ind, Int_t mu[NMU], TH1F *hChi2Track, TH1D *hNMatchedStat, TH1D *hFlightDist, TH1D *hFlightDist_Signif, TH2D *hFlightDistvsP, TH1D *hPtErrOverPt, TH1D *hmassdi, Double_t dimu[NMU], TH1F *hmassQuad, TH1F *hmassQuad_Zero){
    // Fills the histograms after cuts
    for(int k=0; k<NMU; k++){
        hChi2Track->Fill(Muon_innerTrack_normalizedChi2->at(mu[k]));
        hNMatchedStat->Fill(Muon_numberOfMatchedStations->at(mu[k]));
        hPtErrOverPt->Fill(Muon_ptErrOverPt->at(mu[k]));
    }
    hFlightDist->Fill(FlightDistPVSV->at(ind));
    hFlightDist_Signif->Fill(FlightDistPVSV_Significance->at(ind));
    double TripletP = MuonP(Triplet_Pt->at(ind), Triplet_Eta->at(ind), Triplet_Phi->at(ind));
    hFlightDistvsP->Fill(FlightDistPVSV->at(ind), TripletP);
    FillHistoDiMuMass_AC(hmassdi, dimu);
    FillHistoQuadMuMass_AC(hmassQuad, hmassQuad_Zero, mu);
}

void ntupleClass_MC::FillHistoBC(TString type, Int_t ind, TH1D *hMass_tripl, TH1D *hChi2Vertex, TH2D *hMassvsChi2, TH1F *hMass_quad, TH1F *hMass_quad_Zero, TH1D *hMass_di, TH1D *hMass_di2, TH1D *hPtRes, TH1D *hPtRes_mu[NMU], TH1D *hPtResBarrel, TH1D *hPtResBarrel_mu[NMU], TH1D *hPtResEndcap, TH1D *hPtResEndcap_mu[NMU], Int_t IdsummaryDaughter_Gen[NPARTICLES], Int_t IdsummaryMother_Gen[NPARTICLES], Int_t Idsummary2D_Gen[NPARTICLES][NPARTICLES]){
    // Fills the histograms before cuts
    hMass_tripl->Fill(Triplet_Mass->at(ind));
    hChi2Vertex->Fill(TripletVtx_Chi2->at(ind));
    hMassvsChi2->Fill(Triplet_Mass->at(ind), TripletVtx_Chi2->at(ind));
    if(MuonPt->size() > NMU) FillHistoQuadMuMass_BC(hMass_quad, hMass_quad_Zero);
    if(MuonPt->size() > 1) FillHistoDiMuMass_BC(hMass_di, hMass_di2);
    if (strcmp(type, "data") != 0){
        FillHistoResoPt_BC(hPtRes, hPtRes_mu, hPtResBarrel, hPtResBarrel_mu, hPtResEndcap, hPtResEndcap_mu);
        for (int k=0; k<GenParticle_MotherPdgId->size(); k++){
            Fill_particleId(GenParticle_PdgId->at(k), IdsummaryDaughter_Gen);
            Fill_particleId(GenParticle_MotherPdgId->at(k), IdsummaryMother_Gen);
            Fill_particleId_2D(GenParticle_PdgId->at(k), GenParticle_MotherPdgId->at(k), Idsummary2D_Gen);
        }
    }
}
void ntupleClass_MC::FillHistoDiMuMass_AC(TH1D *hist, Double_t dimu[NMU]){
    // This function fills the dimuon mass histogram
    for(int i=0; i<NMU; i++){
        if(dimu[i] != 0) hist->Fill(dimu[i]);
    }
}

void ntupleClass_MC::FillHistoDiMuMass_BC(TH1D *h_Zero, TH1D *h_Zero2){
    // It loops over the muons and fills the histo w/ the inv mass of o.s. GLB muons w/ |DeltaZ| < 0.5 and other conditions...
    for (int k=0; k<(MuonPt->size()-1); k++){
        for (int l=k+1; l<(MuonPt->size()); l++){
            if(Muon_isGlobal->at(k) == 1 && Muon_isGlobal->at(l) == 1 && MuonPt->at(k) > ptmin && MuonPt->at(l) > ptmin && abs(Muon_vz->at(k) - Muon_vz->at(l)) < 0.5){
                double dimass = DimuonMass(MuonCharge->at(k), MuonCharge->at(l), MuonPt->at(k), MuonPt->at(l), MuonEta->at(k), MuonEta->at(l), MuonPhi->at(k), MuonPhi->at(l));
                if(dimass != 0) h_Zero->Fill(dimass);
                if(MuonPt->at(k) >= 5 && MuonPt->at(l) >= 5){
                    if(dimass != 0) h_Zero2->Fill(dimass);
                }
            }
        }
    }
}


void ntupleClass_MC::FillHistoQuadMuMass_AC(TH1F *h, TH1F *h_Zero, Int_t mu[NMU]){
    // Computes the invariant mass of 4 muons (3mu of the triplet + another tracker mu) w/ |DeltaZ| < 0.5 and fills the histo
    for(int k=0; k<MuonPt->size(); k++){
        if(k != mu[0] && k!= mu[1] && k!= mu[2]){
            if(Muon_isTrackerMuon->at(k) == 1 && MuonPt->at(k) > 0.5){
                if (isPairDeltaZGood(Muon_vz->at(k), Muon_vz->at(mu[0]), Muon_vz->at(mu[1]), 0.5) == true && isDeltaZGood(Muon_vz->at(k), Muon_vz->at(mu[2]), 0.5) == true){
                    float quadMass = QuadMuonMass(MuonPt->at(mu[0]), MuonPt->at(mu[1]), MuonPt->at(mu[2]), MuonPt->at(k), MuonEta->at(mu[0]), MuonEta->at(mu[1]), MuonEta->at(mu[2]), MuonEta->at(k), MuonPhi->at(mu[0]), MuonPhi->at(mu[1]), MuonPhi->at(mu[2]), MuonPhi->at(k));
                    h->Fill(quadMass);
                    if((MuonCharge->at(mu[0])+MuonCharge->at(mu[1])+MuonCharge->at(mu[2])+MuonCharge->at(k)) == 0)
                        h_Zero->Fill(quadMass);
                }
            }
        }
    }
}

void ntupleClass_MC::FillHistoQuadMuMass_BC(TH1F *h, TH1F *h_Zero){
    // Computes the invariant mass of 4 muons w/ |DeltaZ| < 0.5 and fills the histo
    for (int k=0; k<(MuonPt->size()-3); k++){
        for (int l=k+1; l <(MuonPt->size()-2); l++){
            for (int m=l+1; m<(MuonPt->size()-1); m++){
                for (int g=m+1; g<(MuonPt->size()); g++){
                    if (isPairDeltaZGood(Muon_vz->at(k), Muon_vz->at(l), Muon_vz->at(m), 0.5) == true && isDeltaZGood(Muon_vz->at(m), Muon_vz->at(g), 0.5) == true && isDeltaZGood(Muon_vz->at(k), Muon_vz->at(g), 0.5) == true && isDeltaZGood(Muon_vz->at(l), Muon_vz->at(g), 0.5) == true){
                        float quadMuMass = QuadMuonMass(MuonPt->at(k), MuonPt->at(l), MuonPt->at(m), MuonPt->at(g), MuonEta->at(k), MuonEta->at(l), MuonEta->at(m), MuonEta->at(g), MuonPhi->at(k), MuonPhi->at(l), MuonPhi->at(m), MuonPhi->at(g));
                        h->Fill(quadMuMass);
                        if((MuonCharge->at(k) + MuonCharge->at(l) + MuonCharge->at(m) + MuonCharge->at(g)) == 0)
                            h_Zero->Fill(quadMuMass);
                    }
                }
            }
        }
    }
}

void ntupleClass_MC::FillHistoResoPt_AC(Int_t muGen[NMU], TH1D *hPtRes, TH1D *hPtRes_mu[NMU], TH1D *hPtResBarrel, TH1D *hPtResBarrel_mu[NMU], TH1D *hPtResEndcap, TH1D *hPtResEndcap_mu[NMU]){
    // Pt Reso After cuts
    double ptResMu[NMU] = {0}, ptGEN[NMU] = {0}, etaGEN[NMU] = {0}, phiGEN[NMU] = {0}, ptSimGEN[NMU] = {0}, etaSimGEN[NMU] = {0}, phiSimGEN[NMU] = {0};
    Fill_MuonVariablesGen(muGen, ptGEN, etaGEN, phiGEN);
    Fill_MuonVariablesGen_Sim(muGen, ptSimGEN, etaSimGEN, phiSimGEN);
    for(int k=0; k<NMU; k++){
        ptResMu[k] = (ptSimGEN[k] - ptGEN[k])/ptSimGEN[k];
        hPtRes_mu[k]->Fill(ptResMu[k]);
        hPtRes->Fill(ptResMu[k]);
        if(etaGEN[k] < 1.4){
            hPtResBarrel_mu[k]->Fill(ptResMu[k]);
            hPtResBarrel->Fill(ptResMu[k]);
        }
        else {
            hPtResEndcap_mu[k]->Fill(ptResMu[k]);
            hPtResEndcap->Fill(ptResMu[k]);
        }
    }
}

void ntupleClass_MC::FillHistoResoPt_BC(TH1D *hPtRes, TH1D *hPtRes_mu[NMU], TH1D *hPtResBarrel, TH1D *hPtResBarrel_mu[NMU], TH1D *hPtResEndcap, TH1D *hPtResEndcap_mu[NMU]){
    // Pt Reso Before cuts
    double ptResMu[NMU] = {0}, ptGEN[NMU] = {0}, etaGEN[NMU] = {0}, phiGEN[NMU] = {0}, ptSimGEN[NMU] = {0}, etaSimGEN[NMU] = {0}, phiSimGEN[NMU] = {0};
    int muGen[NMU] = {0};
    for (int i=0; i<GenMatchMu1_SimPt->size(); i++){
        for (int k=0; k<NMU; k++){
            muGen[k] = i;
        }
        Fill_MuonVariablesGen(muGen, ptGEN, etaGEN, phiGEN);
        Fill_MuonVariablesGen_Sim(muGen, ptSimGEN, etaSimGEN, phiSimGEN);
        for (int k=0; k<NMU; k++){
            ptResMu[k] = (ptSimGEN[k] - ptGEN[k])/ptSimGEN[k];
            hPtRes_mu[k]->Fill(ptResMu[k]);
            hPtRes->Fill(ptResMu[k]);
            if(abs(etaGEN[k]) < 1.4){
                hPtResBarrel_mu[k]->Fill(ptResMu[k]);
                hPtResBarrel->Fill(ptResMu[k]);
            }
            else {
                hPtResEndcap_mu[k]->Fill(ptResMu[k]);
                hPtResEndcap->Fill(ptResMu[k]);
            }
        }
    }
}

void ntupleClass_MC::FillHistoResoTriplMass(Int_t mu_Ind[NMU], Int_t mu[NMU], TH1D *hMassTriRes, TH1D *hMassTriResBarrel, TH1D *hMassTriResEndcap){
    // Triplet mass resolution
    double pt[NMU] = {0}, ptNO[NMU] = {0}, eta[NMU] = {0}, phi[NMU] = {0}, pt_res[NMU] = {0}, pt_bis[NMU] = {0}, dm[NMU] = {0};
    TLorentzVector muon[NMU], muon_bis[NMU], mutot, mutot_bis[NMU];
    
    Fill_MuonVariables(mu_Ind, ptNO, eta, phi);
    for(int k=0; k<NMU; k++){
        pt[k] = Muon_BestTrackPt->at(mu[k]);
        pt_bis[k] = pt[k] + Muon_BestTrackPtErr->at(mu[k]);
        muon[k].SetPtEtaPhiM(pt[k], eta[k], phi[k], mumass);
        muon_bis[k].SetPtEtaPhiM(pt_bis[k], eta[k], phi[k], mumass);
    }
    mutot = muon[0] + muon[1] + muon[2];
    mutot_bis[0] = muon_bis[0] + muon[1] + muon[2];
    mutot_bis[1] = muon[0] + muon_bis[1] + muon[2];
    mutot_bis[2] = muon[0] + muon[1] + muon_bis[2];
    for(int k=0; k<NMU; k++){
        dm[k] = mutot_bis[k].M() - mutot.M();
    }
    double trimassReso = sqrt(pow(dm[0],2)+pow(dm[1],2)+pow(dm[2],2))/mutot.M();
    hMassTriRes->Fill(trimassReso);
    if(abs(eta[0]) < 1.4 && abs(eta[1]) < 1.4 && abs(eta[2]) < 1.4)   hMassTriResBarrel->Fill(trimassReso);
    if(abs(eta[0]) >= 1.4 && abs(eta[1]) >= 1.4 && abs(eta[2]) >= 1.4)    hMassTriResEndcap->Fill(trimassReso);
}

void ntupleClass_MC::FillHistoSingleMu(Int_t mu_Ind[NMU], Int_t mu[NMU], TH1D *hist_pt, TH1D *hist_pt_mu[NMU], TH1D *hist_eta, TH1D *hist_eta_mu[NMU], TH1D *hist_phi, TH1D *hVx, TH1D *hVy, TH1D *hVz){
    // Fills histograms w/ variables of single mu
    double pt[NMU] = {0}, eta[NMU] = {0}, phi[NMU] = {0};
    Fill_MuonVariables(mu_Ind, pt, eta, phi);
    for(int i=0; i<NMU; i++){
        hist_pt->Fill(pt[i]);
        hist_pt_mu[i]->Fill(pt[i]);
        hist_eta->Fill(eta[i]);
        hist_eta_mu[i]->Fill(eta[i]);
        hist_phi->Fill(phi[i]);
        hVx->Fill(Muon_vx->at(mu[i]));
        hVy->Fill(Muon_vy->at(mu[i]));
        hVz->Fill(Muon_vz->at(mu[i]));
    }
}

void ntupleClass_MC::FillHistoStepByStep(TString type, Int_t ind, Int_t mu_Ind[NMU], Int_t mu[NMU], Int_t Ncut, TH1D *hPt[NCUTS], TH1D *hPt_mu[NCUTS][NMU], TH1D *hEta[NCUTS], TH1D *hEta_mu[NCUTS][NMU], TH1D *hPhi[NCUTS], TH1D *hVx[NCUTS], TH1D *hVy[NCUTS], TH1D *hVz[NCUTS], TH1D *hPt_tripl[NCUTS], TH1D *hEta_tripl[NCUTS], TH1D *hPhi_tripl[NCUTS], TH1D *hMass_tripl[NCUTS], Int_t IdsummaryDaughter[NCUTS][NPARTICLES], Int_t IdsummaryMother[NCUTS][NPARTICLES], Int_t Idsummary2D[NCUTS][NPARTICLES][NPARTICLES]){
    // Fills the "StepByStep" histograms
    FillHistoSingleMu(mu_Ind, mu, hPt[Ncut], hPt_mu[Ncut], hEta[Ncut], hEta_mu[Ncut], hPhi[Ncut], hVx[Ncut], hVy[Ncut], hVz[Ncut]);
    FillHistoTriplet(ind, hPt_tripl[Ncut], hEta_tripl[Ncut], hPhi_tripl[Ncut], hMass_tripl[Ncut]);
    if (strcmp(type, "data") != 0) Fill_ParticleIdSummary(mu, IdsummaryDaughter[Ncut], IdsummaryMother[Ncut], Idsummary2D[Ncut]);
}

void ntupleClass_MC::FillHistoTriplet(Int_t ind, TH1D *hist_pt, TH1D *hist_eta, TH1D *hist_phi, TH1D *hist_mass){
    // Fills histograms w/ variables of the triplet
    hist_pt->Fill(Triplet_Pt->at(ind));
    hist_eta->Fill(Triplet_Eta->at(ind));
    hist_phi->Fill(Triplet_Phi->at(ind));
    hist_mass->Fill(Triplet_Mass->at(ind));
}

void ntupleClass_MC::InitHistoAC(TH1I *&hNtripl, TH1F *&hChi2Track, TH1D *&hMassTriRes, TH1D *&hMassTriResBarrel, TH1D *&hMassTriResEndcap, TH1D *&hmassdi, TH1F *&hmassQuad, TH1F *&hmassQuad_Zero, TH1D *&hPtRes, TH1D *hPtRes_mu[NMU], TH1D *&hPtResBarrel, TH1D *hPtResBarrel_mu[NMU], TH1D *&hPtResEndcap, TH1D *hPtResEndcap_mu[NMU], TH1D *&hNMatchedStat, TH1D *&hFlightDist, TH1D *&hFlightDist_Signif, TH2D *&hFlightDistvsP, TH1D *&hPtErrOverPt, TH1D *&hPt_tripl_good, TH1D *&hPt_tripl_fake, TH1D *&hDeltaX, TH1D *&hDeltaY, TH1D *&hDeltaZ, TH1D *&hDeltaX_fake, TH1D *&hDeltaY_fake, TH1D *&hDeltaZ_fake){
    // Init histograms for variables After Cuts
    hNtripl = new TH1I("Ntripl", "Ntripl", 5, -0.5, 4.5);
    hNtripl->GetXaxis()->SetTitle("N. triplets survived per event");
    hNtripl->GetYaxis()->SetTitle("N. entries");
    hChi2Track = new TH1F("Chi2Track", "Chi2Track", 27, -0.3, 5.1); //binning 0.2
    hChi2Track->GetXaxis()->SetTitle("#chi^{2} Muon Inner track");
    hChi2Track->GetYaxis()->SetTitle("N. muons");
    hChi2Track->Sumw2();
    hMassTriRes = new TH1D("TriplMassRes", "TriplMassRes", 600, -0.1, 0.1); //binning 0.00033
    hMassTriRes->GetXaxis()->SetTitle("TripletMass Resolution");
    hMassTriRes->GetYaxis()->SetTitle("N. triplets");
    //    hMassTriRes->Sumw2();
    hMassTriResBarrel = new TH1D("TriplMassRes_Barrel", "TriplMassRes_Barrel", 600, -0.1, 0.1); //binning 0.00033
    hMassTriResBarrel->GetXaxis()->SetTitle("TripletMass Resolution Barrel");
    hMassTriResBarrel->GetYaxis()->SetTitle("N. triplets");
    //    hMassTriResBarrel->Sumw2();
    hMassTriResEndcap = new TH1D("TriplMassRes_Endcap", "TriplMassRes_Endcap", 600, -0.1, 0.1); //binning 0.00033
    hMassTriResEndcap->GetXaxis()->SetTitle("TripletMass Resolution Endcap");
    hMassTriResEndcap->GetYaxis()->SetTitle("N. triplets");
    //    hMassTriResEndcap->Sumw2();
    hmassdi = new TH1D("DimuonMass", "DimuonMass", 60, -0.05, 2.95); //binning 50 MeV
    hmassdi->GetXaxis()->SetTitle("Mass(#mu^{+}#mu^{-}) (GeV/c^{2})");
    hmassdi->GetXaxis()->SetTitle("N. entries");
    hmassdi->Sumw2();
    hmassQuad = new TH1F("QuadMuonMass", "QuadMuonMass", 400, -0.05, 79.95); // binning 200 MeV
    hmassQuad->GetXaxis()->SetTitle("Mass(4#mu) (GeV/c^{2})");
    hmassQuad->GetYaxis()->SetTitle("N. entries");
    hmassQuad->Sumw2();
    hmassQuad_Zero = new TH1F("QuadMuonMass_Zero", "QuadMuonMass_Zero", 400, -0.05, 79.95); // binning 200 MeV
    hmassQuad_Zero->GetXaxis()->SetTitle("Mass(4#mu) (GeV/c^{2})");
    hmassQuad_Zero->GetYaxis()->SetTitle("N. entries");
    hmassQuad_Zero->Sumw2();
    hNMatchedStat = new TH1D("NofMatchedStations", "NofMatchedStations", 6, -0.5, 5.5);
    hNMatchedStat->GetXaxis()->SetTitle("N. of matched muon stations");
    hNMatchedStat->GetYaxis()->SetTitle("N. muons");
    hNMatchedStat->Sumw2();
    hFlightDist = new TH1D("FlightDist", "FlightDist", 90, 0., 3.);
    hFlightDist->GetXaxis()->SetTitle("Decay length (cm)");
    hFlightDist->GetYaxis()->SetTitle("N. triplets");
    hFlightDist->Sumw2();
    hFlightDist_Signif = new TH1D("FlightDist_Signif", "FlightDist_Signif", 100, 0., 100.);
    hFlightDist_Signif->GetXaxis()->SetTitle("Decay length significance");
    hFlightDist_Signif->GetYaxis()->SetTitle("N. triplets");
    hFlightDist_Signif->Sumw2();
    hFlightDistvsP = new TH2D("FlightDistvsP", "FlightDistvsP", 20, 0., 0.3, 20, 0., 45.);
    hFlightDistvsP->GetXaxis()->SetTitle("Decay length (cm)");
    hFlightDistvsP->GetYaxis()->SetTitle("p triplet (GeV/c)");
    hPtErrOverPt = new TH1D("PtErrOverPt", "PtErrOverPt", 120, 0, 1.2);
    hPtErrOverPt->GetXaxis()->SetTitle("p_{T} err/ p_{T}");
    hPtErrOverPt->GetYaxis()->SetTitle("N. muons");
    hPtRes = new TH1D("MuonPtRes", "MuonPtRes", 160, -2., 2.); // binning 0.025
    hPtRes->GetXaxis()->SetTitle("Muon p_{T} resolution");
    hPtRes->GetYaxis()->SetTitle("N. muons");
    hPtRes->Sumw2();
    hPtResBarrel = new TH1D("MuonPtResBarrel", "MuonPtResBarrel", 160, -2., 2.); // binning 0.025
    hPtResBarrel->GetXaxis()->SetTitle("Muon p_{T} resolution Barrel");
    hPtResBarrel->GetYaxis()->SetTitle("N. muons");
    hPtResBarrel->Sumw2();
    hPtResEndcap = new TH1D("MuonPtResEndcap", "MuonPtResEndcap", 160, -2., 2.); // binning 0.025
    hPtResEndcap->GetXaxis()->SetTitle("Muon p_{T} resolution Endcap");
    hPtResEndcap->GetYaxis()->SetTitle("N. muons");
    hPtResEndcap->Sumw2();
    for(int k=0; k<NMU; k++){
        TString hptResMuName = "MuonPtRes_mu"; hptResMuName += k+1;
        TString hptResBarrelMuName = "MuonPtResBarrel_mu"; hptResBarrelMuName += k+1;
        TString hptResEndcapMuName = "MuonPtResEndcap_mu"; hptResEndcapMuName += k+1;
        hPtRes_mu[k] = new TH1D(hptResMuName, hptResMuName, 160, -2., 2.); // binning 0.025
        hPtResBarrel_mu[k] = new TH1D(hptResBarrelMuName, hptResBarrelMuName, 160, -2., 2.); // binning 0.025
        hPtResEndcap_mu[k] = new TH1D(hptResEndcapMuName, hptResEndcapMuName, 160, -2., 2.); // binning 0.025
        hptResMuName = "Mu"; hptResMuName += k+1; hptResMuName += " p_{T} resolution";
        hPtRes_mu[k]->GetXaxis()->SetTitle(hptResMuName);
        hPtRes_mu[k]->GetYaxis()->SetTitle("N. muons");
        hptResBarrelMuName = "Mu"; hptResBarrelMuName += k+1; hptResBarrelMuName += " p_{T} resolution Barrel";
        hPtResBarrel_mu[k]->GetXaxis()->SetTitle(hptResBarrelMuName);
        hPtResBarrel_mu[k]->GetYaxis()->SetTitle("N. muons");
        hptResEndcapMuName = "Mu"; hptResEndcapMuName += k+1; hptResEndcapMuName += " p_{T} resolution Endcap";
        hPtResEndcap_mu[k]->GetXaxis()->SetTitle(hptResEndcapMuName);
        hPtResEndcap_mu[k]->GetYaxis()->SetTitle("N. muons");
        hPtRes_mu[k]->Sumw2();
        hPtResBarrel_mu[k]->Sumw2();
        hPtResEndcap_mu[k]->Sumw2();
    }
    // For study on fake triplets
    hPt_tripl_good = new TH1D("Pt_tripl_good", "Pt_tripl_good", 160, -0.05, 39.95); // binning 250 MeV
    hPt_tripl_fake = new TH1D("Pt_tripl_fake", "Pt_tripl_fake", 160, -0.05, 39.95); // binning 250 MeV
    hDeltaX = new TH1D("DeltaX", "DeltaX", 500, 0., 0.5); // binning 0.01 mm
    hDeltaY = new TH1D("DeltaY", "DeltaY", 500, 0., 0.5); // binning 0.01 mm
    hDeltaZ = new TH1D("DeltaZ", "DeltaZ", 1000, 0., 1.); // binning 0.01 mm
    hDeltaX_fake = new TH1D("DeltaX_fake", "DeltaX_fake", 500, 0., 0.5); // binning 0.01 mm
    hDeltaY_fake = new TH1D("DeltaY_fake", "DeltaY_fake", 500, 0., 0.5); // binning 0.01 mm
    hDeltaZ_fake = new TH1D("DeltaZ_fake", "DeltaZ_fake", 1000, 0., 1.); // binning 0.01 mm
    //
}

void ntupleClass_MC::InitHistoBC(TH1D *&hMass_tripl, TH1D *&hChi2Vertex, TH2D *&hMassvsChi2, TH1F *&hMass_quad, TH1F *&hMass_quad_Zero, TH1D *&hMass_di, TH1D *&hMass_di2, TH1D *&hPtRes, TH1D *hPtRes_mu[NMU], TH1D *&hPtResBarrel, TH1D *hPtResBarrel_mu[NMU], TH1D *&hPtResEndcap, TH1D *hPtResEndcap_mu[NMU], TH1I *&hPdgId_Gen, TH1I *&hMotherPdgId_Gen, TH2I *&hPdgId2D_Gen){
    // Init histograms for variables Before Cuts
    hMass_tripl = new TH1D("TripletMass", "TripletMass", 400, -0.05, 19.95); // binning 50 MeV
    hMass_tripl->GetXaxis()->SetTitle("Triplet Mass (GeV/c^{2})");
    hMass_tripl->GetYaxis()->SetTitle("N. triplets");
    hMass_tripl->Sumw2();
    hChi2Vertex = new TH1D("Chi2Vertex", "Chi2Vertex", 40, -1.5, 18.5); // binning 0.5
    hChi2Vertex->GetXaxis()->SetTitle("#chi^{2} vertex");
    hChi2Vertex->GetYaxis()->SetTitle("N. triplets");
    hChi2Vertex->Sumw2();
    hMassvsChi2 = new TH2D("MassvsChi2", "MassvsChi2", 200, 1., 2.5, 200, -5., 20.);
    hMassvsChi2->GetXaxis()->SetTitle("Mass (GeV/c^{2})");
    hMassvsChi2->GetYaxis()->SetTitle("#chi^{2} vertex");
    hMass_quad = new TH1F("QuadMuonMass", "QuadMuonMass", 400, -0.05, 79.95); // binning 200 MeV
    hMass_quad->GetXaxis()->SetTitle("Mass(4#mu) (GeV/c^{2})");
    hMass_quad->GetYaxis()->SetTitle("N. entries");
    hMass_quad->Sumw2();
    hMass_quad_Zero = new TH1F("QuadMuonMass_Zero", "QuadMuonMass_Zero", 400, -0.05, 79.95); // binning 200 MeV
    hMass_quad_Zero->GetXaxis()->SetTitle("Mass(4#mu) (GeV/c^{2})");
    hMass_quad_Zero->GetYaxis()->SetTitle("N. entries");
    hMass_quad_Zero->Sumw2();
    hMass_di = new TH1D("DiMuon_2glbMu_Zero", "DiMuon_2glbMu_Zero", 500, -0.05, 24.95); // binning 50 MeV
    hMass_di->GetXaxis()->SetTitle("Mass(#mu^{+}#mu^{-}) (GeV/c^{2})");
    hMass_di->GetYaxis()->SetTitle("N. entries");
    hMass_di->Sumw2();
    hMass_di2 = new TH1D("DiMuon_Other", "DiMuon_Other", 500, -0.05, 24.95); // binning 50 MeV
    hMass_di2->GetXaxis()->SetTitle("Mass(#mu^{+}#mu^{-}) (GeV/c^{2})");
    hMass_di2->GetYaxis()->SetTitle("N. entries");
    hMass_di2->Sumw2();
    hPtRes = new TH1D("MuonPtRes", "MuonPtRes", 160, -2., 2.); // binning 0.025
    hPtRes->GetXaxis()->SetTitle("Muon p_{T} resolution");
    hPtRes->GetYaxis()->SetTitle("N. muons");
    hPtRes->Sumw2();
    hPtResBarrel = new TH1D("MuonPtResBarrel", "MuonPtResBarrel", 160, -2., 2.); // binning 0.025
    hPtResBarrel->GetXaxis()->SetTitle("Muon p_{T} resolution Barrel");
    hPtResBarrel->GetYaxis()->SetTitle("N. muons");
    hPtResBarrel->Sumw2();
    hPtResEndcap = new TH1D("MuonPtResEndcap", "MuonPtResEndcap", 160, -2., 2.); // binning 0.025
    hPtResEndcap->GetXaxis()->SetTitle("Muon p_{T} resolution Endcap");
    hPtResEndcap->GetYaxis()->SetTitle("N. muons");
    hPtResEndcap->Sumw2();
    for(int k=0; k<NMU; k++){
        TString hptResMuName = "MuonPtRes_mu"; hptResMuName += k+1;
        TString hptResBarrelMuName = "MuonPtResBarrel_mu"; hptResBarrelMuName += k+1;
        TString hptResEndcapMuName = "MuonPtResEndcap_mu"; hptResEndcapMuName += k+1;
        hPtRes_mu[k] = new TH1D(hptResMuName, hptResMuName, 160, -2., 2.); // binning 0.025
        hPtResBarrel_mu[k] = new TH1D(hptResBarrelMuName, hptResBarrelMuName, 160, -2., 2.); // binning 0.025
        hPtResEndcap_mu[k] = new TH1D(hptResEndcapMuName, hptResEndcapMuName, 160, -2., 2.); // binning 0.025
        hptResMuName = "Mu"; hptResMuName += k+1; hptResMuName += " p_{T} resolution";
        hPtRes_mu[k]->GetXaxis()->SetTitle(hptResMuName);
        hPtRes_mu[k]->GetYaxis()->SetTitle("N. muons");
        hptResBarrelMuName = "Mu"; hptResBarrelMuName += k+1; hptResBarrelMuName += " p_{T} resolution Barrel";
        hPtResBarrel_mu[k]->GetXaxis()->SetTitle(hptResBarrelMuName);
        hPtResBarrel_mu[k]->GetYaxis()->SetTitle("N. muons");
        hptResEndcapMuName = "Mu"; hptResEndcapMuName += k+1; hptResEndcapMuName += " p_{T} resolution Endcap";
        hPtResEndcap_mu[k]->GetXaxis()->SetTitle(hptResEndcapMuName);
        hPtResEndcap_mu[k]->GetYaxis()->SetTitle("N. muons");
        hPtRes_mu[k]->Sumw2();
        hPtResBarrel_mu[k]->Sumw2();
        hPtResEndcap_mu[k]->Sumw2();
    }
    hPdgId_Gen = new TH1I("PdgId_Daughter_Gen", "PdgId_Daughter_Gen", NPARTICLES, -0.5, (NPARTICLES-0.5));
    hPdgId_Gen->GetXaxis()->SetTitle("PdgId Daughter");
    hPdgId_Gen->GetYaxis()->SetTitle("N. particles");
    hMotherPdgId_Gen = new TH1I("PdgId_Mother_Gen", "PdgId_Mother_Gen", NPARTICLES, -0.5, (NPARTICLES-0.5));
    hMotherPdgId_Gen->GetXaxis()->SetTitle("PdgId Mother");
    hMotherPdgId_Gen->GetYaxis()->SetTitle("N. particles");
    hPdgId2D_Gen = new TH2I("PdgId2D_Gen", "PdgId2D_Gen", NPARTICLES, -0.5, (NPARTICLES-0.5), NPARTICLES, -0.5, (NPARTICLES-0.5));
    hPdgId2D_Gen->GetXaxis()->SetTitle("PdgId Daughter");
    hPdgId2D_Gen->GetYaxis()->SetTitle("PdgId Mother");
}

void ntupleClass_MC::InitHistoStepByStep_PdgId(TH1I *hPdgId_cut[NCUTS], TH1I *hMotherPdgId_cut[NCUTS], TH2I *hPdgId2D_cut[NCUTS]){
    // Init histograms StepByStep with PdgId
    for (int i=0; i<NCUTS; i++){
        TString histPdgName = "PdgId_Daughter_cut"; histPdgName += i;
        TString histPdgMotherName = "PdgId_Mother_cut"; histPdgMotherName += i;
        TString histPdg2D = "PdgId2D_cut"; histPdg2D += i;
        hPdgId_cut[i] = new TH1I(histPdgName, histPdgName, NPARTICLES, -0.5, (NPARTICLES-0.5));
        hMotherPdgId_cut[i] = new TH1I(histPdgMotherName, histPdgMotherName, NPARTICLES, -0.5, (NPARTICLES-0.5));
        hPdgId2D_cut[i] = new TH2I(histPdg2D, histPdg2D, NPARTICLES, -0.5, (NPARTICLES-0.5), NPARTICLES, -0.5, (NPARTICLES-0.5));
        hPdgId_cut[i]->GetXaxis()->SetTitle("PdgId Daughter");
        hPdgId_cut[i]->GetYaxis()->SetTitle("N. particles");
        hMotherPdgId_cut[i]->GetXaxis()->SetTitle("PdgId Mother");
        hMotherPdgId_cut[i]->GetYaxis()->SetTitle("N. particles");
        hPdgId2D_cut[i]->GetXaxis()->SetTitle("PdgId Daughter");
        hPdgId2D_cut[i]->GetYaxis()->SetTitle("PdgId Mother");
    }
}

void ntupleClass_MC::InitHistoStepByStep_SingleMu(TH1D *hPt[NCUTS], TH1D *hPt_mu[NCUTS][NMU], TH1D *hEta[NCUTS], TH1D *hEta_mu[NCUTS][NMU], TH1D *hPhi[NCUTS], TH1D *hVx[NCUTS], TH1D *hVy[NCUTS], TH1D *hVz[NCUTS]){
    // Init histograms StepByStep w/ variables of single mu
    for (int i=0; i<NCUTS; i++){
        // General mu histo
        TString hptName = "MuonPt_cut"; hptName += i;
        TString hetaName = "MuonEta_cut"; hetaName += i;
        TString hphiName = "MuonPhi_cut"; hphiName += i;
        TString hVxName = "MuonVx_cut"; hVxName += i;
        TString hVyName = "MuonVy_cut"; hVyName += i;
        TString hVzName = "MuonVz_cut"; hVzName += i;
        hPt[i] = new TH1D(hptName, hptName, 100, -0.05, 24.95); // binning 250 MeV
        hEta[i] = new TH1D(hetaName, hetaName, 100, -2.5, 2.5); // binning 0.05
        hPhi[i] = new TH1D(hphiName, hphiName, 140, -3.5, 3.5); // binning 0.05
        hVx[i] = new TH1D(hVxName, hVxName, 300, -15., 15.); // binning 0.1
        hVy[i] = new TH1D(hVyName, hVyName, 300, -15., 15.); // binning 0.1
        hVz[i] = new TH1D(hVzName, hVzName, 1200, -60., 60.); // binning 0.1
        hPt[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        hPt[i]->GetYaxis()->SetTitle("N. muons");
        hEta[i]->GetXaxis()->SetTitle("#eta");
        hEta[i]->GetYaxis()->SetTitle("N. muons");
        hPhi[i]->GetXaxis()->SetTitle("#Phi");
        hPhi[i]->GetYaxis()->SetTitle("N. muons");
        hVx[i]->GetXaxis()->SetTitle("Muon V_{x} (cm)");
        hVx[i]->GetYaxis()->SetTitle("N. muons");
        hVy[i]->GetXaxis()->SetTitle("Muon V_{y} (cm)");
        hVy[i]->GetYaxis()->SetTitle("N. muons");
        hVz[i]->GetXaxis()->SetTitle("Muon V_{z} (cm)");
        hVz[i]->GetYaxis()->SetTitle("N. muons");
        hPt[i]->Sumw2();
        hEta[i]->Sumw2();
        hPhi[i]->Sumw2();
        hVx[i]->Sumw2();
        hVy[i]->Sumw2();
        hVz[i]->Sumw2();
        // Single mu histo
        for(int k=0; k<NMU; k++){
            TString hptMuName = "MuonPt_mu"; hptMuName += k+1; hptMuName += "_cut"; hptMuName += i;
            TString hetaMuName = "MuonEta_mu"; hetaMuName += k+1; hetaMuName += "_cut"; hetaMuName += i;
            hPt_mu[i][k] = new TH1D(hptMuName, hptMuName, 100, -0.05, 24.95); // binning 250 MeV
            hEta_mu[i][k] = new TH1D(hetaMuName, hetaMuName, 100, -2.5, 2.5); // binning 0.05
            hptMuName = "Mu"; hptMuName += k+1; hptMuName += " p_{T} (GeV/c)";
            hPt_mu[i][k]->GetXaxis()->SetTitle(hptMuName);
            hPt_mu[i][k]->GetYaxis()->SetTitle("N. muons");
            hetaMuName = "Mu"; hetaMuName += k+1; hetaMuName += " #eta";
            hEta_mu[i][k]->GetXaxis()->SetTitle(hetaMuName);
            hEta_mu[i][k]->GetYaxis()->SetTitle("N. muons");
            hPt_mu[i][k]->Sumw2();
            hEta_mu[i][k]->Sumw2();
        }
    }
}

void ntupleClass_MC::InitHistoStepByStep_Triplet(TH1D *hPt_tripl[NCUTS], TH1D *hEta_tripl[NCUTS], TH1D *hPhi_tripl[NCUTS], TH1D *hMass_tripl[NCUTS]){
    for(int i=0; i<NCUTS; i++){
        // Init histograms StepByStep w/ the variables of the triplet
        TString hPtTriplName = "Pt triplet_cut"; hPtTriplName += i;
        TString hEtaTriplName = "Eta triplet_cut"; hEtaTriplName += i;
        TString hPhiTriplName = "Phi triplet_cut"; hPhiTriplName += i;
        TString hMassTriplName = "Mass triplet_cut"; hMassTriplName += i;
        hPt_tripl[i] = new TH1D(hPtTriplName, hPtTriplName, 160, -0.05, 39.95); // binning 250 MeV
        hEta_tripl[i] = new TH1D(hEtaTriplName, hEtaTriplName, 100, -2.5, 2.5); // binning 0.05
        hPhi_tripl[i] = new TH1D(hPhiTriplName, hPhiTriplName, 140, -3.5, 3.5); // binning 0.05
        hMass_tripl[i] = new TH1D(hMassTriplName, hMassTriplName, 42, 1.60, 2.02); // binning 10 MeV
        hPt_tripl[i]->GetXaxis()->SetTitle("p_{T} triplet (GeV/c)");
        hPt_tripl[i]->GetYaxis()->SetTitle("N. triplets");
        hEta_tripl[i]->GetXaxis()->SetTitle("#eta triplet");
        hEta_tripl[i]->GetYaxis()->SetTitle("N. triplets");
        hPhi_tripl[i]->GetXaxis()->SetTitle("#Phi triplet");
        hPhi_tripl[i]->GetYaxis()->SetTitle("N. triplets");
        hMass_tripl[i]->GetXaxis()->SetTitle("Mass triplet (GeV/c^{2})");
        hMass_tripl[i]->GetYaxis()->SetTitle("N. triplets");
        hPt_tripl[i]->Sumw2();
        hEta_tripl[i]->Sumw2();
        hPhi_tripl[i]->Sumw2();
        hMass_tripl[i]->Sumw2();
    }
}

Bool_t ntupleClass_MC::isDeltaRGood(Float_t eta1, Float_t eta2, Float_t phi1, Float_t phi2, Float_t DeltaRmax){
    // Given 2 muons the function returns 'true' if DeltaR < DeltaRmax
    float n = TMath::Sqrt(pow((eta1-eta2),2)+pow((phi1-phi2),2));
    if(n < DeltaRmax) return true;
    else return false;
}

Bool_t ntupleClass_MC::isDeltaZGood(Float_t vz1, Float_t vz2, Float_t DeltaZmax){
    // Given 2 muons the function returns 'true' if |DeltaZ| < DeltaZmax
    float n = TMath::Abs(vz2 - vz1);
    if(n<DeltaZmax) return true;
    else return false;
}

Bool_t ntupleClass_MC::isNotAOmega(Double_t dimumass, Double_t sigma){
    // The function return 'true' if the dimuon mass is NOT compatible w/ the mass of an Omega(782) in 2 sigmas
    if (dimumass < (OmegaMass-2*sigma) || dimumass > (OmegaMass+2*sigma))   return true;
    else return false;
}

Bool_t ntupleClass_MC::isNotAPhi(Double_t dimumass, Double_t sigma){
    // The function return 'true' if the dimuon mass is NOT compatible w/ the mass of a Phi(1020) in 2 sigmas
    if (dimumass < (PhiMass-2*sigma) || dimumass > (PhiMass+2*sigma))   return true;
    else return false;
}

Bool_t ntupleClass_MC::isPairDeltaRGood(Int_t ntriplet, Float_t DeltaRmax){
    // The function returns 'true' if all of the 3 possible pairs of muons of the triplet satisfy isDeltaRGood
    if(isDeltaRGood(Mu1_Eta->at(ntriplet), Mu2_Eta->at(ntriplet), Mu1_Phi->at(ntriplet), Mu2_Phi->at(ntriplet), DeltaRmax) == true && isDeltaRGood(Mu2_Eta->at(ntriplet), Mu3_Eta->at(ntriplet), Mu2_Phi->at(ntriplet), Mu3_Phi->at(ntriplet), DeltaRmax) == true && isDeltaRGood(Mu1_Eta->at(ntriplet), Mu3_Eta->at(ntriplet), Mu1_Phi->at(ntriplet), Mu3_Phi->at(ntriplet), DeltaRmax) == true)
        return true;
    else return false;
}

Bool_t ntupleClass_MC::isPairDeltaZGood(Float_t DeltaZ1, Float_t DeltaZ2, Float_t DeltaZ3, Float_t DeltaZmax){
    // The function returns 'true' if all of the 3 possible pairs of muons of the triplet satisfy isDeltaZGood
    if(isDeltaZGood(DeltaZ1, DeltaZ2, DeltaZmax) == true && isDeltaZGood(DeltaZ2, DeltaZ3, DeltaZmax) == true && isDeltaZGood(DeltaZ1, DeltaZ3, DeltaZmax) == true)
        return true;
    else return false;
}

Bool_t ntupleClass_MC::isPairNotAOmega(Double_t dimu[NMU], Double_t sigma){
    // Given 3 muons it checks, for all the pairs of o.s. muons, if their dimuon mass is NOT compatible w/ the Omega(782)mass
    if(dimu[0] != 0 & isNotAOmega(dimu[0], sigma) == false) return false;
    else if(dimu[1] != 0 & isNotAOmega(dimu[1], sigma) == false)    return false;
    else if(dimu[2] != 0 & isNotAOmega(dimu[2], sigma) == false)    return false;
    else    return true;
}

Bool_t ntupleClass_MC::isPairNotAPhi(Double_t dimu[NMU], Double_t sigma){
    // Given 3 muons it checks, for all the pairs of o.s. muons, if their dimuon mass is NOT compatible w/ the Phi mass(1020)
    if(dimu[0] != 0 & isNotAPhi(dimu[0], sigma) == false)   return false;
    else if(dimu[1] != 0 & isNotAPhi(dimu[1], sigma) == false)  return false;
    else if(dimu[2] != 0 & isNotAPhi(dimu[2], sigma) == false)  return false;
    else    return true;
}

void ntupleClass_MC::MatchIndex(TString type, Int_t ind, Int_t mu_Ind[NMU], Int_t mu[NMU]){
    // This function matches the index of muons in different cases (ID or GEN muons)
    mu_Ind[0] = Mu1_TripletIndex->at(ind);
    mu_Ind[1] = Mu2_TripletIndex->at(ind);
    mu_Ind[2] = Mu3_TripletIndex->at(ind);
    if (mu_Ind[0] != ind || mu_Ind[1] != ind || mu_Ind[2] != ind) cout << "Error : Different triplet mu indices!" << endl;
    double pt[NMU] = {0}, eta[NMU] = {0}, phi[NMU] = {0};
    Fill_MuonVariables(mu_Ind, pt, eta, phi);
    for(int k=0; k<NMU; k++){
        if (strcmp(type, "ID") == 0)    mu[k] = MuonFinder(pt[k], eta[k], phi[k]);
        if (strcmp(type, "Gen") == 0)   mu[k] = MuonFinderGen(k+1, pt[k], eta[k], phi[k]);
    }
}

Double_t ntupleClass_MC::MuonFinder(Double_t pt, Double_t eta, Double_t phi){
    // Given the characteristics of a muon (pt, eta, phi), the function returns the index of the corresponding muon in the event
    
    // PRIORITY [in case there is more than 1 muon that matches the conditions]
    
    // *** if (there aren't muons w/ outerTrachChi2 != -999)    ## CASE 0 ##
    //      -> consider the one w/ the smallest innerTrackChi2
    //      -> if innerTrackChi2 are equal:
    //          -> consider the one w/ the biggest number of matches
    //              (if they are equal, print an error message and return the index of the last one)
    //
    // *** if (outerTrack != -999)                              ## CASE 1,2 ##
    //      -> consider the one w/ the smallest outerTrackChi2
    //      -> if outerTrackChi2 are equal:
    //          -> consider the one w/ the smallest innerTrackChi2
    //              -> consider the one w/ the biggest number of matches
    //                   (if they are equal, print an error message and return the index of the last one)
    
    int n=0, m=0, badOuterChi2=0;
    for(int g=0; g<MuonPt->size(); g++){
        if(pt == MuonPt->at(g) && eta == MuonEta->at(g) && phi == MuonPhi->at(g)){
            n++;
            m = g;
            if(Muon_outerTrack_normalizedChi2->at(g) == -999)
                badOuterChi2++;
        }
    }
    
    // MULTIPLE muons (There is more than 1 muon that matches the conditions)
    if(n>1) {
        
        // ####### CASE 0: there aren't muons w/ outerTrachChi2 != -999
        if((n-badOuterChi2) == 0){
            double Chi2InnerTrackmin[2] = {0}; // Chi2InnerTrackmin[1] is the value of Chi2; Chi2InnerTrackmin[0] is the corresponding index
            Chi2InnerTrackmin[1] = 999;
            // find the muon w/ the smallest innerTrackChi2
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_innerTrack_normalizedChi2->at(k) < Chi2InnerTrackmin[1]){
                    Chi2InnerTrackmin[1] = Muon_innerTrack_normalizedChi2->at(k);
                    Chi2InnerTrackmin[0] = k;
                }
            }
            int nMuonChi2InnerMin = 0;
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k)){
                    nMuonChi2InnerMin++;
                }
            }
            if(nMuonChi2InnerMin == 1){
                return Chi2InnerTrackmin[0];
            }
            // find the mu w/ the biggest number of matches
            else {
                double NumberOfMatches[2] = {0}; // like Chi2InnerTrackmin[]
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k)){
                        if(Muon_numberOfMatches->at(k) > NumberOfMatches[1]){
                            NumberOfMatches[1] = Muon_numberOfMatches->at(k);
                            NumberOfMatches[0] = k;
                        }
                    }
                }
                int nMuonMatchesMax = 0;
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k) && NumberOfMatches[1] == Muon_numberOfMatches->at(k)){
                        nMuonMatchesMax++;
                    }
                }
                if(nMuonMatchesMax == 1){
                    return NumberOfMatches[0];
                }
                // They are equal, therefore print an error message and return the index of the last one
                else{
                    cout << "Multiple muons with outerChi2 = -999! " << endl;
                    return NumberOfMatches[0];
                }
            }
        }
        
        // ####### CASE 1:  // There is ONLY 1 muon that has outerTrachChi2 != -999
        else if((n-badOuterChi2) == 1){
            // Find this muon and return its index
            int indexGoodMu = -1;
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) != -999)
                    indexGoodMu = k;
            }
            if (indexGoodMu == -1)  cout << "There is a BUG for sure !!!" << endl;
            return indexGoodMu;
        }
        
        // ####### CASE 2:  // There is more than 1 muon w/ outerTrachChi2 != -999
        else if((n-badOuterChi2) > 1){
            double Chi2OuterTrackmin[2] = {0}; Chi2OuterTrackmin[1] = 999;
            // find the one w/ the smallest outerTrackChi2
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) < Chi2OuterTrackmin[1]){
                    Chi2OuterTrackmin[1] = Muon_outerTrack_normalizedChi2->at(k);
                    Chi2OuterTrackmin[0] = k;
                }
            }
            int nMuonChi2OuterMin = 0;
            for (int k=0; k<MuonPt->size(); k++){
                if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Chi2OuterTrackmin[1] == Muon_outerTrack_normalizedChi2->at(k))
                    nMuonChi2OuterMin++;
            }
            if(nMuonChi2OuterMin == 1)  return Chi2OuterTrackmin[0];
            // if outerTrackChi2 are equal
            else{
                // consider the one w/ the smallest innerTrackChi2
                double Chi2InnerTrackmin[2] = {0};
                Chi2InnerTrackmin[1] = 999;
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) == Chi2OuterTrackmin[1] && Muon_innerTrack_normalizedChi2->at(k) < Chi2InnerTrackmin[1]){
                        Chi2InnerTrackmin[1] = Muon_innerTrack_normalizedChi2->at(k);
                        Chi2InnerTrackmin[0] = k;
                    }
                }
                int nMuonChi2InnerMin = 0;
                for (int k=0; k<MuonPt->size(); k++){
                    if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) == Chi2OuterTrackmin[1] && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k))
                        nMuonChi2InnerMin++;
                }
                if(nMuonChi2InnerMin == 1)  return Chi2InnerTrackmin[0];
                // find the one w/ the biggest number of matches
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
                    int nMuonMatchesMax = 0;
                    for (int k=0; k<MuonPt->size(); k++){
                        if(pt == MuonPt->at(k) && eta == MuonEta->at(k) && phi == MuonPhi->at(k) && Muon_outerTrack_normalizedChi2->at(k) == Chi2OuterTrackmin[1] && Chi2InnerTrackmin[1] == Muon_innerTrack_normalizedChi2->at(k) && NumberOfMatches[1] == Muon_numberOfMatches->at(k))
                            nMuonMatchesMax++;
                    }
                    if(nMuonMatchesMax == 1)    return NumberOfMatches[0];
                    // They are equal, therefore print an error message and return the index of the last one
                    else{
                        cout << "Multiple muons with outerChi2 != -999! " << endl;
                        return NumberOfMatches[0];
                    }
                }
            }
        }
        else {  // Error
            cout << "ERROR in multiple muon number!" << endl;
            return m;
        }
    }   // end multiple muons section
    
    else
        return m;
}

Double_t ntupleClass_MC::MuonFinderGen(Int_t muind, Double_t pt, Double_t eta, Double_t phi){
    // Given the characteristics of a GEN muon, the function return the index of the corresponding GENmuon in the event
    int n=0, m=-999;
    if(muind == 1){ // muon GEN1
        for(int g=0; g<GenMatchMu1_Pt->size(); g++){
            if(pt == GenMatchMu1_Pt->at(g) && eta == GenMatchMu1_Eta->at(g) && phi == GenMatchMu1_Phi->at(g)){
                n++;
                m = g;
            }
        }
    }
    if(muind == 2){  // muon GEN2
        for(int g=0; g<GenMatchMu2_Pt->size(); g++){
            if(pt == GenMatchMu2_Pt->at(g) && eta == GenMatchMu2_Eta->at(g) && phi == GenMatchMu2_Phi->at(g)){
                n++;
                m = g;
            }
        }
    }
    if(muind == 3){  // muon GEN3
        for(int g=0; g<GenMatchMu3_Pt->size(); g++){
            if(pt == GenMatchMu3_Pt->at(g) && eta == GenMatchMu3_Eta->at(g) && phi == GenMatchMu3_Phi->at(g)){
                n++;
                m = g;
            }
        }
    }
    if(n>1) cout << "Error: There is more than one muonGEN " << muind << "  that matches the conditions!" << endl;
    if(n==0) cout << "Error: There are NO muonGEN " << muind << " that match the conditions!" << endl;
    return m;
}

Double_t ntupleClass_MC::MuonP(Double_t pt, Double_t eta, Double_t phi){
    // Given the energy, eta, phi of a muon, the function returns the momentum of the muon
    TLorentzVector mu;
    mu.SetPtEtaPhiM(pt, eta, phi, mumass);
    return mu.P();
}

Float_t ntupleClass_MC::QuadMuonMass(Float_t pt1, Float_t pt2, Float_t pt3, Float_t pt4, Float_t eta1, Float_t eta2, Float_t eta3, Float_t eta4, Float_t phi1, Float_t phi2, Float_t phi3, Float_t phi4){
    // Given the characteristics of 4 muons it returns their invariant mass
    TLorentzVector mu1, mu2, mu3, mu4, mutot;
    mu1.SetPtEtaPhiM(pt1, eta1, phi1, mumass);
    mu2.SetPtEtaPhiM(pt2, eta2, phi2, mumass);
    mu3.SetPtEtaPhiM(pt3, eta3, phi3, mumass);
    mu4.SetPtEtaPhiM(pt4, eta4, phi4, mumass);
    mutot = mu1 + mu2 + mu3 + mu4;
    return mutot.M();
}

void ntupleClass_MC::StudyOnTriplet(TString type, Int_t ind, Int_t mu[NMU], TH1D *hDeltaX, TH1D *hDeltaY, TH1D *hDeltaZ, TH1D *hPt_tripl){
    // Deep study on triplets
    /*
     if(strcmp(type, "good") == 0)
     cout << "GOOD TRIPLET !!!" << endl << endl;
     if(strcmp(type, "bad") == 0)
     cout << "BAD TRIPLET !!!" << endl << endl;
     
     cout << "Triplet  " << endl;
     cout << "X : " << TripletVtx_x->at(ind) << endl;
     cout << "Y : " << TripletVtx_y->at(ind) << endl;
     cout << "Z : " << TripletVtx_z->at(ind) << endl << endl;
     cout << "Mu1 Triplet  " << endl;
     cout << "Pt : " << MuonPt->at(mu1) << endl;
     cout << "X : " << Muon_vx->at(mu1) << endl;
     cout << "Y : " << Muon_vy->at(mu1) << endl;
     cout << "Z : " << Muon_vz->at(mu1) << endl << endl;
     cout << "Mu2 Triplet  " << endl;
     cout << "Pt : " << MuonPt->at(mu2) << endl;
     cout << "X : " << Muon_vx->at(mu2) << endl;
     cout << "Y : " << Muon_vy->at(mu2) << endl;
     cout << "Z : " << Muon_vz->at(mu2) << endl << endl;
     cout << "Mu3 Triplet  " << endl;
     cout << "Pt : " << MuonPt->at(mu3) << endl;
     cout << "X : " << Muon_vx->at(mu3) << endl;
     cout << "Y : " << Muon_vy->at(mu3) << endl;
     cout << "Z : " << Muon_vz->at(mu3) << endl << endl;
     if(strcmp(type, "bad") == 0){
     cout << "Coordinates Other muons" << endl;
     for(int k=0; k<MuonPt->size(); k++){
     if(k != mu1 && k!= mu2 && k!= mu3){
     cout << "muon n. " << k+1 << endl;
     cout << "Pt : " << MuonPt->at(k) << endl;
     cout << "X : " << Muon_vx->at(k) << endl;
     cout << "Y : " << Muon_vy->at(k) << endl;
     cout << "Z : " << Muon_vz->at(k) << endl << endl;
     }
     }
     }
     */
    hDeltaX->Fill(abs(Muon_vx->at(mu[0]) - Muon_vx->at(mu[1])));
    hDeltaX->Fill(abs(Muon_vx->at(mu[0]) - Muon_vx->at(mu[2])));
    hDeltaX->Fill(abs(Muon_vx->at(mu[1]) - Muon_vx->at(mu[2])));
    hDeltaY->Fill(abs(Muon_vy->at(mu[0]) - Muon_vy->at(mu[1])));
    hDeltaY->Fill(abs(Muon_vy->at(mu[0]) - Muon_vy->at(mu[2])));
    hDeltaY->Fill(abs(Muon_vy->at(mu[1]) - Muon_vy->at(mu[2])));
    hDeltaZ->Fill(abs(Muon_vz->at(mu[0]) - Muon_vz->at(mu[1])));
    hDeltaZ->Fill(abs(Muon_vz->at(mu[0]) - Muon_vz->at(mu[2])));
    hDeltaZ->Fill(abs(Muon_vz->at(mu[1]) - Muon_vz->at(mu[2])));
    hPt_tripl->Fill(Triplet_Pt->at(ind));
}

void ntupleClass_MC::TreeFin_Fill(TTree *tree, Int_t ind, Int_t mu_Ind[NMU], Int_t mu[NMU], Double_t &Pmu3, Double_t &cLP, Float_t &tKink, Double_t &segmComp, Double_t &fv_nC, Double_t &fv_dphi3D, Double_t &fv_d3Dsig, Double_t &d0sig, Double_t &mindca_iso, Double_t &trkRel){
    // Fills the tree branches
    TLorentzVector trMu;
    trMu.SetPtEtaPhiM(Mu3_Pt->at(mu_Ind[2]), Mu3_Eta->at(mu_Ind[2]), Mu3_Phi->at(mu_Ind[2]), mumass);
    Pmu3 = trMu.P();
    cLP = 0; tKink = 0; segmComp = 1; double temp[NMU] = {0};
    temp[0] = abs(dxy_mu1->at(mu_Ind[0])/ dxyErr_mu1->at(mu_Ind[0]));
    temp[1] = abs(dxy_mu2->at(mu_Ind[1])/ dxyErr_mu2->at(mu_Ind[1]));
    temp[2] = abs(dxy_mu3->at(mu_Ind[2])/ dxyErr_mu3->at(mu_Ind[2]));
    d0sig = temp[0];
    for (int k=0; k<NMU; k++){
        //  * cLP MAX
        //  * kink MAX
        //  * segmComp MIN
        //  * d0sig MIN
        if (Muon_combinedQuality_chi2LocalPosition->at(mu[k]) > cLP) cLP = Muon_combinedQuality_chi2LocalPosition->at(mu[k]);
        if (MuonTrkKink->at(mu[k]) > tKink) tKink = MuonTrkKink->at(mu[k]);
        if (Muon_segmentCompatibility->at(mu[k]) < segmComp) segmComp = Muon_segmentCompatibility->at(mu[k]);
        if (temp[k] < d0sig) d0sig = temp[k];
    }
    fv_nC = TripletVtx_Chi2->at(ind)/3;
    fv_d3Dsig = FlightDistPVSV_Significance->at(ind);
    tree->Fill();
}

void ntupleClass_MC::TreeFin_Init(TTree *&tree, UInt_t &evt, Double_t &Pmu3, Double_t &cLP, Float_t &tKink, Double_t &segmComp, Double_t &fv_nC, Double_t &fv_dphi3D, Double_t &fv_d3Dsig, Double_t &d0sig, Double_t &mindca_iso, Double_t &trkRel){
    // Set tree branches
    tree->Branch("Event", &evt);
    tree->Branch("P", &Pmu3);
    tree->Branch("cLP", &cLP);
    tree->Branch("tKink", &tKink);
    tree->Branch("segmComp", &segmComp);
    tree->Branch("fv_nC", &fv_nC);
    tree->Branch("fv_dphi3D", &fv_dphi3D);
    tree->Branch("fv_d3Dsig", &fv_d3Dsig);
    tree->Branch("d0sig", &d0sig);
    tree->Branch("mindca_iso", &mindca_iso);
    tree->Branch("trkRel", &trkRel);
}
