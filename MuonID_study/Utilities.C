
#define ntupleClass_muonid_cxx
#define NCUTS 19
#define NPARTICLES 560
#define NMU 3
#define mumass 0.1056583715 // Muon mass in GeV
#define PhiMass 1.019461 // Phi mass in GeV
#define sigmaPhiMass 0.011 //sigma of Phi mass in GeV 
#define OmegaMass 0.78265 // Omega mass in GeV
#define sigmaOmegaMass 0.0085 //sigma of Omega mass in GeV
#define ptmin 2.0

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TMVA/Reader.h>

std::vector<double> pileup_weight;
double pileupFactor = 1;


void ntupleClass_muonid::MatchIndex(TString type, Int_t ind, Int_t mu_Ind[NMU], Int_t mu[NMU]){
    // This function matches the index of muons in different cases (ID or GEN muons)
    mu_Ind[0] = ind; 
    mu_Ind[1] = ind; 
    mu_Ind[2] = ind; 
    if (mu_Ind[0] != ind || mu_Ind[1] != ind || mu_Ind[2] != ind) cout << "Error : Different triplet mu indices!" << endl;
    double pt[NMU] = {0}, eta[NMU] = {0}, phi[NMU] = {0};
    Get_MuonVariables(mu_Ind, pt, eta, phi);
    for(int k=0; k<NMU; k++){
        if (strcmp(type, "ID") == 0)    mu[k] = MuonFinder(pt[k], eta[k], phi[k]);
    }
}


Int_t ntupleClass_muonid::BestTripletFinder(std::vector< Int_t > triplIndex){
    // Given the index of all the triplets of an event that passed all the cuts, it returns the index of the one with the smallest Chi2 of the vertex
    int index = 0; double bestChi2 = 100000;
    int dim = triplIndex.size();
    for(int i=0; i<dim; i++){
        if(TripletVtx_Chi2->at(triplIndex[i])<0) continue;
        if(TripletVtx_Chi2->at(triplIndex[i]) < bestChi2){
            bestChi2 = TripletVtx_Chi2->at(triplIndex[i]);
            index = triplIndex[i];
        }
    }
    return index;
}

Double_t ntupleClass_muonid::MuonFinder(Double_t pt, Double_t eta, Double_t phi){
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

Double_t ntupleClass_muonid::DimuonMass(Int_t mu_index1, Int_t mu_index2){
    // Given the characteristics of 2 muons, if their charge is opposite the function returns their invariant mass, otherwise it returns 0
    double inv = 0;
    double charge1 = MuonCharge->at(mu_index1);
    double charge2 = MuonCharge->at(mu_index2);
    double pt1 = MuonPt->at(mu_index1);
    double pt2 = MuonPt->at(mu_index2);
    double eta1 = MuonEta->at(mu_index1);
    double eta2 = MuonEta->at(mu_index2);
    double phi1 = MuonPhi->at(mu_index1);
    double phi2 = MuonPhi->at(mu_index2);
    double en1 = MuonEnergy->at(mu_index1);
    double en2 = MuonEnergy->at(mu_index2);
    if(charge1 + charge2 != 0)  return inv;
    else {
        TLorentzVector mu1, mu2, mutot;
        mu1.SetPtEtaPhiE(pt1, eta1, phi1, en1);
        mu2.SetPtEtaPhiE(pt2, eta2, phi2, en2);
        mutot = mu1 + mu2;
        return mutot.M();
    }
}

void ntupleClass_muonid::Draw_PdgIdCanvas(TCanvas *canv, TH1I *hist, Int_t Idsummary[NPARTICLES], TString pIdList[NPARTICLES]){
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

void ntupleClass_muonid::Draw_PdgIdCanvas_2D(TCanvas *canv, TH2I *hist, Int_t Idsummary2D[NPARTICLES][NPARTICLES], TString pIdList[NPARTICLES]){
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


void ntupleClass_muonid::Draw_PdgIdCanvasGen(TCanvas *PdgIdCanvas_Gen, TH1I *hPdgId_Gen, Int_t IdsummaryDaughter_Gen[NPARTICLES], TCanvas *PdgIdMotherCanvas_Gen, TH1I *hMotherPdgId_Gen, Int_t IdsummaryMother_Gen[NPARTICLES], TCanvas *PdgIdCanvas2D_Gen, TH2I *hPdgId2D_Gen, Int_t Idsummary2D_Gen[NPARTICLES][NPARTICLES], TString pId[NPARTICLES]){
    // This function draws the Pdg histograms for Gen particles
    PdgIdCanvas_Gen = new TCanvas("PdgId_Daughter_Gen", "PdgId_Daughter_Gen", 0, 0, 1600, 1000);
    Draw_PdgIdCanvas(PdgIdCanvas_Gen, hPdgId_Gen, IdsummaryDaughter_Gen, pId);
    PdgIdMotherCanvas_Gen = new TCanvas("PdgId_Mother_Gen", "PdgId_Mother_Gen", 0, 0, 1600, 1000);
    Draw_PdgIdCanvas(PdgIdMotherCanvas_Gen, hMotherPdgId_Gen, IdsummaryMother_Gen, pId);
    PdgIdCanvas2D_Gen = new TCanvas("PdgId2D_Gen", "PdgId2D_Gen", 0, 0, 1600, 1000);
    Draw_PdgIdCanvas_2D(PdgIdCanvas2D_Gen, hPdgId2D_Gen, Idsummary2D_Gen, pId);
}


void ntupleClass_muonid::FillHistoSingleMu(Int_t mu_Ind[NMU], Int_t mu[NMU], TH1D *hist_pt, TH1D *hist_pt_mu[NMU], TH1D *hist_eta, TH1D *hist_eta_mu[NMU], TH1D *hist_phi, TH1D *hVx, TH1D *hVy, TH1D *hVz){
    // Fills histograms w/ variables of single mu
    double pt[NMU] = {0}, eta[NMU] = {0}, phi[NMU] = {0};
//    Get_MuonVariables(mu_Ind, pt, eta, phi);
    for(int i=0; i<NMU; i++){
        hist_pt->Fill(pt[i], pileupFactor);
        hist_pt_mu[i]->Fill(pt[i], pileupFactor);
        hist_eta->Fill(eta[i], pileupFactor);
        hist_eta_mu[i]->Fill(eta[i], pileupFactor);
        hist_phi->Fill(phi[i], pileupFactor);
        hVx->Fill(Muon_vx->at(mu[i]), pileupFactor);
        hVy->Fill(Muon_vy->at(mu[i]), pileupFactor);
        hVz->Fill(Muon_vz->at(mu[i]), pileupFactor);
    }
}


void ntupleClass_muonid::InitHistoStepByStep_SingleMu(TH1D *hPt[NCUTS], TH1D *hPt_mu[NCUTS][NMU], TH1D *hEta[NCUTS], TH1D *hEta_mu[NCUTS][NMU], TH1D *hPhi[NCUTS], TH1D *hVx[NCUTS], TH1D *hVy[NCUTS], TH1D *hVz[NCUTS]){
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


Float_t ntupleClass_muonid::dR(Float_t eta1, Float_t eta2, Float_t phi1, Float_t phi2){
    auto dp = std::abs(phi1 - phi2);
    auto deta = std::abs(eta1 - eta2);
    if (dp > Float_t(M_PI))
        dp -= Float_t(2 * M_PI);
    Float_t n = TMath::Sqrt(dp*dp + deta*deta);
    return n;
}

Int_t ntupleClass_muonid::genParticleMatch(Int_t muonIndex, bool isVerbose){
    // For each reco muon, it returns the index of matching gen particle
    Int_t genIndex = 999;

    if( !GenParticle_Pt->size()>0 || !GenParticle_Eta->size()>0 || !GenParticle_Phi->size()>0 ) return genIndex;
    else{
        std::vector< std::array<double, 2> > dR_dP_list;
        std::vector< double> dR_dP_prod_list;
        //compute dR vnd dP alue for each genParticle in the event
        for( std::size_t k=0; k<GenParticle_Pt->size(); k++ ){
            std::array<double, 2> tmp;
            Float_t dR_temp = dR( MuonEta->at(muonIndex), GenParticle_Eta->at(k), MuonPhi->at(muonIndex), GenParticle_Phi->at(k));
            Float_t dP_temp = std::abs(GenParticle_Pt->at(k) - MuonPt->at(muonIndex)) / GenParticle_Pt->at(k);
            tmp[0] = dR_temp;
            tmp[1] = dP_temp;
            dR_dP_list.push_back(tmp);
            dR_dP_prod_list.push_back(dR_temp*dP_temp);
        }

        //print dR and dP values
        if(isVerbose){
            cout<<"\n dR | dP | "<<endl;
            for(std::size_t i=0; i<GenParticle_Pt->size(); i++){
                std::array<double, 2> tmp = dR_dP_list.at(i);
                cout<<" "<< tmp[0] <<" | "<<tmp[1] <<" | "<<endl;
            } cout<<"\n";
        }

        //gets the row and column location of the min element
        std::size_t min = (min_element(dR_dP_prod_list.begin(),dR_dP_prod_list.end()) - dR_dP_prod_list.begin());
        // gets the value of the min element
        Float_t val = *min_element(dR_dP_prod_list.begin(),dR_dP_prod_list.end());
        if(isVerbose){
           cout<<" Min element is located at: "<<min<<" and the value is "<<val<<endl;
           cout<<" Matching genParticle "<<min<<": Pt "<<GenParticle_Pt->at(min)<<": Eta "<<GenParticle_Eta->at(min)<<": Phi "<<GenParticle_Phi->at(min)<<endl;
           cout<<" recoMuon "<<muonIndex<<": Pt "<<MuonPt->at(muonIndex)<<": Eta "<<MuonEta->at(muonIndex)<<": Phi "<<MuonPhi->at(muonIndex)<<endl;
        } 
        //store index
        genIndex = min;
        return genIndex;
    }
}


void ntupleClass_muonid::Get_MuonVariables(Int_t mu_Ind[NMU], Double_t pt[NMU], Double_t eta[NMU], Double_t phi[NMU]){
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

Double_t ntupleClass_muonid::MuonP(Double_t pt, Double_t eta, Double_t phi){
    // Given the energy, eta, phi of a muon, the function returns the momentum of the muon
    TVector3 muon;
    muon.SetPtEtaPhi(pt, eta, phi);
    return muon.Mag();
}


void ntupleClass_muonid::TreeFin_Fill(TTree *tree, Int_t genIndex, Int_t muIndex, Double_t &run, Double_t &lumi, Double_t &evt, Double_t &ptetaWeight, Double_t &genP_PdgId, Double_t &genP_Pt, Double_t &genP_Eta, Double_t &genP_Phi, Double_t &genP_MotherPdgId, Double_t &mu_pt, Double_t &mu_eta, Double_t &mu_phi, Double_t &mu_energy, Double_t &mu_charge, Double_t &mu_simPdgId, Double_t &mu_simMotherPdgId, Double_t &mu_simFlavour, Double_t &mu_simType, Double_t &mu_simBX, Double_t &mu_isGlobal, Double_t &mu_isSoft, Double_t &mu_isLoose, Double_t &mu_isTight, Double_t &mu_isPF, Double_t &mu_isRPC, Double_t &mu_isStandAlone, Double_t &mu_isTracker, Double_t &mu_isCalo, Double_t &mu_isQualityValid, Double_t &mu_SoftMVA, Double_t &mu_isTimeValid, Double_t &mu_isIsolationValid, Double_t &mu_numberOfMatchedStations, Double_t &mu_numberOfMatches, Double_t &mu_timeAtIpInOut, Double_t &mu_timeAtIpInOutErr, Double_t &mu_GLnormChi2, Double_t &mu_GLhitPattern_numberOfValidMuonHits, Double_t &mu_trackerLayersWithMeasurement, Double_t &mu_Numberofvalidpixelhits, Double_t &mu_Numberofvalidtrackerhits, Double_t &mu_outerTrack_p, Double_t &mu_outerTrack_eta, Double_t &mu_outerTrack_phi, Double_t &mu_outerTrack_normalizedChi2, Double_t &mu_outerTrack_muonStationsWithValidHits, Double_t &mu_innerTrack_p, Double_t &mu_innerTrack_eta, Double_t &mu_innerTrack_phi, Double_t &mu_innerTrack_validFraction, Double_t &mu_innerTrack_highPurity, Double_t &mu_innerTrack_normalizedChi2, Double_t &mu_QInnerOuter, Double_t &mu_combinedQuality_updatedSta, Double_t &mu_combinedQuality_trkKink, Double_t &mu_combinedQuality_glbKink, Double_t &mu_combinedQuality_trkRelChi2, Double_t &mu_combinedQuality_staRelChi2, Double_t &mu_combinedQuality_chi2LocalPosition, Double_t &mu_combinedQuality_chi2LocalMomentum, Double_t &mu_combinedQuality_localDistance, Double_t &mu_combinedQuality_globalDeltaEtaPhi, Double_t &mu_combinedQuality_tightMatch, Double_t &mu_combinedQuality_glbTrackProbability, Double_t &mu_IP3D_BS, Double_t &mu_IP2D_BS, Double_t &mu_IP3D_PV, Double_t &mu_IP2D_PV, Double_t &mu_validMuonHitComb, Double_t &mu_calEnergy_em, Double_t &mu_calEnergy_emS9, Double_t &mu_calEnergy_emS25, Double_t &mu_calEnergy_had, Double_t &mu_calEnergy_hadS9, Double_t &mu_segmentCompatibility, Double_t &mu_caloCompatibility, Double_t &mu_ptErrOverPt, Double_t &mu_BestTrackPt, Double_t &mu_BestTrackPtErr, Double_t &mu_BestTrackEta, Double_t &mu_BestTrackEtaErr, Double_t &mu_BestTrackPhi, Double_t &mu_BestTrackPhiErr, Double_t &mu_emEt03, Double_t &mu_hadEt03, Double_t &mu_nJets03, Double_t &mu_nTracks03, Double_t &mu_sumPt03, Double_t &mu_hadVetoEt03, Double_t &mu_emVetoEt03, Double_t &mu_trackerVetoPt03, Double_t &mu_emEt05, Double_t &mu_hadEt05, Double_t &mu_nJets05, Double_t &mu_nTracks05, Double_t &mu_sumPt05, Double_t &mu_hadVetoEt05, Double_t &mu_emVetoEt05, Double_t &mu_trackerVetoPt05){
  
    genP_PdgId       = GenParticle_PdgId      ->at(genIndex);
    genP_Pt          = GenParticle_Pt         ->at(genIndex);
    genP_Eta         = GenParticle_Eta        ->at(genIndex);
    genP_Phi         = GenParticle_Phi        ->at(genIndex);
    genP_MotherPdgId = GenParticle_MotherPdgId->at(genIndex);

    mu_pt                                  = MuonPt->at(muIndex);                                              
    mu_eta                                 = MuonEta->at(muIndex);
    mu_phi                                 = MuonPhi->at(muIndex);
    mu_energy                              = MuonEnergy->at(muIndex);
    mu_charge                              = MuonCharge->at(muIndex);
    mu_simPdgId                            = Muon_simPdgId->at(muIndex);
    mu_simMotherPdgId                      = Muon_simMotherPdgId->at(muIndex);
    mu_simFlavour                          = Muon_simFlavour->at(muIndex);
    mu_simType                             = Muon_simType->at(muIndex);
    mu_simBX                               = Muon_simBX->at(muIndex);
                                                                                                  
    mu_isGlobal                            = Muon_isGlobal                ->at(muIndex);   
    mu_isSoft                              = Muon_isSoft                  ->at(muIndex);
    mu_isLoose                             = Muon_isLoose                 ->at(muIndex);
    mu_isTight                             = Muon_isTight                 ->at(muIndex);
    mu_isPF                                = Muon_isPF->at(muIndex);
    mu_isRPC                               = Muon_isRPCMuon->at(muIndex);
    mu_isStandAlone                        = Muon_isStandAloneMuon->at(muIndex);
    mu_isTracker                           = Muon_isTrackerMuon->at(muIndex);
    mu_isCalo                              = Muon_isCaloMuon->at(muIndex);
    mu_isQualityValid                      = Muon_isQualityValid->at(muIndex);
    mu_SoftMVA                             = Muon_SoftMVA_Val->at(muIndex);                   //new
    mu_isTimeValid                         = Muon_isTimeValid->at(muIndex);
    mu_isIsolationValid                    = Muon_isIsolationValid->at(muIndex);
    mu_numberOfMatchedStations             = Muon_numberOfMatchedStations->at(muIndex);
    mu_numberOfMatches                     = Muon_numberOfMatches->at(muIndex);
                                                                                                  
    mu_timeAtIpInOut                       = Muon_timeAtIpInOut                       ->at(muIndex);
    mu_timeAtIpInOutErr                    = Muon_timeAtIpInOutErr                    ->at(muIndex);
    mu_GLnormChi2                          = Muon_GLnormChi2                          ->at(muIndex);
    mu_GLhitPattern_numberOfValidMuonHits  = Muon_GLhitPattern_numberOfValidMuonHits    ->at(muIndex);
    mu_trackerLayersWithMeasurement        = Muon_trackerLayersWithMeasurement        ->at(muIndex);
    mu_Numberofvalidpixelhits              = Muon_Numberofvalidpixelhits              ->at(muIndex);
    mu_Numberofvalidtrackerhits            = Muon_Numberofvalidtrackerhits            ->at(muIndex);  //new
                                            
    mu_outerTrack_p                        = Muon_outerTrack_p                        ->at(muIndex); 
    mu_outerTrack_eta                      = Muon_outerTrack_eta                      ->at(muIndex);
    mu_outerTrack_phi                      = Muon_outerTrack_phi                      ->at(muIndex);
    mu_outerTrack_normalizedChi2           = Muon_outerTrack_normalizedChi2           ->at(muIndex);
    mu_outerTrack_muonStationsWithValidHits= Muon_outerTrack_muonStationsWithValidHits->at(muIndex);
                                                                                                
    mu_innerTrack_p                        = Muon_innerTrack_p->at(muIndex);
    mu_innerTrack_eta                      = Muon_innerTrack_eta->at(muIndex);
    mu_innerTrack_phi                      = Muon_innerTrack_phi->at(muIndex);
    mu_innerTrack_validFraction            = Muon_innerTrack_ValidFraction->at(muIndex); //new
    mu_innerTrack_highPurity               = Muon_innerTrack_highPurity->at(muIndex);    //new
    mu_innerTrack_normalizedChi2           = Muon_innerTrack_normalizedChi2->at(muIndex);
    mu_QInnerOuter                         = Muon_QInnerOuter->at(muIndex);
                                           
                                           
    mu_combinedQuality_updatedSta          = Muon_combinedQuality_updatedSta->at(muIndex);
    mu_combinedQuality_trkKink             = Muon_combinedQuality_trkKink->at(muIndex);
    mu_combinedQuality_glbKink             = Muon_combinedQuality_glbKink->at(muIndex);
    mu_combinedQuality_trkRelChi2          = Muon_combinedQuality_trkRelChi2->at(muIndex);
    mu_combinedQuality_staRelChi2          = Muon_combinedQuality_staRelChi2->at(muIndex);
    mu_combinedQuality_chi2LocalPosition   = Muon_combinedQuality_chi2LocalPosition->at(muIndex);
    mu_combinedQuality_chi2LocalMomentum   = Muon_combinedQuality_chi2LocalMomentum->at(muIndex);
    mu_combinedQuality_localDistance       = Muon_combinedQuality_localDistance->at(muIndex);
    mu_combinedQuality_globalDeltaEtaPhi   = Muon_combinedQuality_globalDeltaEtaPhi->at(muIndex);
    mu_combinedQuality_tightMatch          = Muon_combinedQuality_tightMatch->at(muIndex);
    mu_combinedQuality_glbTrackProbability = Muon_combinedQuality_glbTrackProbability->at(muIndex);

  //  mu_IP3D_BS                             = Muon_IP3D_BS->at(muIndex);           //new  
  //  mu_IP2D_BS                             = Muon_IP2D_BS->at(muIndex);           //new
  //  mu_IP3D_PV                             = Muon_IP3D_PV->at(muIndex);           //new
  //  mu_IP2D_PV                             = Muon_IP2D_PV->at(muIndex);           //new

    mu_validMuonHitComb                    = Muon_validMuonHitComb->at(muIndex);  //new
                                          
    mu_calEnergy_em                        = Muon_calEnergy_em->at(muIndex);
    mu_calEnergy_emS9                      = Muon_calEnergy_emS9->at(muIndex);
    mu_calEnergy_emS25                     = Muon_calEnergy_emS25->at(muIndex);
    mu_calEnergy_had                       = Muon_calEnergy_had->at(muIndex);
    mu_calEnergy_hadS9                     = Muon_calEnergy_hadS9->at(muIndex);
                                          
    mu_segmentCompatibility                = Muon_segmentCompatibility->at(muIndex);
    mu_caloCompatibility                   = Muon_caloCompatibility->at(muIndex);
                                           
    mu_ptErrOverPt                         = Muon_ptErrOverPt->at(muIndex);
    mu_BestTrackPt                         = Muon_BestTrackPt->at(muIndex);
    mu_BestTrackPtErr                      = Muon_BestTrackPtErr->at(muIndex);
    mu_BestTrackEta                        = Muon_BestTrackEta->at(muIndex);
    mu_BestTrackEtaErr                     = Muon_BestTrackEtaErr->at(muIndex);
    mu_BestTrackPhi                        = Muon_BestTrackPhi->at(muIndex);
    mu_BestTrackPhiErr                     = Muon_BestTrackPhiErr->at(muIndex);
                                                              
    mu_emEt03                              = Muon_emEt03->at(muIndex);
    mu_hadEt03                             = Muon_hadEt03->at(muIndex);
    mu_nJets03                             = Muon_nJets03->at(muIndex);
    mu_nTracks03                           = Muon_nTracks03->at(muIndex);
    mu_sumPt03                             = Muon_sumPt03->at(muIndex);
    mu_hadVetoEt03                         = Muon_hadVetoEt03->at(muIndex);
    mu_emVetoEt03                          = Muon_emVetoEt03->at(muIndex);
    mu_trackerVetoPt03                     = Muon_trackerVetoPt03->at(muIndex);

    mu_emEt05                              = Muon_emEt05->at(muIndex);
    mu_hadEt05                             = Muon_hadEt05->at(muIndex);
    mu_nJets05                             = Muon_nJets05->at(muIndex);
    mu_nTracks05                           = Muon_nTracks05->at(muIndex);
    mu_sumPt05                             = Muon_sumPt05->at(muIndex);
    mu_hadVetoEt05                         = Muon_hadVetoEt05->at(muIndex);
    mu_emVetoEt05                          = Muon_emVetoEt05->at(muIndex);
    mu_trackerVetoPt05                     = Muon_trackerVetoPt05->at(muIndex);

    tree->Fill();
}

void ntupleClass_muonid::TreeFin_Init(TTree *&tree_, Double_t &run, Double_t &lumi, Double_t &evt, Double_t &ptetaWeight, Double_t &genP_PdgId, Double_t &genP_Pt, Double_t &genP_Eta, Double_t &genP_Phi, Double_t &genP_MotherPdgId, Double_t &mu_pt, Double_t &mu_eta, Double_t &mu_phi, Double_t &mu_energy, Double_t &mu_charge, Double_t &mu_simPdgId, Double_t &mu_simMotherPdgId, Double_t &mu_simFlavour, Double_t &mu_simType, Double_t &mu_simBX, Double_t &mu_isGlobal, Double_t &mu_isSoft, Double_t &mu_isLoose, Double_t &mu_isTight, Double_t &mu_isPF, Double_t &mu_isRPC, Double_t &mu_isStandAlone, Double_t &mu_isTracker, Double_t &mu_isCalo, Double_t &mu_isQualityValid, Double_t &mu_SoftMVA, Double_t &mu_isTimeValid, Double_t &mu_isIsolationValid, Double_t &mu_numberOfMatchedStations, Double_t &mu_numberOfMatches, Double_t &mu_timeAtIpInOut, Double_t &mu_timeAtIpInOutErr, Double_t &mu_GLnormChi2, Double_t &mu_GLhitPattern_numberOfValidMuonHits, Double_t &mu_trackerLayersWithMeasurement, Double_t &mu_Numberofvalidpixelhits, Double_t &mu_Numberofvalidtrackerhits, Double_t &mu_outerTrack_p, Double_t &mu_outerTrack_eta, Double_t &mu_outerTrack_phi, Double_t &mu_outerTrack_normalizedChi2, Double_t &mu_outerTrack_muonStationsWithValidHits, Double_t &mu_innerTrack_p, Double_t &mu_innerTrack_eta, Double_t &mu_innerTrack_phi, Double_t &mu_innerTrack_validFraction, Double_t &mu_innerTrack_highPurity, Double_t &mu_innerTrack_normalizedChi2, Double_t &mu_QInnerOuter, Double_t &mu_combinedQuality_updatedSta, Double_t &mu_combinedQuality_trkKink, Double_t &mu_combinedQuality_glbKink, Double_t &mu_combinedQuality_trkRelChi2, Double_t &mu_combinedQuality_staRelChi2, Double_t &mu_combinedQuality_chi2LocalPosition, Double_t &mu_combinedQuality_chi2LocalMomentum, Double_t &mu_combinedQuality_localDistance, Double_t &mu_combinedQuality_globalDeltaEtaPhi, Double_t &mu_combinedQuality_tightMatch, Double_t &mu_combinedQuality_glbTrackProbability, Double_t &mu_IP3D_BS, Double_t &mu_IP2D_BS, Double_t &mu_IP3D_PV, Double_t &mu_IP2D_PV, Double_t &mu_validMuonHitComb,  Double_t &mu_calEnergy_em, Double_t &mu_calEnergy_emS9, Double_t &mu_calEnergy_emS25, Double_t &mu_calEnergy_had, Double_t &mu_calEnergy_hadS9, Double_t &mu_segmentCompatibility, Double_t &mu_caloCompatibility, Double_t &mu_ptErrOverPt, Double_t &mu_BestTrackPt, Double_t &mu_BestTrackPtErr, Double_t &mu_BestTrackEta, Double_t &mu_BestTrackEtaErr, Double_t &mu_BestTrackPhi, Double_t &mu_BestTrackPhiErr, Double_t &mu_emEt03, Double_t &mu_hadEt03, Double_t &mu_nJets03, Double_t &mu_nTracks03, Double_t &mu_sumPt03, Double_t &mu_hadVetoEt03, Double_t &mu_emVetoEt03, Double_t &mu_trackerVetoPt03, Double_t &mu_emEt05, Double_t &mu_hadEt05, Double_t &mu_nJets05, Double_t &mu_nTracks05, Double_t &mu_sumPt05, Double_t &mu_hadVetoEt05, Double_t &mu_emVetoEt05, Double_t &mu_trackerVetoPt05){
        // Set tree branches
        tree_->Branch("run", &run);
        tree_->Branch("lumi", &lumi);
        tree_->Branch("evt", &evt);
        tree_->Branch("ptetaWeight", &ptetaWeight);

        tree_->Branch("genP_PdgId", &genP_PdgId);
        tree_->Branch("genP_Pt", &genP_Pt);
        tree_->Branch("genP_Eta", &genP_Eta);
        tree_->Branch("genP_Phi", &genP_Phi);
        tree_->Branch("genP_MotherPdgId", &genP_MotherPdgId);

        tree_->Branch("mu_pt",&mu_pt);
        tree_->Branch("mu_eta",&mu_eta);
        tree_->Branch("mu_phi",&mu_phi);
        tree_->Branch("mu_energy", &mu_energy);
        tree_->Branch("mu_charge", &mu_charge);
        tree_->Branch("mu_simPdgId", &mu_simPdgId);
        tree_->Branch("mu_simMotherPdgId", &mu_simMotherPdgId);
        tree_->Branch("mu_simFlavour", &mu_simFlavour);
	tree_->Branch("mu_simType", &mu_simType);
	tree_->Branch("mu_simBX", &mu_simBX);

        tree_->Branch("mu_isGlobal", &mu_isGlobal);
        tree_->Branch("mu_isSoft", &mu_isSoft);
        tree_->Branch("mu_isLoose", &mu_isLoose);
        tree_->Branch("mu_isTight", &mu_isTight);
        tree_->Branch("mu_isPF", &mu_isPF);
        tree_->Branch("mu_isRPC", &mu_isRPC);
        tree_->Branch("mu_isStandAlone", &mu_isStandAlone);
        tree_->Branch("mu_isTracker", &mu_isTracker);
        tree_->Branch("mu_isCalo", &mu_isCalo);
        tree_->Branch("mu_isQualityValid", &mu_isQualityValid);
        tree_->Branch("mu_SoftMVA", &mu_SoftMVA);
        tree_->Branch("mu_isTimeValid", &mu_isTimeValid);
        tree_->Branch("mu_isIsolationValid", &mu_isIsolationValid);
        tree_->Branch("mu_numberOfMatchedStations", &mu_numberOfMatchedStations);
        tree_->Branch("mu_numberOfMatches", &mu_numberOfMatches);

        tree_->Branch("mu_timeAtIpInOut",&mu_timeAtIpInOut);
        tree_->Branch("mu_timeAtIpInOutErr",&mu_timeAtIpInOutErr);
        tree_->Branch("mu_GLnormChi2", &mu_GLnormChi2);
        tree_->Branch("mu_GLhitPattern_numberOfValidMuonHits", &mu_GLhitPattern_numberOfValidMuonHits);

        tree_->Branch("mu_trackerLayersWithMeasurement", &mu_trackerLayersWithMeasurement);
        tree_->Branch("mu_Numberofvalidpixelhits", &mu_Numberofvalidpixelhits);
        tree_->Branch("mu_Numberofvalidtrackerhits", &mu_Numberofvalidtrackerhits);
        
        tree_->Branch("mu_outerTrack_p", &mu_outerTrack_p);
        tree_->Branch("mu_outerTrack_eta", &mu_outerTrack_eta);
        tree_->Branch("mu_outerTrack_phi", &mu_outerTrack_phi);
        tree_->Branch("mu_outerTrack_normalizedChi2", &mu_outerTrack_normalizedChi2);
        tree_->Branch("mu_outerTrack_muonStationsWithValidHits", &mu_outerTrack_muonStationsWithValidHits);

        tree_->Branch("mu_innerTrack_p", &mu_innerTrack_p);
        tree_->Branch("mu_innerTrack_eta", &mu_innerTrack_eta);
        tree_->Branch("mu_innerTrack_phi", &mu_innerTrack_phi);
        tree_->Branch("mu_innerTrack_validFraction", &mu_innerTrack_validFraction);
        tree_->Branch("mu_innerTrack_highPurity", &mu_innerTrack_highPurity);
        tree_->Branch("mu_innerTrack_normalizedChi2", &mu_innerTrack_normalizedChi2);
        tree_->Branch("mu_QInnerOuter", &mu_QInnerOuter);
        
        
        tree_->Branch("mu_combinedQuality_updatedSta", &mu_combinedQuality_updatedSta);
        tree_->Branch("mu_combinedQuality_trkKink", &mu_combinedQuality_trkKink);
        tree_->Branch("mu_combinedQuality_glbKink", &mu_combinedQuality_glbKink);
        tree_->Branch("mu_combinedQuality_trkRelChi2", &mu_combinedQuality_trkRelChi2);
        tree_->Branch("mu_combinedQuality_staRelChi2", &mu_combinedQuality_staRelChi2);
        tree_->Branch("mu_combinedQuality_chi2LocalPosition", &mu_combinedQuality_chi2LocalPosition);
        tree_->Branch("mu_combinedQuality_chi2LocalMomentum", &mu_combinedQuality_chi2LocalMomentum);
        tree_->Branch("mu_combinedQuality_localDistance", &mu_combinedQuality_localDistance);
        tree_->Branch("mu_combinedQuality_globalDeltaEtaPhi", &mu_combinedQuality_globalDeltaEtaPhi);
        tree_->Branch("mu_combinedQuality_tightMatch", &mu_combinedQuality_tightMatch); 
        tree_->Branch("mu_combinedQuality_glbTrackProbability", &mu_combinedQuality_glbTrackProbability);

        tree_->Branch("mu_IP3D_BS", &mu_IP3D_BS);
        tree_->Branch("mu_IP2D_BS", &mu_IP2D_BS);
        tree_->Branch("mu_IP3D_PV", &mu_IP3D_PV);
        tree_->Branch("mu_IP2D_PV", &mu_IP2D_PV);

        tree_->Branch("mu_validMuonHitComb", &mu_validMuonHitComb);
 
        tree_->Branch("mu_calEnergy_em", &mu_calEnergy_em);
        tree_->Branch("mu_calEnergy_emS9", &mu_calEnergy_emS9);
        tree_->Branch("mu_calEnergy_emS25", &mu_calEnergy_emS25);
        tree_->Branch("mu_calEnergy_had", &mu_calEnergy_had);
        tree_->Branch("mu_calEnergy_hadS9", &mu_calEnergy_hadS9);
        
        tree_->Branch("mu_segmentCompatibility", &mu_segmentCompatibility);
        tree_->Branch("mu_caloCompatibility", &mu_caloCompatibility);
        
        tree_->Branch("mu_ptErrOverPt", &mu_ptErrOverPt);
	tree_->Branch("mu_BestTrackPt", &mu_BestTrackPt);
        tree_->Branch("mu_BestTrackPtErr", &mu_BestTrackPtErr);
        tree_->Branch("mu_BestTrackEta", &mu_BestTrackEta);
	tree_->Branch("mu_BestTrackEtaErr", &mu_BestTrackEtaErr);
	tree_->Branch("mu_BestTrackPhi", &mu_BestTrackPhi);
	tree_->Branch("mu_BestTrackPhiErr", &mu_BestTrackPhiErr);

        tree_->Branch("mu_emEt03", &mu_emEt03);
        tree_->Branch("mu_hadEt03", &mu_hadEt03);
        tree_->Branch("mu_nJets03", &mu_nJets03);
        tree_->Branch("mu_nTracks03", &mu_nTracks03);
        tree_->Branch("mu_sumPt03", &mu_sumPt03);
        tree_->Branch("mu_hadVetoEt03", &mu_hadVetoEt03);
        tree_->Branch("mu_emVetoEt03", &mu_emVetoEt03);
        tree_->Branch("mu_trackerVetoPt03", &mu_trackerVetoPt03);

        tree_->Branch("mu_emEt05", &mu_emEt05);
        tree_->Branch("mu_hadEt05", &mu_hadEt05);
        tree_->Branch("mu_nJets05", &mu_nJets05);
        tree_->Branch("mu_nTracks05", &mu_nTracks05);
        tree_->Branch("mu_sumPt05", &mu_sumPt05);
        tree_->Branch("mu_hadVetoEt05", &mu_hadVetoEt05);
        tree_->Branch("mu_emVetoEt05", &mu_emVetoEt05);
        tree_->Branch("mu_trackerVetoPt05", &mu_trackerVetoPt05);
    
}

