
#ifndef ResultsComparison_new_h
#define ResultsComparison_new_h

#define Nhist_BC 17
#define Nhist_StepByStep 128
#define Nhist_AC 30

void HistName_BC(TString hname[Nhist_BC]){
    hname[0] = "BeforeCuts/TripletMass";
    hname[1] = "BeforeCuts/Chi2Vertex";
    hname[2] = "BeforeCuts/QuadMuonMass";
    hname[3] = "BeforeCuts/QuadMuonMass_Zero";
    hname[4] = "BeforeCuts/DiMuon_2glbMu_Zero";
    hname[5] = "BeforeCuts/DiMuon_Other";
    hname[6] = "BeforeCuts/MuonPtRes";
    hname[7] = "BeforeCuts/MuonPtResBarrel";
    hname[8] = "BeforeCuts/MuonPtResEndcap";
    hname[9] = "BeforeCuts/MuonPtRes_mu1";
    hname[10] = "BeforeCuts/MuonPtResBarrel_mu1";
    hname[11] = "BeforeCuts/MuonPtResEndcap_mu1";
    hname[12] = "BeforeCuts/MuonPtRes_mu2";
    hname[13] = "BeforeCuts/MuonPtResBarrel_mu2";
    hname[14] = "BeforeCuts/MuonPtResEndcap_mu2";
    hname[15] = "BeforeCuts/hNPileUp";
    hname[16] = "BeforeCuts/hNPrimaryVertices";
}

void HistName_StepByStep(TString hname[Nhist_StepByStep]){
    // SingleMu
    TString name[12];
    int i=0;
    for (int k=0; k<8; k++){
        name[0] = "StepByStep/SingleMu/MuonPt_cut";
        name[1] = "StepByStep/SingleMu/MuonEta_cut";
        name[2] = "StepByStep/SingleMu/MuonPhi_cut";
        name[3] = "StepByStep/SingleMu/MuonVx_cut";
        name[4] = "StepByStep/SingleMu/MuonVy_cut";
        name[5] = "StepByStep/SingleMu/MuonVz_cut";
        name[6] = "StepByStep/SingleMu/MuonPt_mu1_cut";
        name[7] = "StepByStep/SingleMu/MuonEta_mu1_cut";
        name[8] = "StepByStep/SingleMu/MuonPt_mu2_cut";
        name[9] = "StepByStep/SingleMu/MuonEta_mu2_cut";
        name[10] = "StepByStep/SingleMu/MuonPt_mu3_cut";
        name[11] = "StepByStep/SingleMu/MuonEta_mu3_cut";
        for (int j=0; j<12; j++){
            name[j] += k;
            hname[i] = name[j];
            i++;
        }
    }
    // Triplet
    TString tname[4];
    i=96;
    for (int k=0; k<8; k++){
        name[0] = "StepByStep/Triplet/Pt triplet_cut";
        name[1] = "StepByStep/Triplet/Eta triplet_cut";
        name[2] = "StepByStep/Triplet/Phi triplet_cut";
        name[3] = "StepByStep/Triplet/Mass triplet_cut";
        for (int j=0; j<4; j++){
            name[j] += k;
            hname[i] = name[j];
            i++;
        }
    }
}

void HistName_AC(TString hname[Nhist_AC]){
    hname[0] = "AfterCuts/Chi2Track";
    hname[1] = "AfterCuts/TriplMassRes";
    hname[2] = "AfterCuts/TriplMassRes_Barrel";
    hname[3] = "AfterCuts/TriplMassRes_Endcap";
    hname[4] = "AfterCuts/DimuonMass";
    hname[5] = "AfterCuts/QuadMuonMass";
    hname[6] = "AfterCuts/QuadMuonMass_Zero";
    hname[7] = "AfterCuts/NofMatchedStations";
    hname[8] = "AfterCuts/FlightDist";
    hname[9] = "AfterCuts/FlightDist_Signif";
    hname[10] = "AfterCuts/PtErrOverPt";
    hname[11] = "AfterCuts/MuonPtRes";
    hname[12] = "AfterCuts/MuonPtResBarrel";
    hname[13] = "AfterCuts/MuonPtResEndcap";
    hname[14] = "AfterCuts/MuonPtRes_mu1";
    hname[15] = "AfterCuts/MuonPtResBarrel_mu1";
    hname[16] = "AfterCuts/MuonPtResEndcap_mu1";
    hname[17] = "AfterCuts/MuonPtRes_mu2";
    hname[18] = "AfterCuts/MuonPtResBarrel_mu2";
    hname[19] = "AfterCuts/MuonPtResEndcap_mu2";
    hname[20] = "AfterCuts/Pt_tripl_good";
    hname[21] = "AfterCuts/Pt_tripl_fake";
    hname[22] = "AfterCuts/DeltaX";
    hname[23] = "AfterCuts/DeltaY";
    hname[24] = "AfterCuts/DeltaZ";
    hname[25] = "AfterCuts/DeltaX_fake";
    hname[26] = "AfterCuts/DeltaY_fake";
    hname[27] = "AfterCuts/DeltaZ_fake";
    hname[28] = "AfterCuts/hNPileUp";
    hname[29] = "AfterCuts/hNPrimaryVertices";
}


void XaxisName_BC(TString hname[Nhist_BC]){
    hname[0] = "Triplet Mass (GeV/c^{2})";
    hname[1] = "#chi^{2} vertex";
    hname[2] = "Mass(4#mu) (GeV/c^{2})";
    hname[3] = "Mass(4#mu) (GeV/c^{2})";
    hname[4] = "Mass(#mu^{+}#mu^{-}) (GeV/c^{2})";
    hname[5] = "Mass(#mu^{+}#mu^{-}) (GeV/c^{2})";
    hname[6] = "Muon p_{T} resolution";
    hname[7] = "Muon p_{T} resolution Barrel";
    hname[8] = "Muon p_{T} resolution Endcap";
    hname[9] = "MuonPtRes_mu1";
    hname[10] = "MuonPtResBarrel_mu1";
    hname[11] = "MuonPtResEndcap_mu1";
    hname[12] = "MuonPtRes_mu2";
    hname[13] = "MuonPtResBarrel_mu2";
    hname[14] = "MuonPtResEndcap_mu2";
    hname[15] = "N. pile-up interactions";
    hname[16] = "N. primary vertices";
}

void XaxisName_StepByStep(TString hname[Nhist_StepByStep]){
    // SingleMu
    TString xname[12];
    xname[0] = "p_{T} (GeV/c)";
    xname[1] = "#eta";
    xname[2] = "#Phi";
    xname[3] = "Muon V_{x} (cm)";
    xname[4] = "Muon V_{y} (cm)";
    xname[5] = "Muon V_{z} (cm)";
    xname[6] = "Mu1 p_{T} (GeV/c)";
    xname[7] = "Mu1 #eta";
    xname[8] = "Mu2 p_{T} (GeV/c)";
    xname[9] = "Mu2 #eta";
    xname[10] = "Mu3 p_{T} (GeV/c)";
    xname[11] = "Mu3 #eta";
    int i=0;
    for (int k=0; k<8; k++){
        for (int j=0; j<12; j++){
            hname[i] = xname[j];
            i++;
        }
    }
    //Triplet
    TString txname[4];
    txname[0] = "p_{T} triplet (GeV/c)";
    txname[1] = "#eta triplet";
    txname[2] = "#Phi triplet";
    txname[3] = "Mass triplet (GeV/c^{2})";
    i=96;
    for (int k=0; k<8; k++){
        for (int j=0; j<4; j++){
            hname[i] = txname[j];
            i++;
        }
    }
}

void XaxisName_AC(TString hname[Nhist_AC]){
    hname[0] = "#chi^{2} Muon Inner track";
    hname[1] = "TripletMass Resolution";
    hname[2] = "TripletMass Resolution Barrel";
    hname[3] = "TripletMass Resolution Endcap";
    hname[4] = "Mass(#mu^{+}#mu^{-}) (GeV/c^{2})";
    hname[5] = "Mass(4#mu) (GeV/c^{2})";
    hname[6] = "Mass(4#mu) (GeV/c^{2})";
    hname[7] = "N. of matched muon stations";
    hname[8] = "Decay length (cm)";
    hname[9] = "Decay length significance";
    hname[10] = "p_{T} err/ p_{T}";
    hname[11] = "Muon p_{T} resolution";
    hname[12] = "Muon p_{T} resolution Barrel";
    hname[13] = "Muon p_{T} resolution Endcap";
    hname[14] = "MuonPtRes_mu1";
    hname[15] = "MuonPtResBarrel_mu1";
    hname[16] = "MuonPtResEndcap_mu1";
    hname[17] = "MuonPtRes_mu2";
    hname[18] = "MuonPtResBarrel_mu2";
    hname[19] = "MuonPtResEndcap_mu2";
    hname[20] = "p_{T} triplet GOOD (GeV/c)";
    hname[21] = "p_{T} triplet FAKE (GeV/c)";
    hname[22] = "#Delta X";
    hname[23] = "#Delta Y";
    hname[24] = "#Delta Z";
    hname[25] = "#Delta X fake";
    hname[26] = "#Delta Y fake";
    hname[27] = "#Delta Z fake";
    hname[28] = "N. pile-up interactions";
    hname[29] = "N. primary vertices";
}


#endif /* ResultsComparison_new_h */
