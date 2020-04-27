#include <stdio.h>
#include "../T3M_common.h"

//root -l PT_reweighting.cpp\(\"barrel\"\)

void PT_reweighting(TString categ){

    //Check on input argument
    if(!categ.Contains("barrel") && !categ.Contains("endcap")){
        cout << "Please choose between any combination of 'barrel' and 'endcap'" << endl;
        return;
    }
    
    //open input files
    size_t n_bkg = sizeof(inputpath_bkg)/sizeof(inputpath_bkg[0]);
    //Bkg samples
    TChain *tbkg = new TChain("FinalTree");
    for(auto i=0; i<n_bkg; i++){
        tbkg->Add(inputpath_bkg[i]);
        std::cout<<"Opened input file: "<<inputpath_bkg[i]<<std::endl;
    }
    //MC DsTau3Mu and B0Tau3Mu
    TChain *tsgn = new TChain("FinalTree");
    tsgn->Add(inputpath_Ds);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    tsgn->Add(inputpath_B0);
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;

    //Make sure TChain points to firts event
    tbkg->LoadTree(0);
    tsgn->LoadTree(0);

    //Book output histogram or MiniTree
    TString fout_path = "PT_reweighting_muonid_"+categ+".root";
    TFile *fout = new TFile(fout_path, "recreate");
    fout->cd();

    //signal muons
    TCut cutS = "abs(mu_simPdgId)==13 && abs(mu_simMotherPdgId)==15"; 
    //bkg muons = pions, kaon, muons frmo decay in flight
    TCut cutB_pi = "( abs(mu_simPdgId) == 211 )"; //pion
    TCut cutB_k = "( abs(mu_simPdgId) == 321 )"; //kaon
    TCut cutB_mu_from_pi = "( abs(mu_simPdgId) == 13 && abs(mu_simMotherPdgId) == 211 )"; //true mu from pion decay
    TCut cutB_mu_from_k = "( abs(mu_simPdgId) == 13 && abs(mu_simMotherPdgId) == 321 )"; //true mu from kaon decay
    //running on global muons with associated simInfo
    TCut preselCut = "mu_simType != 0 && mu_isGlobal == 1 && mu_pt > 2 && abs(mu_eta)<2.4";
    //training separately on endcap and barrel
    TCut etarange = "abs(mu_eta)<0"; //always false
    TCut barrel = "abs(mu_eta)<1.2";
    TCut endcap = "abs(mu_eta)>=1.2";
    if(categ.Contains("endcap")) etarange = etarange || endcap;
    if(categ.Contains("barrel")) etarange = etarange || barrel;
 
    //Take muon PT histograms 
    TH1D *hbkg_pt, *hsgn_pt;
    TString varname = "mu_pt";
    TString binning = "(125,0,50)";
    tbkg->Draw(varname+">>hbkg_pt"+binning, preselCut && etarange && (cutB_pi||cutB_k||cutB_mu_from_pi||cutB_mu_from_k));
    tsgn->Draw(varname+">>hsgn_pt"+binning, preselCut && etarange && cutS );

    hbkg_pt = (TH1D *)gDirectory->Get("hbkg_pt");
    hsgn_pt = (TH1D *)gDirectory->Get("hsgn_pt");
    
    //Normalise distributions
    Double_t scale_bkg = hbkg_pt->Integral();
    hbkg_pt->Scale(1/scale_bkg);
    Double_t scale_sgn = hsgn_pt->Integral();
    hsgn_pt->Scale(1/scale_sgn);
    
    //Compute ratio
    TH1D *hweight_pt = (TH1D*)hsgn_pt->Clone("hweight_pt");
    hweight_pt->Sumw2();
    hweight_pt->Divide(hbkg_pt);


    //Take muon ETA histograms 
    TH1D *hbkg_eta, *hsgn_eta;
    varname = "mu_eta";
    binning = "(120,-2.4,2.4)";
    tbkg->Draw(varname+">>hbkg_eta"+binning, preselCut && etarange && (cutB_pi||cutB_k||cutB_mu_from_pi||cutB_mu_from_k));
    tsgn->Draw(varname+">>hsgn_eta"+binning, preselCut && etarange && cutS );

    hbkg_eta = (TH1D *)gDirectory->Get("hbkg_eta");
    hsgn_eta = (TH1D *)gDirectory->Get("hsgn_eta");
    
    //Normalise distributions
    scale_bkg = hbkg_eta->Integral();
    hbkg_eta->Scale(1/scale_bkg);
    scale_sgn = hsgn_eta->Integral();
    hsgn_eta->Scale(1/scale_sgn);
    
    //Compute ratio
    TH1D *hweight_eta = (TH1D*)hsgn_eta->Clone("hweight_eta");
    hweight_eta->Sumw2();
    hweight_eta->Divide(hbkg_eta);
    
    //Take muon PT-ETA 2D histograms
    TH2D *hbkg_pteta = new TH2D("hbkg_pteta","hbkg_pteta",250,0,50,240,-2.4,2.4);
    TH2D *hsgn_pteta = new TH2D("hsgn_pteta","hsgn_pteta",250,0,50,240,-2.4,2.4);
    varname = "mu_eta:mu_pt";
    tbkg->Draw(varname+">>hbkg_pteta", preselCut && etarange && (cutB_pi||cutB_k||cutB_mu_from_pi||cutB_mu_from_k));
    tsgn->Draw(varname+">>hsgn_pteta", preselCut && etarange && cutS );

    //Normalise distributions
    scale_bkg = hbkg_pteta->Integral();
    hbkg_pteta->Scale(1/scale_bkg);
    scale_sgn = hsgn_pteta->Integral();
    hsgn_pteta->Scale(1/scale_sgn);
    
    //Compute ratio
    TH2D *hweight_pteta = (TH2D*)hsgn_pteta->Clone("hweight_pteta");
    hweight_pteta->Sumw2();
    hweight_pteta->Divide(hbkg_pteta);
     // Write and close the file
    fout->Write();

}
