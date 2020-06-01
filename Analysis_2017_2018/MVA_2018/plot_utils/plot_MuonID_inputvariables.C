#include "TH1F.h"
#include <cmath>
#include <string> 
#include "../T3M_common.h"

void plot_MuonID_inputvariables(TString muon) 
{
    if(! (muon=="Mu1" || muon=="Mu2" || muon=="Mu3") ) {cout<<"Argument Mu1 Mu2 Mu3"<<endl; return;}

    TString var[]  = {
                                      //Muon reconstruction
                                     "mu_combinedQuality_chi2LocalMomentum",
                                     "mu_combinedQuality_chi2LocalPosition",
                                     "mu_combinedQuality_staRelChi2",
                                     "mu_combinedQuality_trkRelChi2",
                                     "mu_combinedQuality_globalDeltaEtaPhi",
                                     "mu_combinedQuality_trkKink", //"log_mu_combinedQuality_trkKink",
                                     "mu_combinedQuality_glbKink", //"log_mu_combinedQuality_glbKink",
                                     "mu_combinedQuality_glbTrackProbability",
                                   
                                      //collection of hits in the HitPattern
                                     "mu_Numberofvalidtrackerhits", //Valid Tracker Hits
                                     "mu_Numberofvalidpixelhits",
                                     "mu_trackerLayersWithMeasurement",
                                     "mu_GLhitPattern_numberOfValidMuonHits",
                                     "mu_validMuonHitComb", //Hits in DT, CSC, RPC
                                   
                                      //muon track reconstruction
                                     "mu_numberOfMatchedStations",
                                     "mu_segmentCompatibility",
                                     "mu_timeAtIpInOut",
                                     "mu_timeAtIpInOutErr",
                                   
                                     //general track properties
                                      "mu_GLnormChi2",
                                      "mu_innerTrack_normalizedChi2",
                                      "mu_outerTrack_normalizedChi2",
                                      "mu_innerTrack_validFraction", //Inner Valid Fraction
                                      "mu_QInnerOuter",
                                      //"", //dxyRef
                                      //"", //dzRef
                                      
                                      //custom variables track multiplicity

                                      //isolation variables
                                      "mu_emEt03",
                                      "mu_hadEt03",
                                      "mu_nJets03",
                                      "mu_nTracks03",
                                      "mu_sumPt03",
                                      "mu_hadVetoEt03", 
                                      "mu_emVetoEt03",
                                      "mu_trackerVetoPt03" 
                                      };

    size_t nrun = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    //open input files
    TString treename = "Tree"+muon;
    //data
    cout<<"Data"<<endl;
    TChain *tdata = new TChain(treename);
    for(auto i=0; i<nrun; i++){
        tdata->Add(inputpath_datarun[i]);
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC 
    cout<<"MC"<<endl;
    TChain *tmc = new TChain(treename);
    tmc->Add(inputpath_Ds);
    tmc->Add(inputpath_B0);
    tmc->Add(inputpath_Bp);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;

    int n = sizeof(var)/sizeof(var[0]);
    TH1F *hdata[n];
    TH1F *hmc[n];

    TString binning;

    //Loop on variables
    for(int i = 0; i<n; i++){
        TString varname = var[i];
        cout<<"Input variable "<<varname<<endl;
        TString s = std::to_string(i);

        binning = "";
        if(varname=="mu_combinedQuality_chi2LocalMomentum") binning = "(400,0,8000)";
        if(varname=="mu_combinedQuality_glbKink") binning = "(500,0,10000)";
        if(varname=="mu_GLnormChi2") binning = "(200,0,6000)";
        if(varname=="mu_emEt03") binning = "(100,0,300)";
        if(varname=="mu_hadEt03") binning = "(110,0,350)";
        if(varname=="mu_nJets03") binning = "(5,0,5)";
        if(varname=="mu_nTracks03") binning = "(40,0,40)";
        if(varname=="mu_sumPt03") binning = "(150,0,1000)";
        if(varname=="mu_emVetoEt03") binning = "(50,0,200)";
        if(varname=="mu_hadVetoEt03") binning = "(50,0,200)";
        if(varname=="mu_trackerVetoPt03") binning = "(60,0,300)";

        tdata->Draw(varname+">>hdata"+s+binning);
        hdata[i] = (TH1F *)gDirectory->Get("hdata"+s);

        tmc->Draw(varname+">>hmc"+s+binning);
        hmc[i] = (TH1F *)gDirectory->Get("hmc"+s);

 
        //plot categories on same canvas - Data
        TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        hdata[i]->SetLineColor(kRed);
        hmc[i]->SetLineColor(kBlue);
    
        Double_t min = std::min(hdata[i]->GetXaxis()->GetXmin(),hmc[i]->GetXaxis()->GetXmin() );
        min = min - min*0.1;
        Double_t max = std::max(hdata[i]->GetXaxis()->GetXmax(),hmc[i]->GetXaxis()->GetXmax() );
        max = max + max*0.1;
    
        //hdata[i]->GetXaxis()->SetRangeUser(min, max);
        hdata[i]->GetXaxis()->SetTitle(varname);
    
        hdata[i]->Draw("");
        hmc[i]->Draw("same");
    
        TLegend*leg = new TLegend(0.4,0.7,0.7,0.9);
        leg->AddEntry(hdata[i],"2018_"+varname+"_"+muon+"_data","f");
        leg->AddEntry(hmc[i],"2018_"+varname+"_"+muon+"_mc","f");
        leg->Draw();
    
        c1->SetLogy();
        c1->Update();
        c1->SaveAs("../plots/"+TMVA_inputpath+varname+muon+".png");

    }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}

