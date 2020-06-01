#include "TH1F.h"
#include <cmath>
#include <string> 
#include "../T3M_common.h"

void plot_BDT_inputvariables() 
{

    TString cat_label[] = {"A", "B", "C"};
    TString var[] = {
                     "Pmu3","cLP","tKink","segmComp","fv_nC","fv_dphi3D","fv_d3Dsig","d0sig","mindca_iso","trkRel","nMatchesMu3",
                     "fv_d3D", "abs(dxy1/dxyErr1)", "abs(dxy2/dxyErr2)", "abs(dxy3/dxyErr3)",
                     "TreeMu1.mu_SoftMVA",
                     "TreeMu2.mu_SoftMVA",
                     "TreeMu3.mu_SoftMVA",
                     "MuonIDeval_Mu1.MuonID",
                     "MuonIDeval_Mu2.MuonID",
                     "MuonIDeval_Mu3.MuonID"
                     };

    int nrun = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    //open input files
    //data
    TChain *tdata = new TChain("FinalTree");
    for(auto i=0; i<nrun; i++){
        tdata->Add(inputpath_datarun[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC Ds, B0, Bp
    TChain *tmc = new TChain("FinalTree");
    tmc->Add(inputpath_Ds); 
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    tmc->Add(inputpath_B0); 
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    tmc->Add(inputpath_Bp); 
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;

    //data
    TChain *tdata_mu1 = new TChain("TreeMu1");
    TChain *tdata_mu2 = new TChain("TreeMu2");
    TChain *tdata_mu3 = new TChain("TreeMu3");
    for(auto i=0; i<nrun; i++){
        tdata_mu1->Add(inputpath_datarun[i]);
        tdata_mu2->Add(inputpath_datarun[i]);
        tdata_mu3->Add(inputpath_datarun[i]);
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC Ds, B0, Bp
    TChain *tmc_mu1 = new TChain("TreeMu1");
    tmc_mu1->Add(inputpath_Ds);
    tmc_mu1->Add(inputpath_B0);
    tmc_mu1->Add(inputpath_Bp);
    TChain *tmc_mu2 = new TChain("TreeMu2");
    tmc_mu2->Add(inputpath_Ds);
    tmc_mu2->Add(inputpath_B0);
    tmc_mu2->Add(inputpath_Bp);
    TChain *tmc_mu3 = new TChain("TreeMu3");
    tmc_mu3->Add(inputpath_Ds);
    tmc_mu3->Add(inputpath_B0);
    tmc_mu3->Add(inputpath_Bp);

    //data
    TChain *tdata_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tdata_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tdata_muid3 = new TChain("MuonIDeval_Mu3");
    for(auto i=0; i<nrun; i++){
        TString inputpath_data_muId = inputpath_datarun[i].ReplaceAll(".root", "_MuonID.root");
        tdata_muid1->Add(inputpath_datarun[i]);
        tdata_muid2->Add(inputpath_datarun[i]);
        tdata_muid3->Add(inputpath_datarun[i]);
        std::cout<<"Opened input file: "<<inputpath_data_muId<<std::endl;
    }
    //MC Ds, B0, Bp
    TString inputpath_Ds_muId = inputpath_Ds.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_B0_muId = inputpath_B0.ReplaceAll(".root", "_MuonID.root");
    TString inputpath_Bp_muId = inputpath_Bp.ReplaceAll(".root", "_MuonID.root");
    TChain *tmc_muid1 = new TChain("MuonIDeval_Mu1");
    tmc_muid1->Add(inputpath_Ds_muId);
    tmc_muid1->Add(inputpath_B0_muId);
    tmc_muid1->Add(inputpath_Bp_muId);
    TChain *tmc_muid2 = new TChain("MuonIDeval_Mu2");
    tmc_muid2->Add(inputpath_Ds_muId);
    tmc_muid2->Add(inputpath_B0_muId);
    tmc_muid2->Add(inputpath_Bp_muId);
    TChain *tmc_muid3 = new TChain("MuonIDeval_Mu3");
    tmc_muid3->Add(inputpath_Ds_muId);
    tmc_muid3->Add(inputpath_B0_muId);
    tmc_muid3->Add(inputpath_Bp_muId);
    std::cout<<"Opened input file: "<<inputpath_Ds_muId<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_B0_muId<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_Bp_muId<<std::endl;

    tdata->AddFriend(tdata_mu1);
    tdata->AddFriend(tdata_mu2);
    tdata->AddFriend(tdata_mu3);
    tdata->AddFriend(tdata_muid1);
    tdata->AddFriend(tdata_muid2);
    tdata->AddFriend(tdata_muid3);

    tmc->AddFriend(tmc_mu1);
    tmc->AddFriend(tmc_mu2);
    tmc->AddFriend(tmc_mu3);
    tmc->AddFriend(tmc_muid1);
    tmc->AddFriend(tmc_muid2);
    tmc->AddFriend(tmc_muid3);

    int n = sizeof(var)/sizeof(var[0]);
    TH1F *hdata[3];
    TH1F *hmc[3];
    TH1F *hdata_merged;
    TH1F *hmc_merged;

    TString binning;

    TCut cutS = "tripletMass<2.0 && tripletMass>1.62"; //Signal -> MC full range 
    TCut cutB = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)"; //Background -> data sidebands

    TCut reso_A = "tripletMassReso < 0.007";
    TCut reso_B = "tripletMassReso >= 0.007 && tripletMassReso <= 0.0105";
    TCut reso_C = "tripletMassReso > 0.0105";
    TCut reso_cat = "tripletMassReso < 0"; //always false 

    //Loop on variables
    for(int i = 0; i<n; i++){
        TString varname = var[i];
        cout<<"Input variable "<<varname<<endl;
        binning = "";
        if(varname=="Pmu3") binning = "(250,0,50)";
        if(varname=="cLP") binning = "(250,0,1000)";
        if(varname=="segmComp") binning = "(100,-0.1,1.1)";
        if(varname=="tKink") binning = "(250,0,2000)";
        if(varname=="fv_nC") binning = "(250,-0.5,50.5)";
        if(varname=="fv_dphi3D") binning = "(100,-0.1,2)";
        if(varname=="fv_d3Dsig") binning = "(250,-0.1,1000)";
        if(varname=="d0sig") binning = "(150,-0.1,65)";
        if(varname=="mindca_iso") binning = "(150,-0.1,2)";
        if(varname=="trkRel") binning = "(150,-0.5,50)";
        if(varname=="nMatchesMu3") binning = "(20,0,10)";
        if(varname=="fv_d3D") binning = "(25,0,10)";
        if(varname.Contains("mu_SoftMVA")) binning = "(200,-1,1)";
        if(varname.Contains("MuonID"))     binning = "(200,-1,1)";
 
        //Loop on categories A, B, C
        for(auto j = 0; j<3; j++){
            TString s = std::to_string(j);
            TString category = cat_label[j];
            cout<<"Category "<<category<<endl;
            if(category == "A") reso_cat = reso_cat || reso_A;
            if(category == "B") reso_cat = reso_cat || reso_B;
            if(category == "C") reso_cat = reso_cat || reso_C;

            tdata->Draw(varname+">>hdata"+s+binning, cutB&&reso_cat);
            hdata[j] = (TH1F *)gDirectory->Get("hdata"+s);

            tmc->Draw(varname+">>hmc"+s+binning, cutS&&reso_cat);
            hmc[j] = (TH1F *)gDirectory->Get("hmc"+s);

            reso_cat = "tripletMassReso < 0";
        }
        tdata->Draw(varname+">>hdata_merged"+binning, cutB);
        hdata_merged = (TH1F *)gDirectory->Get("hdata_merged");

        tmc->Draw(varname+">>hmc_merged"+binning, cutS);
        hmc_merged = (TH1F *)gDirectory->Get("hmc_merged");
 
        //plot categories on same canvas - Data
        TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        hdata[0]->SetLineColor(kBlue);
        hdata[1]->SetLineColor(kRed);
        hdata[2]->SetLineColor(kBlack);
    
        Double_t min_data = std::min( std::min(hdata[0]->GetXaxis()->GetXmin(),hdata[1]->GetXaxis()->GetXmin()) , hdata[2]->GetXaxis()->GetXmin() );
        min_data = min_data - min_data*0.1;
        Double_t max_data = std::max( std::max(hdata[0]->GetXaxis()->GetXmax(),hdata[1]->GetXaxis()->GetXmax()) , hdata[2]->GetXaxis()->GetXmax() );
        max_data = max_data + max_data*0.1;
        if(varname == "mindca_iso") max_data = 4;
    
        hdata[0]->GetXaxis()->SetRangeUser(min_data, max_data);
        hdata[0]->GetXaxis()->SetTitle(varname);
    
        hdata[0]->Draw("");
        hdata[1]->Draw("same");
        hdata[2]->Draw("same");
    
        TLegend*leg = new TLegend(0.4,0.7,0.7,0.9);
        leg->AddEntry(hdata[0],"2018 "+varname+" data_A","f");
        leg->AddEntry(hdata[1],"2018 "+varname+" data_B","f");
        leg->AddEntry(hdata[2],"2018 "+varname+" data_C","f");
        leg->Draw();
    
        //c1->SetLogy();
        c1->Update();
        varname = varname.ReplaceAll("/","_");
        varname = varname.ReplaceAll("(","_");
        varname = varname.ReplaceAll(")","_");
        c1->SaveAs("../plots/"+TMVA_inputpath+"_Data_"+varname+".png");
    
    
        //plot categories on same canvas - MC
        TCanvas *c2 = new TCanvas("c2","c2",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        hmc[0]->SetLineColor(kBlue);
        hmc[1]->SetLineColor(kRed);
        hmc[2]->SetLineColor(kBlack);
    
        Double_t min_mc = std::min( std::min(hmc[0]->GetXaxis()->GetXmin(),hmc[1]->GetXaxis()->GetXmin()) , hmc[2]->GetXaxis()->GetXmin() );
        min_mc = min_mc - min_mc*0.1;
        Double_t max_mc = std::max( std::max(hmc[0]->GetXaxis()->GetXmax(),hmc[1]->GetXaxis()->GetXmax()) , hmc[2]->GetXaxis()->GetXmax() );
        max_mc = max_mc + max_mc*0.1;
        if(varname == "mindca_iso") max_mc = 4;
    
        hmc[0]->GetXaxis()->SetRangeUser(min_mc, max_mc);
        hmc[0]->GetXaxis()->SetTitle(varname);
    
        hmc[0]->Draw("");
        hmc[1]->Draw("same");
        hmc[2]->Draw("same");
    
        TLegend*leg2 = new TLegend(0.4,0.7,0.7,0.9);
        leg2->AddEntry(hmc[0],"2018 "+varname+" mc_A","f");
        leg2->AddEntry(hmc[1],"2018 "+varname+" mc_B","f");
        leg2->AddEntry(hmc[2],"2018 "+varname+" mc_C","f");
        leg2->Draw();
    
        //c2->SetLogy();
        c2->Update();
        c2->SaveAs("../plots/"+TMVA_inputpath+"_MC_"+varname+".png");

        TCanvas *c3 = new TCanvas("c3","c3",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        hmc_merged->SetLineColor(kBlue);
        hmc_merged->SetLineWidth(3);
        hmc_merged->SetFillStyle(3004);
        hmc_merged->SetFillColor(kBlue);
        hdata_merged->SetLineColor(kRed);
        hdata_merged->SetLineWidth(3);
        hdata_merged->SetFillStyle(3005);
        hdata_merged->SetFillColor(kRed);

        Double_t min_compare = std::min(min_data, min_mc);
        Double_t max_compare = std::max(max_data, max_mc);

        hdata_merged->GetXaxis()->SetRangeUser(min_compare, max_compare);
        hdata_merged->GetXaxis()->SetTitle(varname);
        hdata_merged->Scale(1.0/(hdata_merged->Integral()));

        hmc_merged->GetXaxis()->SetRangeUser(min_compare, max_compare);
        hmc_merged->GetXaxis()->SetTitle(varname);
        hmc_merged->Scale(1.0/(hmc_merged->Integral()));

        c3->Update();
        hmc_merged->Draw();
        hdata_merged->Draw("same");

        TLegend*leg3 = new TLegend(0.2,0.7,0.45,0.9);
        leg3->AddEntry(hmc_merged,"2018 "+varname+" mc","f");
        leg3->AddEntry(hdata_merged,"2018 "+varname+" data","f");
        leg3->Draw();

        //c3->SetLogy();
        c3->Update();
        c3->SaveAs("../plots/"+TMVA_inputpath+"_compare_"+varname+".png");

    }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}

