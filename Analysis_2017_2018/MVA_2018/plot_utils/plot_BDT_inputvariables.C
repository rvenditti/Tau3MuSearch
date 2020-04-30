#include "TH1F.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"

void plot_BDT_inputvariables() 
{

    TString cat_label[] = {"A", "B", "C"};
    TString var[] = {"Pmu3","cLP","tKink","segmComp","fv_nC","fv_dphi3D","fv_d3Dsig","d0sig","mindca_iso","trkRel","nMatchesMu3"};
    TString run[] = {"A", "B", "C", "D"};
    size_t nrun = sizeof(run)/sizeof(run[0]);
    //open input files
    //data
    TChain *tdata = new TChain("FinalTree");
    for(auto i=0; i<nrun; i++){
        tdata->Add(inputpath_datarun[i]);
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC 
    TChain *tmc = new TChain("FinalTree");
    tmc->Add(inputpath_Ds);
    tmc->Add(inputpath_B0);
    tmc->Add(inputpath_Bp);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;

    int n = sizeof(var)/sizeof(var[0]);
    TH1F *hdata[3];
    TH1F *hmc[3];

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
    
        c1->SetLogy();
        c1->Update();
        c1->SaveAs("plots/"+TMVA_inputpath+"_Data_"+varname+".png");
    
    
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
    
        c2->SetLogy();
        c2->Update();
        c2->SaveAs("plots/"+TMVA_inputpath+"_MC_"+varname+".png");

    }
    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}

