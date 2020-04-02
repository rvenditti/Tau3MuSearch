#include "TH1F.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"

void plot_BDT_inputvariables() 
{
    TString cat[3] = {"A","B","C"};
    TString var[] = {"Pmu3","cLP","tKink","segmComp","fv_nC","fv_dphi3D","fv_d3Dsig","d0sig","mindca_iso","trkRel","nMatchesMu3"};
    cout<<"Opening data file"<<endl;
    TFile *f_data = new TFile(inputpath_data,"READ");
    cout<<"Opened data file"<<endl;
    cout<<"Opening MC file"<<endl;
    TFile *f_mc = new TFile(inputpath_MC,"READ");
    cout<<"Opened MC file"<<endl;
    
    TTree *tdata_A = (TTree*)f_data->Get("FinalTreeA_sgn");
    TTree *tdata_B = (TTree*)f_data->Get("FinalTreeB_sgn");
    TTree *tdata_C = (TTree*)f_data->Get("FinalTreeC_sgn");

    TTree *tmc_A = (TTree*)f_mc->Get("FinalTreeA_sgn");
    TTree *tmc_B = (TTree*)f_mc->Get("FinalTreeB_sgn");
    TTree *tmc_C = (TTree*)f_mc->Get("FinalTreeC_sgn");

    int n = sizeof(var)/sizeof(var[0]);
    TH1F *hdata_A[n];
    TH1F *hdata_B[n];
    TH1F *hdata_C[n];

    TH1F *hmc_A[n];
    TH1F *hmc_B[n];
    TH1F *hmc_C[n];
    
    TString binning;
    TString invmass_selection = "(tripletMass<1.75 && tripletMass>1.62) || (tripletMass<2.0 && tripletMass>1.80)";
 
   for(int i = 0; i<n; i++){
        TString varname = var[i];

        TString s = std::to_string(i);
        cout<<"Input variable "<<varname<<endl;

        TString binning;
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
        tdata_A->Draw(varname+">>hdata_A"+s+binning, invmass_selection);
        tdata_B->Draw(varname+">>hdata_B"+s+binning, invmass_selection);
        //if(varname=="mindca_iso")
        //    tdata_C->Draw(varname+">>hdata_C"+s+binning, "mindca_iso<10");
        //else if(varname=="trkRel")
        //    tdata_C->Draw(varname+">>hdata_C"+s+binning, "trkRel<20");
        //else
            tdata_C->Draw(varname+">>hdata_C"+s+binning, invmass_selection);

        hdata_A[i] = (TH1F *)gDirectory->Get("hdata_A"+s); 
        hdata_B[i] = (TH1F *)gDirectory->Get("hdata_B"+s); 
        hdata_C[i] = (TH1F *)gDirectory->Get("hdata_C"+s); 

        tmc_A->Draw(varname+">>hmc_A"+s+binning);
        tmc_B->Draw(varname+">>hmc_B"+s+binning);
        //if(varname=="mindca_iso")
        //    tmc_C->Draw(varname+">>hmc_C"+s+binning, "mindca_iso<10");
        //else if(varname=="trkRel")
        //    tmc_C->Draw(varname+">>hmc_C"+s+binning, "trkRel<20");
        //else
        tmc_C->Draw(varname+">>hmc_C"+s+binning);

        hmc_A[i] = (TH1F *)gDirectory->Get("hmc_A"+s); 
        hmc_B[i] = (TH1F *)gDirectory->Get("hmc_B"+s); 
        hmc_C[i] = (TH1F *)gDirectory->Get("hmc_C"+s); 

        //Data
        TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        hdata_A[i]->SetLineColor(kBlue);
        hdata_B[i]->SetLineColor(kRed);
        hdata_C[i]->SetLineColor(kBlack);
        cout<<"A max "<<hdata_A[i]->GetXaxis()->GetXmax()<<" min "<<hdata_A[i]->GetXaxis()->GetXmin()<<" mean "<<hdata_A[i]->GetMean()<<" RMS "<<hdata_A[i]->GetRMS()<<endl;
        cout<<"B max "<<hdata_B[i]->GetXaxis()->GetXmax()<<" min "<<hdata_B[i]->GetXaxis()->GetXmin()<<" mean "<<hdata_B[i]->GetMean()<<" RMS "<<hdata_B[i]->GetRMS()<<endl;
        cout<<"C max "<<hdata_C[i]->GetXaxis()->GetXmax()<<" min "<<hdata_C[i]->GetXaxis()->GetXmin()<<" mean "<<hdata_C[i]->GetMean()<<" RMS "<<hdata_C[i]->GetRMS()<<endl;
        
        Double_t min_data = std::min( std::min(hdata_A[i]->GetXaxis()->GetXmin(),hdata_B[i]->GetXaxis()->GetXmin()) , hdata_C[i]->GetXaxis()->GetXmin() );
        min_data = min_data - min_data*0.1;
        Double_t max_data = std::max( std::max(hdata_A[i]->GetXaxis()->GetXmax(),hdata_B[i]->GetXaxis()->GetXmax()) , hdata_C[i]->GetXaxis()->GetXmax() );
        max_data = max_data + max_data*0.1;
        if(varname == "mindca_iso") max_data = 4;
        cout<<"max_data "<<max_data<<endl; 
        hdata_A[i]->GetXaxis()->SetRangeUser(min_data, max_data);
        hdata_A[i]->GetXaxis()->SetTitle(varname);
        //hdata_A[i]->GetXaxis()->SetLimits(min_data, max_data);

        hdata_A[i]->Draw("");
        hdata_B[i]->Draw("same");
        hdata_C[i]->Draw("same");
 
        TLegend*leg = new TLegend(0.4,0.7,0.7,0.9);
        leg->AddEntry(hdata_A[i],"2018 "+varname+" dataA","f");
        leg->AddEntry(hdata_B[i],"2018 "+varname+" dataB","f");
        leg->AddEntry(hdata_C[i],"2018 "+varname+" dataC","f");
        leg->Draw();

        c1->SetLogy();
        c1->Update();
        c1->SaveAs(TMVA_inputpath+"A/Data_"+varname+".png");

        //MC
        TCanvas *c2 = new TCanvas("c2","c2",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        hmc_A[i]->SetLineColor(kBlue);
        hmc_B[i]->SetLineColor(kRed);
        hmc_C[i]->SetLineColor(kBlack);
        cout<<"A max "<<hmc_A[i]->GetXaxis()->GetXmax()<<" min "<<hmc_A[i]->GetXaxis()->GetXmin()<<" mean "<<hmc_A[i]->GetMean()<<" RMS "<<hmc_A[i]->GetRMS()<<endl;
        cout<<"B max "<<hmc_B[i]->GetXaxis()->GetXmax()<<" min "<<hmc_B[i]->GetXaxis()->GetXmin()<<" mean "<<hmc_B[i]->GetMean()<<" RMS "<<hmc_B[i]->GetRMS()<<endl;
        cout<<"C max "<<hmc_C[i]->GetXaxis()->GetXmax()<<" min "<<hmc_C[i]->GetXaxis()->GetXmin()<<" mean "<<hmc_C[i]->GetMean()<<" RMS "<<hmc_C[i]->GetRMS()<<endl;
        
        Double_t min_mc = std::min( std::max(hmc_A[i]->GetXaxis()->GetXmin(),hmc_B[i]->GetXaxis()->GetXmin()) , hmc_C[i]->GetXaxis()->GetXmin() );
        min_mc = min_mc - min_mc*0.1;
        Double_t max_mc = std::max( std::min(hmc_A[i]->GetXaxis()->GetXmax(),hmc_B[i]->GetXaxis()->GetXmax()) , hmc_C[i]->GetXaxis()->GetXmax() );
        max_mc = max_mc + max_mc*0.1;
        cout<<"max_mc "<<max_mc<<endl; 
        hmc_B[i]->GetXaxis()->SetRangeUser(min_mc, max_mc);
        hmc_B[i]->GetXaxis()->SetTitle(varname);
        //hmc_A[i]->GetXaxis()->SetLimits(min_mc, max_mc);
 
        hmc_B[i]->Draw("");
        hmc_A[i]->Draw("same");
        hmc_C[i]->Draw("same");

//        if(i==5) hmc_A[i]->GetXaxis()->SetRangeUser(-0.1, 5);
        TLegend*leg2 = new TLegend(0.4,0.7,0.7,0.9);
        leg2->AddEntry(hmc_A[i],"2018 "+varname+" mcA","f");
        leg2->AddEntry(hmc_B[i],"2018 "+varname+" mcB","f");
        leg2->AddEntry(hmc_C[i],"2018 "+varname+" mcC","f");
        leg2->Draw();

        c2->SetLogy();
        c2->Update();
        c2->SaveAs(TMVA_inputpath+"A/MC_"+varname+".png");
        
    }
}
