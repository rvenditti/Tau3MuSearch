//root -l -b plot_TMVA_inputvariables.C\(\"C\"\)
#include "TH1F.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"

double TH1_integral(TH1F *h, float xmin, float xmax){
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/ axis->GetBinWidth(bmax);

    return integral;
}

void plot_TMVA_inputvariables(TString category) 
{
    TString TMVA_filename = TMVA_inputpath+category;
    TString var[] = {"Pmu3","cLP","tKink","segmComp","fv_nC","fv_dphi3D","fv_d3Dsig","d0sig","mindca_iso","trkRel", "nMatchesMu3"};
    cout<<"Opening data file"<<endl;
    TFile *f_tmva = new TFile("TMVA_"+TMVA_filename+".root","READ");
    cout<<"Opened TMVA file"<<endl;
    
    TTree *tTrain = (TTree*)f_tmva->Get(TMVA_filename+"/TrainTree");
    //TTree *tTest = (TTree*)f_tmva->Get(TMVA_filename+"/TestTree");

    int n = sizeof(var)/sizeof(var[0]);
    TH1F *hTrain_signal[n];
    TH1F *hTrain_bkg[n];

    TString binning;
 
   for(int i = 0; i<n; i++){
        TString varname = var[i];

        TString s = std::to_string(i);
        cout<<"Input variable "<<varname<<endl;

        TString binning;
        if(varname=="Pmu3") binning = "(40, 3, 45)"; //binning AN
        if(varname=="cLP") binning = "(40,0,30)"; //binning AN
        if(varname=="segmComp") binning = "(60,0,1)"; //binning AN
        if(varname=="tKink") binning = "(40,0,80)"; //binning AN
        if(varname=="fv_nC") binning = "(40, 0, 25)"; //binning AN
        if(varname=="fv_dphi3D") binning = "(40, 0, 0.16)"; //binning AN
        if(varname=="fv_d3Dsig") binning = "(40, 0, 100)"; //binning AN
        if(varname=="d0sig") binning = "(40, 0, 10)"; //binning AN
        if(varname=="mindca_iso") binning = "(40,0,0.5)"; //binning AN
        if(varname=="trkRel") binning = "(40,0,10)"; //binning AN 
        if(varname=="nMatchesMu3") binning = "(40,0,10)"; //binning AN 
        tTrain->Draw(varname+">>hTrain_signal"+s+binning, "classID==0");
        tTrain->Draw(varname+">>hTrain_bkg"+s+binning, "classID==1");

        hTrain_signal[i] = (TH1F *)gDirectory->Get("hTrain_signal"+s); 
        hTrain_bkg[i] = (TH1F *)gDirectory->Get("hTrain_bkg"+s); 

        //Data
        TCanvas *c1 = new TCanvas("c1","c1",150,10,990,660);
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(111111);
        hTrain_signal[i]->SetLineColor(kBlue);
        hTrain_signal[i]->SetLineWidth(3);
        hTrain_signal[i]->SetFillStyle(3004);
        hTrain_signal[i]->SetFillColor(kBlue);
        hTrain_bkg[i]->SetLineColor(kRed);
        hTrain_bkg[i]->SetLineWidth(3);
        hTrain_bkg[i]->SetFillStyle(3005);
        hTrain_bkg[i]->SetFillColor(kRed);
       
        Double_t norm = 1; 
        double X_min = std::min(hTrain_signal[i]->GetXaxis()->GetXmin(), hTrain_bkg[i]->GetXaxis()->GetXmin());
        double X_max = std::max(hTrain_signal[i]->GetXaxis()->GetXmax(), hTrain_bkg[i]->GetXaxis()->GetXmax());
        hTrain_signal[i]->Scale(norm / TH1_integral(hTrain_signal[i], X_min, X_max));
        hTrain_bkg[i]->Scale(norm / TH1_integral(hTrain_bkg[i], X_min, X_max));

        hTrain_signal[i]->GetXaxis()->SetTitle(varname);

        hTrain_signal[i]->GetXaxis()->SetRange(1, hTrain_signal[i]->GetNbinsX() + 1); // will draw with the overflow bin
        hTrain_bkg[i]->GetXaxis()->SetRange(1, hTrain_bkg[i]->GetNbinsX() + 1);    // will draw with the overflow bin

        ////Adjust Y range
        //if(varname=="Pmu3")      hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.17);
        //if(varname=="cLP")       hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.275); 
        //if(varname=="segmComp")  hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 3.6);
        //if(varname=="tKink")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.11);
        //if(varname=="fv_nC")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.55);
        //if(varname=="fv_dphi3D") hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 50);
        //if(varname=="fv_d3Dsig") hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.3);
        //if(varname=="d0sig")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 0.8);
        //if(varname=="mindca_iso") hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 31);
        //if(varname=="trkRel")     hTrain_signal[i]->GetYaxis()->SetRangeUser(0, 1.5);

        THStack hs("hs",varname);
        hs.Add(hTrain_signal[i]);
        hs.Add(hTrain_bkg[i]);
        hs.Draw("hist nostack");

        //hTrain_signal[i]->Draw("hist");
        //hTrain_bkg[i]->Draw("same hist");
 
        TLegend*leg = new TLegend(0.4,0.7,0.7,0.9);
        leg->AddEntry(hTrain_signal[i],varname+"_signal","f");
        leg->AddEntry(hTrain_bkg[i],varname+"_bkg","f");
        leg->Draw();

        //c1->SetLogy();
        c1->Update();
        c1->SaveAs(TMVA_inputpath+category+"/TMVA_train_"+varname+".png");

    }
}
