
#define Nhist_BC 17
#define Nhist_StepByStep 128
#define Nhist_AC 30
#define Nsets 3
#define NDs 86407
#define NBkg 80000

#include "ResultsComparison_Control.h"
#include <cmath>


void ResultsComparison_Control(TString outfile){
    
    double scale[Nsets] = {0}, ymax[Nsets] = {0}, ymax_fin = 0, ymax_fin2 = 0;
    TString dataName[Nsets], filename[Nsets];
    TFile *f[Nsets];
    outfile += ".root";
    
    // Input & Output files
    filename[0] = "../AnalysedTree/Control/AnalysedTree_2glb_DsPhiPi_Data2017B_C_D_E_F_sgn.root"; // Data DsPhiPi - signal (Ds) region [1.93-2.01]
    filename[1] = "../AnalysedTree/Control/AnalysedTree_2glb_DsPhiPi_Data2017B_C_D_E_F_bkg.root"; // Data DsPhiPi - bkg region [1.7-1.8]
    filename[2] = "../AnalysedTree/Control/AnalysedTree_2glb_MC_DsPhiPi.root"; // MC DsPhiPi
    int k=0;
    for(k=0; k<Nsets; k++){
        f[k] = new TFile(filename[k]);
    }
    TFile *fout = new TFile(outfile, "RECREATE");
    fout->cd();
    
    // Legend labels
    TString legendLabel[Nsets];
     legendLabel[0] = "Control data";
     legendLabel[2] = "Control MC DsPhiPi";
    
    //######################    Histo BEFORE cuts
    TDirectory *dirBeforeCuts = fout->mkdir("BeforeCuts");
    dirBeforeCuts->cd();
    // Creo e inizializzo cose utili per isto BC
    TString hname_BC[Nhist_BC], Xaxisname_BC[Nhist_BC];
    HistName_BC(hname_BC);
    XaxisName_BC(Xaxisname_BC);
    TH1 *h_BC[Nsets][Nhist_BC];
    TCanvas *canv_BC[Nhist_BC];
    TLegend *legend_BC[Nhist_BC];
    //  Get, draw & save histos BC
    for (int i=0; i<Nhist_BC; i++){
        for (k=0; k<Nsets; k++){
            h_BC[k][i] = (TH1*)f[k]->Get(hname_BC[i]);
            if(h_BC[k][i]->GetEntries() != 0)
                scale[k] = 1/(h_BC[k][i]->Integral());
            if (k==1)    h_BC[k][i]->Scale(NBkg * scale[k]);
            if (k==2){
                h_BC[k][i]->Scale(NDs/h_BC[k][i]->Integral());
                h_BC[k][i]->Rebin(4);
                h_BC[k][i]->Write();
            }
            ymax[k] =  (h_BC[k][i]->GetMaximum()) * 1.05;
            if (ymax[k] > ymax_fin) ymax_fin = ymax[k];
        }
        
        h_BC[0][i]->SetMarkerStyle(31);
        h_BC[0][i]->SetMarkerColor(1);
        h_BC[1][i]->SetMarkerStyle(31);
        h_BC[1][i]->SetMarkerColor(1);
        h_BC[2][i]->SetFillColorAlpha(kOrange+7, 0.5);
        h_BC[2][i]->SetLineColor(kOrange+7);
        
        h_BC[0][i]->Add(h_BC[1][i], -1);
        h_BC[0][i]->Rebin(4);
        h_BC[0][i]->Write();
        
        canv_BC[i] = new TCanvas(hname_BC[i], hname_BC[i], 0, 0, 1200, 1000);
        legend_BC[i] = new TLegend(0.80,0.80,0.99,0.94);
        
        legend_BC[i]->AddEntry(h_BC[2][i], legendLabel[2], "f");
        legend_BC[i]->AddEntry(h_BC[0][i], legendLabel[0], "p");
        
        h_BC[2][i]->SetStats(0);
        h_BC[2][i]->Draw("HIST");
        h_BC[0][i]->Draw("same");
        h_BC[2][i]->GetXaxis()->SetTitle(Xaxisname_BC[i]);
        h_BC[2][i]->GetYaxis()->SetTitle("Entries");
        
        legend_BC[i]->Draw();
//        canv_BC[i]->SetLogy();
        canv_BC[i]->Write();
        canv_BC[i]->Close();
    }
    fout->cd();
    
    //######################    Histo StepByStep
    TDirectory *dirInt = fout->mkdir("StepByStep");
    dirInt->cd();
    // Creo e inizializzo cose utili per isto intermediate
    TString hname_int[Nhist_StepByStep], Xaxisname_int[Nhist_StepByStep];
    HistName_StepByStep(hname_int);
    XaxisName_StepByStep(Xaxisname_int);
    TH1 *h_int[Nsets][Nhist_StepByStep];
    TCanvas *canv_int[Nhist_StepByStep];
    TLegend *legend_int[Nhist_StepByStep];
    
    //  Get, draw & save histos intermediate
    for (int i=0; i<Nhist_StepByStep; i++){
        for (k=0; k<Nsets; k++){
            h_int[k][i] = (TH1*)f[k]->Get(hname_int[i]);
            if(h_int[k][i]->GetEntries() != 0)
                scale[k] = 1/(h_int[k][i]->Integral());
            if (k==1)    {
//                h_int[k][i]->Write();
                h_int[k][i]->Scale(NBkg * scale[k]);
            }
            if (k==2){
//                h_int[k][i]->Write();
                h_int[k][i]->Scale(NDs/h_int[k][i]->Integral());
                h_int[k][i]->Rebin(4);
                
            }
        }
        
        h_int[0][i]->SetMarkerStyle(31);
        h_int[0][i]->SetMarkerColor(1);
        h_int[1][i]->SetMarkerStyle(31);
        h_int[1][i]->SetMarkerColor(1);
        h_int[2][i]->SetFillColorAlpha(kOrange+7, 0.5);
        h_int[2][i]->SetLineColor(kOrange+7);
        
        h_int[0][i]->Add(h_int[1][i], -1);
//        h_int[0][i]->Write();
        h_int[0][i]->Rebin(4);
        
        canv_int[i] = new TCanvas(hname_int[i], hname_int[i], 0, 0, 1200, 1000);
        legend_int[i] = new TLegend(0.80,0.80,0.99,0.94);
        
        legend_int[i]->AddEntry(h_int[2][i], legendLabel[2], "f");
        legend_int[i]->AddEntry(h_int[0][i], legendLabel[0], "p");
        
        h_int[2][i]->SetStats(0);
        h_int[2][i]->Draw("HIST");
        h_int[0][i]->Draw("same");
        h_int[2][i]->GetXaxis()->SetTitle(Xaxisname_int[i]);
        h_int[2][i]->GetYaxis()->SetTitle("Entries");
        
        legend_int[i]->Draw();
//        canv_int[i]->SetLogy();
        canv_int[i]->Write();
        canv_int[i]->Close();
    }
    fout->cd();
    
    //######################    Histo AFTER cuts
    TDirectory *dirAfterCuts = fout->mkdir("AfterCuts");
    dirAfterCuts->cd();
    // Creo e inizializzo cose utili per isto AC
    TString hname_AC[Nhist_AC], Xaxisname_AC[Nhist_AC];
    HistName_AC(hname_AC);
    XaxisName_AC(Xaxisname_AC);
    TH1 *h_AC[Nsets][Nhist_AC];
    TCanvas *canv_AC[Nhist_AC];
    TLegend *legend_AC[Nhist_AC];
    //  Get, draw & save histos AC
    for (int i=0; i<Nhist_AC; i++){
        for (k=0; k<Nsets; k++){
            h_AC[k][i] = (TH1*)f[k]->Get(hname_AC[i]);
            if(h_AC[k][i]->GetEntries() != 0)
                scale[k] = 1/(h_AC[k][i]->Integral());
            if (k==1)    h_AC[k][i]->Scale(NBkg * scale[k]);
            if (k==2){
                h_AC[k][i]->Scale(NDs/h_AC[k][i]->Integral());
                h_AC[k][i]->Rebin(4);
//                h_AC[k][i]->Write();
            }
        }
        
        h_AC[0][i]->SetMarkerStyle(31);
        h_AC[0][i]->SetMarkerColor(1);
        h_AC[1][i]->SetMarkerStyle(31);
        h_AC[1][i]->SetMarkerColor(1);
        h_AC[2][i]->SetFillColorAlpha(kOrange+7, 0.5);
        h_AC[2][i]->SetLineColor(kOrange+7);
        
        h_AC[0][i]->Add(h_AC[1][i], -1);
        h_AC[0][i]->Rebin(4);
        h_AC[0][i]->Write();
        
        canv_AC[i] = new TCanvas(hname_AC[i], hname_AC[i], 0, 0, 1200, 1000);
        legend_AC[i] = new TLegend(0.80,0.80,0.99,0.94);
        
        legend_AC[i]->AddEntry(h_AC[2][i], legendLabel[2], "f");
        legend_AC[i]->AddEntry(h_AC[0][i], legendLabel[0], "p");
        
        h_AC[2][i]->SetStats(0);
        h_AC[2][i]->Draw("HIST");
        h_AC[0][i]->Draw("same");
        h_AC[2][i]->GetXaxis()->SetTitle(Xaxisname_AC[i]);
        h_AC[2][i]->GetYaxis()->SetTitle("Entries");
        
        legend_AC[i]->Draw();
//        canv_AC[i]->SetLogy();
        canv_AC[i]->Write();
        canv_AC[i]->Close();
        
    }
    fout->cd();

    //Write and close the file
    fout->Write();
    fout->Close();
    for(k=0; k<Nsets; k++)
        f[k]->Close();
}
