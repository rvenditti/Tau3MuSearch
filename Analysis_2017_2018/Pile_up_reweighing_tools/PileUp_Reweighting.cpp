#include <stdio.h>

//root -l PileUp_Reweighting.cpp\(\"2017C\"\)

void PileUp_Reweighting(TString datasetName){
    
    TTree *tree;

    //TString filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20191116_1520/AnalysedTree_MC_Ds_tau3mu_16nov.root"; // MC_DsTau3Mu
    //TString filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20191116_1521/AnalysedTree_MC_B0_tau3mu_16nov.root"; // MC_B0Tau3Mu
    //TString filename1 = "/lustre/cms/store/user/fsimone/Tau3Mu/Analysis/20191116_1523/AnalysedTree_MC_Bp_tau3mu_16nov.root"; // MC_BpTau3Mu
    TString filename1 = "AnalysedTree_2global_MC_DsPhiPi_control_09oct.root"; // MC_DsPhiPi
    TString filename2 = "MyDataPileupHistogram_";
    filename2 += datasetName;   filename2 += ".root"; // Data
    
    TFile *f1 = new TFile(filename1);
    TFile *f2 = new TFile(filename2);
    
    // Creation of the output file
    TString foutName = "PileUp_ReweightingStudy_";
    foutName += datasetName; foutName += ".root";
    TFile *fout = new TFile(foutName, "RECREATE");
    fout->cd();
    
    
    TH1 *hPileUp_MC, *hPileUp_data_temp;
    double weight = 0, scale1 = 0, scale2 = 0;
    // MC pile-up histo
    hPileUp_MC = (TH1*)f1->Get("BeforeCuts/hNPileUp");
    if(hPileUp_MC->GetEntries() != 0)
        scale1 = 1/(hPileUp_MC->Integral());
    else scale1 = 1;
    cout << "Scale factor MC : " << scale1 << endl;
    hPileUp_MC->Scale(scale1);
    // Data pile-up histo
    hPileUp_data_temp = (TH1*)f2->Get("pileup");
    TH1 *hPileUp_data = new TH1D("Pileup distributions", "Pileup distributions", 80, -0.5, 79.5);
    for(int m=1; m<hPileUp_MC->GetXaxis()->GetXmax(); m++){
        hPileUp_data->SetBinContent(m+1, hPileUp_data_temp->GetBinContent(m));
    }
        
    if(hPileUp_data->GetEntries() != 0)
        scale2 = 1/(hPileUp_data->Integral());
    else scale2 = 1;
    cout << "Scale factor data : " << scale2 << endl;
    hPileUp_data->Scale(scale2);
    
    int nBins1 = hPileUp_MC->GetNbinsX();
    cout << "MC pile-up histo has " << nBins1 << " bins." << endl;
    int nBins2 = hPileUp_data->GetNbinsX();
    cout << "Data pile-up histo has " << nBins2 << " bins." << endl;
    double Xmin = hPileUp_MC->GetXaxis()->GetXmin();
    cout << "Xmin : " << Xmin << endl;
    double Xmax = hPileUp_MC->GetXaxis()->GetXmax();
    cout << "Xmax : " << Xmax << endl;
    
    double val1[100] = {0}, val2[100] = {0};
    for(int k=0; k<nBins1; k++){
        val1[k] = hPileUp_MC->GetBinContent(k);
        val2[k] = hPileUp_data->GetBinContent(k);
    }
    
    // Canvas with PileUp distributions normalized to 1
    TString canvasName = "PileUp_distr_"; canvasName += datasetName;
    TCanvas *canv = new TCanvas(canvasName, canvasName, 0, 0, 1200, 1000);
    hPileUp_data->SetStats(0);
    hPileUp_data->SetMarkerStyle(22);
    hPileUp_data->SetMarkerColor(1);
    hPileUp_data->GetXaxis()->SetTitle("nPileUp");
    hPileUp_data->GetYaxis()->SetTitle("Entries");
    hPileUp_data->Draw();

    hPileUp_MC->SetLineColor(kBlue+1);
    hPileUp_MC->SetLineWidth(2);
    hPileUp_MC->SetFillColorAlpha(kBlue, 0.35);
    hPileUp_MC->Draw("bar same");
    
    TLegend *leg = new TLegend(0.80,0.80,0.99,0.94);
    TString legLabel[2];
    legLabel[0] = "Pileup_MC";
    legLabel[1] = "Pileup_data2017"; //legLabel[0] += datasetName;
    leg->AddEntry(hPileUp_MC, legLabel[0], "f");
    leg->AddEntry(hPileUp_data, legLabel[1], "p");
    leg->Draw();
    
    canv->Write();
    canv->Close();
    
    //Pile-up Reweighting
    TH1D *hweight = new TH1D("PileUp_Reweighting", "PileUp_Reweighting", 80, -0.5, 79.5);
    for(int i=1; i<81; i++){
        cout << " BIN n. " << i << endl;
        cout << " val_MC : " << val1[i] << endl;
        cout << " val_data : " << val2[i] << endl;
        if(val1[i] == 0 || val2[i] < 1e-7) weight = 0;
        else weight = val2[i]/val1[i];
        cout << "weight : " << weight << endl << endl;
        hweight->SetBinContent(i, weight);
    }
    
     // Write and close the file
    fout->Write();
//    fout->Close();
    f1->Close();
    f2->Close();
    
    
}
