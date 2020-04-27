#include "TH1F.h"
#include <cmath>

double TH1_integral (TH1F *h, float xmin, float xmax){
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/ axis->GetBinWidth(bmax);

    return integral;
}

void entries_invmass()
{
  TString run[] = {"2017B", "2017C", "2017D", "2017E", "2017F"};
  TH1F *h_mass_14;
  TH1F *h_mass_15;
  TH1F *h_mass_16;

  for(int i = 0; i<5; i++){
       //file containing invariant mass plot
      TString filename_tau3mu = "AnalysedTree_data_"+run[i]+"_tau3mu_29oct.root";
      TFile *f_tau3mu = new TFile(filename_tau3mu,"READ");
      h_mass_14 = (TH1F*)f_tau3mu->Get("StepByStep/Triplet/Mass triplet_cut14");
      h_mass_15 = (TH1F*)f_tau3mu->Get("StepByStep/Triplet/Mass triplet_cut15");
      h_mass_16 = (TH1F*)f_tau3mu->Get("StepByStep/Triplet/Mass triplet_cut16");

      cout<<run[i]<<" Double "<<h_mass_14->GetEntries()<<endl;
      cout<<run[i]<<" Triple "<<h_mass_15->GetEntries()<<endl;
      Double_t SB_or = TH1_integral(h_mass_16, 1.65, 1.73) + TH1_integral(h_mass_16, 1.82, 1.90);
      cout<<run[i]<<" OR "<<SB_or<<endl;
  }
}
