#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

double TH1_integral(TH1F *h, float xmin, float xmax){
    TAxis *axis = h->GetXaxis();
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    double integral = h->Integral(bmin,bmax);
    integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/ axis->GetBinWidth(bmax);

    return integral;
}

void count_datacard() 
{
    TString category[] = {"A", "B", "C"};
    TString fv_nC_cut[] = {
                     "fv_nC<10000",
                     "fv_nC<5000",
                     "fv_nC<2500",
                     "fv_nC<1000",
                     "fv_nC<500",
                     "fv_nC<100",
                     "fv_nC<50",
                     "fv_nC<40",
                     "fv_nC<30",
                     "fv_nC<20",
                     "fv_nC<15",
                     "fv_nC<10",
                     "fv_nC<8",
                     "fv_nC<6",
                     "fv_nC<4",
    };
    size_t n = sizeof(fv_nC_cut)/sizeof(fv_nC_cut[0]);
    float sensitivity1[n];
    float sensitivity2[n];
    float signal1, signal2, data1, data2;

    //count reparately for each category
    for(auto j = 0; j<3; j++){
        cout<<"category "<<category[j]<<endl;
        cout<<"fv_nC_cut\tS1/sqrt(B1)\tB1\tS1\tS2/sqrt(B2)\tB2\tS2"<<endl; 
        //for each fv_nC_cut i open datacard and compute the sensitivity after BDT
        for(auto i = 0; i<n; i++){

            TString filename = "datacardT3Mu_dataset_2018_27april_Chi2_"+fv_nC_cut[i]+"_"+category[j]+".root";
            TFile *f = TFile::Open(filename, "READ");
            if (!f || !f->IsOpen()) std::cout<<"Error"<<endl;

            TH1F *hdata1 = (TH1F*)f->Get("data_obs"+category[j]+"1");
            TH1F *hdata2 = (TH1F*)f->Get("data_obs"+category[j]+"2");
            TH1F *hsgn1 =  (TH1F*)f->Get("signal"+category[j]+"1");
            TH1F *hsgn2 =  (TH1F*)f->Get("signal"+category[j]+"2");

            //cout<<hdata0->GetEntries()<<"\t"<<hdata1->GetEntries()<<"\t"<<hdata2->GetEntries()<<"\t"<<hdata3->GetEntries()<<endl; 
            //cout<<hmcDs->GetEntries()<<"\t"<<hmcB0->GetEntries()<<"\t"<<hmcBp->GetEntries()<<endl; 

            signal1 = TH1_integral(hsgn1, 1.75, 1.80);
            signal2 = TH1_integral(hsgn2, 1.75, 1.80);
            data1   = TH1_integral(hdata1, 1.75, 1.80);
            data2   = TH1_integral(hdata2, 1.75, 1.80);
            sensitivity1[i] = signal1 / sqrt(data1); 
            sensitivity2[i] = signal2 / sqrt(data2); 
            cout<<fv_nC_cut[i]<<"\t"<<sensitivity1[i]<<"\t"<<data1<<"\t"<<signal1<<"\t"<<sensitivity2[i]<<"\t"<<data2<<"\t"<<signal2<<endl; 
            
            f->Close();
        }
    cout<<"__________________________________________"<<endl;
    }
}
