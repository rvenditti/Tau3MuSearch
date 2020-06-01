#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "../T3M_common.h"
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

double log_significance(double S, double B){
    double significance = 0;
    significance = sqrt(2*( (S+B) * log( 1+S/B ) - S));
    //cout<<"log sign is "<<significance<<" while S/sqrt(S + B) gives "<< S/sqrt(S + B) <<endl;
    return significance;
}

void count_datacard() 
{
    int ncat = sizeof(cat_label)/sizeof(cat_label[0]);
    //TString cat_label[] = {"A", "B", "C"};
    TString additional_cut[] = {
        //             "v2",
                     "ptsplit",
                   //  "fv_nC<10000",
                   //  "fv_nC<5000",
                   //  "fv_nC<2500",
                   //  "fv_nC<1000",
                   //  "fv_nC<500",
                   //  "fv_nC<100",
                   //  "fv_nC<50",
                   //  "fv_nC<40",
                   //  "fv_nC<30",
                   //  "fv_nC<20",
                   //  "fv_nC<15",
                   //  "fv_nC<10",
                   //  "fv_nC<8",
                   //  "fv_nC<6",
                   //  "fv_nC<4",
    };
    size_t n = sizeof(additional_cut)/sizeof(additional_cut[0]);
    float sensitivity1[n];
    float sensitivity2[n];
    float signal1, signal2, data1, data2;

    //count reparately for each cat_label
    for(auto j = 0; j<ncat; j++){
        cout<<"category "<<cat_label[j]<<endl;
        cout<<"additional_cut\tsignif1\tbkg1\tsignal1\tsignif2\tbkg2\tsignal2"<<endl; 
        TString cat_histo = cat_label[j]; 
        cat_histo =cat_histo.ReplaceAll("lowpt",""); cat_histo = cat_histo.ReplaceAll("highpt","");

        //range for counting events
        double x_low = 1.60; // 1.75; //1.6;
        double x_high = 2.02; //1.80; //2.02;
        //for each additional_cut i open datacard and compute the sensitivity after BDT
        for(auto i = 0; i<n; i++){

            //TString filename = "datacardT3Mu_dataset_2018_27april_Chi2_"+additional_cut[i]+"_"+cat_label[j]+".root";
            TString filename = "/lustrehome/fsimone/MVA_2018/datacardT3Mu_dataset_2018_6may_optimised_"+additional_cut[i]+"_noptlabel.root";
            //TString filename = "/lustrehome/fsimone/MVA_2018/datacardT3Mu_dataset_2018_6may_optimised_"+additional_cut[i]+"_"+cat_label[j]+".root";
            TFile *f = TFile::Open(filename, "READ");
            if (!f || !f->IsOpen()) { std::cout<<"Error"<<endl; continue;}

            TH1F *hdata1 = (TH1F*)f->Get("data_obs"+cat_histo+"1");
            TH1F *hdata2 = (TH1F*)f->Get("data_obs"+cat_histo+"2");
            TH1F *hsgn1 =  (TH1F*)f->Get("signal"+cat_histo+"1");
            TH1F *hsgn2 =  (TH1F*)f->Get("signal"+cat_histo+"2");

            //cout<<hdata0->GetEntries()<<"\t"<<hdata1->GetEntries()<<"\t"<<hdata2->GetEntries()<<"\t"<<hdata3->GetEntries()<<endl; 
            //cout<<hmcDs->GetEntries()<<"\t"<<hmcB0->GetEntries()<<"\t"<<hmcBp->GetEntries()<<endl; 

            signal1 = TH1_integral(hsgn1,  x_low, x_high);
            signal2 = TH1_integral(hsgn2,  x_low, x_high);
            data1   = TH1_integral(hdata1, x_low, x_high);
            data2   = TH1_integral(hdata2, x_low, x_high);
            sensitivity1[i] = signal1 / sqrt(data1); 
            sensitivity2[i] = signal2 / sqrt(data2); 
            double signif1 = log_significance(signal1, data1);
            double signif2 = log_significance(signal2, data2);
            cout<<additional_cut[i]<<"\t"<<signif1<<"\t"<<data1<<"\t"<<signal1<<"\t"<<signif2<<"\t"<<data2<<"\t"<<signal2<<endl; 
            
            f->Close();
        }
    cout<<"__________________________________________"<<endl;
    }
}
