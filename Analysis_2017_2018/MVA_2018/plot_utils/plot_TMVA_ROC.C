//root -l -b plot_TMVA_bdtcorr.C\(\"C\"\)

#include "TH1F.h"
#include <cmath>
#include <string> 
#include "/lustrehome/fsimone/MVA_2018/T3M_common.h"
#include "/lustrehome/fsimone/MVA_2018/BDT_optimal_cut.h"

void plot_TMVA_ROC(TString category) 
{
    bool compute_optBDT = false;
    TString TMVA_filename_list[25];
    TMVA_filename_list[0] = "dataset_2018_5may_baseline_";
    for(int i = 1; i<10; i++) {
        auto s = std::to_string(i);
        TMVA_filename_list[i] = "dataset_2018_5may_v"+s+"_";
    }
    for(int i = 10; i<25; i++) {
        auto s = std::to_string(i);
        TMVA_filename_list[i] = "dataset_2018_6may_v"+s+"_";
    }

    TCanvas *c1 = new TCanvas("c1","ROC "+category,150,10,990,660);
    TLegend*leg = new TLegend(0.1,0.7,0.1,0.4);
    TH1F *ROC[25];
    THStack hs("ROC "+category,"ROC "+category);
    for(int j = 0; j<25; j++){
        TString TMVA_filename = TMVA_filename_list[j]+category;

        TFile *f_tmva = new TFile("/lustrehome/fsimone/MVA_2018/TMVA_"+TMVA_filename+".root","READ");
        //cout<<"Opened TMVA file "<<"TMVA_"+TMVA_filename+".root"<<endl;

        ROC[j] = (TH1F*)f_tmva->Get(TMVA_filename+"/Method_BDT/BDT/MVA_BDT_rejBvsS");
        ROC[j]->GetYaxis()->SetRangeUser(0.5,1.05);
        ROC[j]->SetLineColor(1+j);
        leg->AddEntry(ROC[j],TMVA_filename,"l");
        hs.Add(ROC[j]);
        ROC[j]->Draw("same");
        cout<<TMVA_filename<<" ROC[0,1] "<<TH1_integral(ROC[j], 0., 1.)<<" ROC[0,.9] "<<TH1_integral(ROC[j], 0., 0.9)<< " ROC[0,0.8] "<<TH1_integral(ROC[j], 0., 0.8)<< " ROC[0,0.7] "<<TH1_integral(ROC[j], 0., 0.7)<<" ROC[0,0.6] "<<TH1_integral(ROC[j], 0., 0.6)<<" ROC[0,0.5] "<<TH1_integral(ROC[j], 0., 0.5)<<" ROC[0,0.4] "<<TH1_integral(ROC[j], 0., 0.4)<<" ROC[0,0.3] "<<TH1_integral(ROC[j], 0., 0.3)<< endl;

        if(compute_optBDT) {
            //use BDT score distributions to set optimal categorisation
            TString file_name = "/lustrehome/fsimone/MVA_2018/"+TMVA_inputpath+category+"/BDTdecision_"+category+".root";
            TFile *f = new TFile(file_name,"READ");
            TH1F *h_test_signal;
            TH1F *h_test_bkg;
            h_test_signal = (TH1F*)f->Get("BDTdecision_signal"+category);
            h_test_bkg = (TH1F*)f->Get("BDTdecision_data_obs"+category);

            BDTcut BDTcutvalues = Get_BDT_cut(category, h_test_signal, h_test_bkg, false);
            cout<<"BDT cut set based on S and B distribution in "<<file_name<<endl;

            TTree *tTrain = (TTree*)f_tmva->Get(TMVA_filename+"/TrainTree");
            TTreeReader reader (tTrain);

            TTreeReaderValue<float> reader_bdt (reader, "BDT");
            TTreeReaderValue<int> reader_classid (reader, "classID");

            int nentries = tTrain->GetEntries();
            int signal_class = 0;

            double n_b_true_positive = 0;
            double n_b_false_positive = 0;
            double n_b_true_negative_rej = 0;
            double n_b_false_negative_rej = 0;
            double n_a_true_positive = 0;
            double n_a_false_positive = 0;
            double n_a_true_negative_rej = 0;
            double n_a_false_negative_rej = 0;
            
            while (reader.Next()) {
                if ((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 0) {
                    n_b_true_positive += 1.;
                } else if((*reader_bdt) >= BDTcutvalues.b and (*reader_classid) == 1){
                    n_b_false_positive += 1.;
                }

                if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 1) {
                    n_b_true_negative_rej += 1.;
                } else if ((*reader_bdt) < BDTcutvalues.b and (*reader_classid) == 0) {
                    n_b_false_negative_rej += 1.;
                }

                if ((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 0) {
                    n_a_true_positive += 1.;
                } else if((*reader_bdt) >= BDTcutvalues.a and (*reader_classid) == 1){
                    n_a_false_positive += 1.;
                }

                if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 1) {
                    n_a_true_negative_rej += 1.;
                } else if ((*reader_bdt) < BDTcutvalues.a and (*reader_classid) == 0) {
                    n_a_false_negative_rej += 1.;
                }
             }
            cout<<"BDT working point b "<<BDTcutvalues.b<<endl;
            cout<<"n_b_true_positive "<<n_b_true_positive<<" n_b_false_positive "<<n_b_false_positive<<endl;
            cout<<"n_b_true_negative_rej "<<n_b_true_negative_rej<<" n_b_false_negative_rej "<<n_b_false_negative_rej<<endl;
            cout<<"signal efficiency = TP/P "<<n_b_true_positive/(n_b_true_positive+n_b_false_negative_rej)<<endl;
            cout<<"bkg rejection = TN_rej/N "<<n_b_true_negative_rej/(n_b_false_positive+n_b_true_negative_rej)<<endl;
            cout<<"-------\nBDT working point a "<<BDTcutvalues.a<<endl;
            cout<<"n_a_true_positive "<<n_a_true_positive<<" n_a_false_positive "<<n_a_false_positive<<endl;
            cout<<"n_a_true_negative_rej "<<n_a_true_negative_rej<<" n_a_false_negative_rej "<<n_a_false_negative_rej<<endl;
            cout<<"signal efficiency = TP/P "<<n_a_true_positive/(n_a_true_positive+n_a_false_negative_rej)<<endl;
            cout<<"bkg rejection = TN_rej/N "<<n_a_true_negative_rej/(n_a_false_positive+n_a_true_negative_rej)<<endl;
            Double_t x[2], y[2];
            x[0] = n_b_true_positive/(n_b_true_positive+n_b_false_negative_rej);
            y[0] = n_b_true_negative_rej/(n_b_false_positive+n_b_true_negative_rej);
            x[1] = n_a_true_positive/(n_a_true_positive+n_a_false_negative_rej);
            y[1] = n_a_true_negative_rej/(n_a_false_positive+n_a_true_negative_rej);
            TGraph* gr = new TGraph(2,x,y);
            TLatex *l = new TLatex(0.5, 0.5, "label");
            l->SetTextSize(0.025);
            l->SetTextFont(42);
            l->SetTextAlign(21);
            l->SetTextColor(kBlue);
            l->DrawLatex(x[0],y[0]+0.01,Form("bdt.b %4.2f",BDTcutvalues.b));
            l->DrawLatex(x[1],y[1]+0.01,Form("bdt.a %4.2f",BDTcutvalues.a));
            gr->Draw("same *p");
        }
    }
    leg->Draw();
    //hs.Draw("hist nostack");
    //hs.GetXaxis()->SetTitle("Signal efficiency");
    //hs.GetYaxis()->SetTitle("Background rejection");
    //leg->Draw();
    gPad->Modified();
    c1->Update();
}
