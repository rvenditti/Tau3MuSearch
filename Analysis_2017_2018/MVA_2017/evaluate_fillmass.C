#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "BDT_optimal_cut.h"

void fill_mass_categories(TChain *t, TH1F* hTriplMass1, TH1F* hTriplMass2, TString category, BDTcut BDTcutvalues, bool isMC){

    //Create TTreeReader
    cout<<"Accessing to input tree"<<endl;
    TTreeReader treeReader(t);

    //create the TMVA::Reader
    TMVA::Tools::Instance();
    TMVA::Reader *MVAreader = new TMVA::Reader("!Color:!Silent:!V");

    //number of spectators
    size_t n_spec = var_spec_names.size();
    std::vector<Float_t> var_spec;
    for(int j = 0; j<n_spec; j++) var_spec.push_back(0);
    // Spectators declaration
    for(int k = 0; k<n_spec; k++){
        MVAreader->TMVA::Reader::AddSpectator(var_spec_names.at(k), &var_spec.at(k));
    }
    //number of variables used for training
    size_t n_train = var_BDT.size();
    std::vector<Float_t> var_train;
    for(int j = 0; j<n_train; j++) var_train.push_back(0);
    // Variables declaration
    for(int k = 0; k<n_train; k++){
        MVAreader->TMVA::Reader::AddVariable(var_BDT.at(k), &var_train.at(k));
    }

    //Book the MVA method(s)
    TString weight_path = TMVA_inputpath+category+TMVA_weightfilename;
    Bool_t weightfileExists = (gSystem->AccessPathName(weight_path) == kFALSE);
    if (weightfileExists) {
       MVAreader->TMVA::Reader::BookMVA(method, weight_path);
       cout<<"Using weights in "<<weight_path<<endl;
    } else {
       std::cout << "Weightfile " <<weight_path<<" for method " << method << " not found."
                    " Did you run TMVACrossValidation with a specified splitExpr?" << std::endl;
       exit(0);
    }

    //Read branches
    std::vector<TTreeReaderValue<double>> reader_spec;
    for(int k = 0; k<n_spec; k++){
        reader_spec.emplace_back(treeReader, var_spec_names.at(k));
    }
    std::vector<TTreeReaderValue<double>> reader_train;
    for(int k = 0; k<n_train; k++){
        reader_train.emplace_back(treeReader, var_names.at(k));
    }
    TTreeReaderValue<double> reader_triplMassReso(treeReader, "tripletMassReso");
    TTreeReaderValue<double> reader_triplMass(treeReader, "tripletMass");

    //Loop on input Tree
    while (treeReader.Next()) {
        for(int k = 0; k<n_spec; k++){
            var_spec.at(k) = *reader_spec.at(k);
        }
        for(int k = 0; k<n_train; k++){
            var_train.at(k) = *reader_train.at(k);
        }
        //Evaluate method(s) and fill histogram or MiniTree
        auto reso = *reader_triplMassReso;
        if(((category == "A" && reso < 0.007 ) ||
            (category == "B" && reso>= 0.007 && reso <= 0.0105 ) ||
            (category == "C" && reso > 0.0105 )) ){
            auto BDTscore = MVAreader->EvaluateMVA(method);
            auto mass = *reader_triplMass; 
            if( BDTscore >= BDTcutvalues.a ) //category 1
                hTriplMass1->Fill(mass);
            else if ( BDTscore < BDTcutvalues.a && BDTscore >= BDTcutvalues.b) //category 2
                hTriplMass2->Fill(mass);
        }
    }
    delete MVAreader;
}


void evaluate_fillmass() 
{
    TString cat_label[] = {"A", "B", "C"};
    TString run[] = {"B", "C", "D", "E", "F"};
    size_t nrun = sizeof(run)/sizeof(run[0]);
    //open input files
    //data
    TChain *tdata = new TChain("FinalTree");
    for(auto i=0; i<nrun; i++){
        tdata->Add(inputpath_datarun[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC Ds
    TChain *tmc1 = new TChain("FinalTree");
    tmc1->Add(inputpath_Ds); 
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    //MC B0
    TChain *tmc2 = new TChain("FinalTree");
    tmc2->Add(inputpath_B0); 
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    //MC Bp
    TChain *tmc3 = new TChain("FinalTree");
    tmc3->Add(inputpath_Bp); 
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;

    //Loop on categories A, B, C
    for(auto i = 0; i<3; i++){
        TString category = cat_label[i];
        cout<<"Category "<<category<<endl;

        //Make sure TChain points to firts event
        tdata->LoadTree(0);
        tmc1->LoadTree(0);
        tmc2->LoadTree(0);
        tmc3->LoadTree(0);

        //use BDT score distributions to set optimal categorisation
        TString file_name = TMVA_inputpath+category+"/BDTdecision_"+category+".root";
        TFile *f = new TFile(file_name,"READ");
        TH1F *h_test_signal;
        TH1F *h_test_bkg;
        h_test_signal = (TH1F*)f->Get("BDTdecision_signal"+category);
        h_test_bkg = (TH1F*)f->Get("BDTdecision_data_obs"+category);
    
        BDTcut BDTcutvalues = Get_BDT_cut(category, h_test_signal, h_test_bkg, false);
        cout<<"BDT cut set based on S and B distribution in "<<file_name<<endl;
        //cout<<"Category "<< category<<" a cut "<<BDTcutvalues.a<<" b cut "<<BDTcutvalues.b<<endl;

        //Book output histogram or MiniTree
        TString fout_path = "datacardT3Mu_"+TMVA_inputpath+category+".root";
        TFile *fout = new TFile(fout_path,"recreate");
        fout->cd();
        TH1F * hTriplMass_data1 = new TH1F ("data_obs"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_data2 = new TH1F ("data_obs"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Ds1   = new TH1F ("signalDs"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Ds2   = new TH1F ("signalDs"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);
        TH1F * hTriplMass_B01   = new TH1F ("signalB0"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_B02   = new TH1F ("signalB0"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Bp1   = new TH1F ("signalBp"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Bp2   = new TH1F ("signalBp"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);
        
        //Read branches, loop on tree, perform evaluate and fill histo
        fill_mass_categories(tdata, hTriplMass_data1, hTriplMass_data2, category, BDTcutvalues, false); 
        fill_mass_categories(tmc1,  hTriplMass_Ds1,   hTriplMass_Ds2,   category, BDTcutvalues, true); 
        fill_mass_categories(tmc2,  hTriplMass_B01,   hTriplMass_B02,   category, BDTcutvalues, true); 
        fill_mass_categories(tmc3,  hTriplMass_Bp1,   hTriplMass_Bp2,   category, BDTcutvalues, true); 

        //rescale signal
        hTriplMass_Ds1->Scale(wNormDs * Ds_correction * Dplus_correction);
        hTriplMass_Ds2->Scale(wNormDs * Ds_correction * Dplus_correction);
        hTriplMass_B01->Scale(wNormB0 * Ds_correction * Bs_correction * f_correction);
        hTriplMass_B02->Scale(wNormB0 * Ds_correction * Bs_correction * f_correction);
        hTriplMass_Bp1->Scale(wNormBp * Ds_correction * Bs_correction * f_correction);
        hTriplMass_Bp2->Scale(wNormBp * Ds_correction * Bs_correction * f_correction);

        //add signals to common histo
        TH1F * hSignal1 = new TH1F ("signal"+category+"1","Triplet mass "+category+"1",42, 1.600000, 2.020000);
        TH1F * hSignal2 = new TH1F ("signal"+category+"2","Triplet mass "+category+"2",42, 1.600000, 2.020000);        

        hSignal1->Add(hTriplMass_Ds1);
        hSignal1->Add(hTriplMass_B01);
        hSignal1->Add(hTriplMass_Bp1);

        hSignal2->Add(hTriplMass_Ds2);
        hSignal2->Add(hTriplMass_B02);
        hSignal2->Add(hTriplMass_Bp2);

        //Write and close the file
        fout->cd();
        hTriplMass_data1->Write();
        hTriplMass_data2->Write();
        hTriplMass_data1->Write("background"+category+"1");
        hTriplMass_data2->Write("background"+category+"2");
        hTriplMass_Ds1->Write();
        hTriplMass_Ds2->Write();
        hTriplMass_B01->Write();
        hTriplMass_B02->Write();
        hTriplMass_Bp1->Write();
        hTriplMass_Bp2->Write();
        hSignal1->Write();
        hSignal2->Write();

        fout->Close();
        cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    }
    return 0;
}
