#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void fill_BDT_score(TChain *t, TH1F* hBDTdecision, TString category, bool isMC){

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
            var_spec.at(k)= *reader_spec.at(k);
        }
        for(int k = 0; k<n_train; k++){
            var_train.at(k) = *reader_train.at(k);
        }
        //Evaluate method(s) and fill histogram or MiniTree
        auto BDTscore = MVAreader->EvaluateMVA(method);
        auto reso = *reader_triplMassReso;
        auto mass = *reader_triplMass;
        bool isSB = true;
        if(!isMC){
            if( ( mass >= 1.62 && mass <= 1.75 ) || ( mass >= 1.80 && mass <= 2.0 ) ) isSB = true; 
            else isSB = false;
        }
        if(((category == "A" && reso < 0.007 ) ||
            (category == "B" && reso>= 0.007 && reso <= 0.0105 ) ||
            (category == "C" && reso > 0.0105 )) && isSB ){
            hBDTdecision->Fill(BDTscore);
           // for(int k = 0; k<n_spec; k++){
           //     cout<<"   spect. "<<var_spec.at(k)<<endl;
           // }
           // for(int k = 0; k<n_train; k++){
           //     cout<<"   train "<<var_train.at(k)<<endl;
           // }
           // cout<<"reso "<<reso<<" bdt "<<MVAreader->EvaluateMVA(method)<<endl;        
        }
    }
    delete MVAreader;
}


void evaluate_fillbdt() 
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
  
        //Book output histogram or MiniTree
        TString fout_path = TMVA_inputpath+category+"/BDTdecision_"+category+".root";
        TFile *fout = new TFile(fout_path, "recreate");
        fout->cd();

        TH1F * hBDTdecision_data = new TH1F ("BDTdecision_data_obs"+category,"BDTdecision_data_obs"+category, 240, -0.6, 0.6);
        TH1F * hBDTdecision_Ds = new TH1F ("BDTdecision_signalDs"+category,"BDTdecision_signalDs"+category, 240, -0.6, 0.6);
        TH1F * hBDTdecision_B0 = new TH1F ("BDTdecision_signalB0"+category,"BDTdecision_signalB0"+category, 240, -0.6, 0.6);
        TH1F * hBDTdecision_Bp = new TH1F ("BDTdecision_signalBp"+category,"BDTdecision_signalBp"+category, 240, -0.6, 0.6);

        //Read branches, loop on tree, perform evaluate and fill histo
        fill_BDT_score(tdata, hBDTdecision_data, category, false); 
        fill_BDT_score(tmc1,  hBDTdecision_Ds,   category, true); 
        fill_BDT_score(tmc2,  hBDTdecision_B0,   category, true); 
        fill_BDT_score(tmc3,  hBDTdecision_Bp,   category, true); 

        //rescale signal
        hBDTdecision_Ds->Scale(wNormDs * Ds_correction * Dplus_correction);
        hBDTdecision_B0->Scale(wNormB0 * Ds_correction * Bs_correction * f_correction);
        hBDTdecision_Bp->Scale(wNormBp * Ds_correction * Bs_correction * f_correction);

        //add signals to common histo
        TH1F * hBDTdecisionSignal = new TH1F ("BDTdecision_signal"+category, "BDTdecision_signal"+category, 240, -0.6, 0.6);
        hBDTdecisionSignal->Add(hBDTdecision_Ds);
        hBDTdecisionSignal->Add(hBDTdecision_B0);
        hBDTdecisionSignal->Add(hBDTdecision_Bp);

        //Write and close the file
        fout->cd();
        hBDTdecision_data->Write();
        hBDTdecisionSignal->Write();
        hBDTdecision_Ds->Write();
        hBDTdecision_B0->Write();
        hBDTdecision_Bp->Write();

        fout->Close();
        cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;
    }
    return 0;
}
