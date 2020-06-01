#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "../MuonID_study/MuonID_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void fill_MuonID(TChain *t, TH1F* hMuonIDeval, bool isMC){

    //Create TTreeReader
    cout<<"Accessing to input tree"<<endl;
    TTreeReader treeReader(t);

    //create the TMVA::Reader
    TMVA::Tools::Instance();
    TMVA::Reader *MVAreader_barrel = new TMVA::Reader("!Color:!Silent:!V");
    TMVA::Reader *MVAreader_endcap = new TMVA::Reader("!Color:!Silent:!V");

    //number of spectators
    size_t n_spec = var_MuonID_spec_names.size();
    std::vector<Float_t> var_spec;
    for(int j = 0; j<n_spec; j++) var_spec.push_back(0);
    // Spectators declaration
    for(int k = 0; k<n_spec; k++){
        MVAreader_barrel->TMVA::Reader::AddSpectator(var_MuonID_spec_names.at(k), &var_spec.at(k));
        MVAreader_endcap->TMVA::Reader::AddSpectator(var_MuonID_spec_names.at(k), &var_spec.at(k));
        cout<<var_MuonID_spec_names.at(k)<<endl;
    }
    //number of variables used for training
    size_t n_train = var_MuonID_train_def.size();
    std::vector<Float_t> var_train;
    for(int j = 0; j<n_train; j++) var_train.push_back(0);
    // Variables declaration
    for(int k = 0; k<n_train; k++){
        MVAreader_barrel->TMVA::Reader::AddVariable(var_MuonID_train_def.at(k), &var_train.at(k));
        MVAreader_endcap->TMVA::Reader::AddVariable(var_MuonID_train_def.at(k), &var_train.at(k));
        cout<<var_MuonID_train_def.at(k)<<endl;
    }
    TString category = "barrel";
    //Book the MVA method used for MuonID
    TString weight_path_barrel = "/lustrehome/fsimone/MuonID_study/"+TMVA_MuonID_inputpath+category+TMVA_MuonID_weightfilename;
    Bool_t weightfileExists = (gSystem->AccessPathName(weight_path_barrel) == kFALSE);
    if (weightfileExists) {
       MVAreader_barrel->TMVA::Reader::BookMVA(method, weight_path_barrel);
       cout<<"Using weights in "<<weight_path_barrel<<endl;
    } else {
       std::cout << "Weightfile " <<weight_path_barrel<<" for method " << method << " not found."
                    " Did you run TMVACrossValidation with a specified splitExpr?" << std::endl;
       exit(0);
    }
    category = "endcap";
    //Book the MVA method used for MuonID
    TString weight_path_endcap = "/lustrehome/fsimone/MuonID_study/"+TMVA_MuonID_inputpath+category+TMVA_MuonID_weightfilename;
    weightfileExists = (gSystem->AccessPathName(weight_path_endcap) == kFALSE);
    if (weightfileExists) {
       MVAreader_endcap->TMVA::Reader::BookMVA(method, weight_path_endcap);
       cout<<"Using weights in "<<weight_path_endcap<<endl;
    } else {
       std::cout << "Weightfile " <<weight_path_endcap<<" for method " << method << " not found."
                    " Did you run TMVACrossValidation with a specified splitExpr?" << std::endl;
       exit(0);
    }

    //Read branches
    std::vector<TTreeReaderValue<double>> reader_spec;
    for(int k = 0; k<n_spec; k++){
        if(var_MuonID_spec_names.at(k)=="abs(mu_eta)")       {reader_spec.emplace_back(treeReader, "mu_eta"); continue;}
        if(var_MuonID_spec_names.at(k)=="mu_simPdgId")       {reader_spec.emplace_back(treeReader, "mu_eta"); continue;}
        if(var_MuonID_spec_names.at(k)=="mu_simMotherPdgId") {reader_spec.emplace_back(treeReader, "mu_eta"); continue;}
        if(var_MuonID_spec_names.at(k)=="mu_simType")        {reader_spec.emplace_back(treeReader, "mu_eta"); continue;}
        reader_spec.emplace_back(treeReader, var_MuonID_spec_names.at(k));
    }
    std::vector<TTreeReaderValue<double>> reader_train;
    for(int k = 0; k<n_train; k++){
        reader_train.emplace_back(treeReader, var_MuonID_train_names.at(k));
    }
    TTreeReaderValue<double> reader_mueta(treeReader, "mu_eta");

    //Loop on input Tree
    while (treeReader.Next()) {
        for(int k = 0; k<n_spec; k++){
            var_spec.at(k)= *reader_spec.at(k);
        }
        for(int k = 0; k<n_train; k++){
            var_train.at(k) = *reader_train.at(k);
        }
        //Evaluate method(s) and fill histogram or MiniTree
        double BDTscore;
        if(abs(*reader_mueta) < 1.2) BDTscore = MVAreader_barrel->EvaluateMVA(method);
        else BDTscore = MVAreader_endcap->EvaluateMVA(method);
        hMuonIDeval->Fill(BDTscore);
    }
    delete MVAreader_barrel;
    delete MVAreader_endcap;
}


void evaluate_fillMuonID(TString Mu) 
{
    if( !( Mu=="Mu1" || Mu=="Mu2" || Mu=="Mu3" ) ){cout<<"Wrong input arg"<<endl; return;}
    TString treeName = "Tree"+Mu;
    //open input files
    //data
    size_t n_run = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    TChain *tdata = new TChain(treeName);
    for(auto i=0; i<n_run; i++){
        tdata->Add(inputpath_datarun[i]);
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }

    //MC Ds
    TChain *tmc1 = new TChain(treeName);
    tmc1->Add(inputpath_Ds);
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    //MC B0
    TChain *tmc2 = new TChain(treeName);
    tmc2->Add(inputpath_B0); 
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    //MC Bp
    TChain *tmc3 = new TChain("FinalTree");
    tmc3->Add(inputpath_Bp); 
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;


        //Make sure TChain points to firts event
        tdata->LoadTree(0);
        tmc1->LoadTree(0);
        tmc2->LoadTree(0);
        tmc3->LoadTree(0);
        TString category = "A";
        //Book output histogram or MiniTree
        TString fout_path = TMVA_inputpath+category+"/"+Mu+"_MuonIDeval.root";
        TFile *fout = new TFile(fout_path, "recreate");
        fout->cd();

        TH1F * hMuonIDeval_data = new TH1F ("MuonIDeval_data","MuonIDeval_data", 400, -1.0, 1.0);
        TH1F * hMuonIDeval_Ds = new TH1F ("MuonIDeval_signalDs","MuonIDeval_signalDs", 400, -1.0, 1.0);
        TH1F * hMuonIDeval_B0 = new TH1F ("MuonIDeval_signalB0","MuonIDeval_signalB0", 400, -1.0, 1.0);
        TH1F * hMuonIDeval_Bp = new TH1F ("MuonIDeval_signalBp","MuonIDeval_signalBp", 400, -1.0, 1.0);

        //Read branches, loop on tree, perform evaluate and fill histo
        fill_MuonID(tdata, hMuonIDeval_data, false); 
        fill_MuonID(tmc1,  hMuonIDeval_Ds, true); 
        fill_MuonID(tmc2,  hMuonIDeval_B0, true); 
        fill_MuonID(tmc3,  hMuonIDeval_Bp, true); 

        hMuonIDeval_data->Scale(1/hMuonIDeval_data->Integral());
        //rescale signal
        hMuonIDeval_Ds->Scale(wNormDs * Ds_correction * Dplus_correction);
        hMuonIDeval_B0->Scale(wNormB0 * Ds_correction * Bs_correction * f_correction);
        hMuonIDeval_Bp->Scale(wNormBp * Ds_correction * Bs_correction * f_correction);

        //add signals to common histo
        TH1F * hMuonIDevalSignal = new TH1F ("MuonIDeval_signal", "MuonIDeval_signal", 400, -1.0, 1.0);
        hMuonIDevalSignal->Add(hMuonIDeval_Ds);
        hMuonIDevalSignal->Add(hMuonIDeval_B0);
        hMuonIDevalSignal->Add(hMuonIDeval_Bp);
        hMuonIDevalSignal->Scale(1/hMuonIDevalSignal->Integral());

        //Write and close the file
        fout->cd();
        hMuonIDeval_data->Write();
        hMuonIDevalSignal->Write();
        hMuonIDeval_Ds->Write();
        hMuonIDeval_B0->Write();
        hMuonIDeval_Bp->Write();

        fout->Close();
        cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;

    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}
