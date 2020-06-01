#include "TH1F.h"
#include "TFile.h"
#include <cmath>
#include <string> 
#include "T3M_common.h"
#include "../MuonID_study/MuonID_common.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void filltree_MuonID(TTree *tin, TTree *tout, Double_t &value_MuonID){

    //Create TTreeReader
    cout<<"Accessing to input tree"<<endl;
    TTreeReader treeReader(tin);

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
        reader_spec.emplace_back(treeReader, var_MuonID_spec_names.at(k));
    }
    std::vector<TTreeReaderValue<double>> reader_train;
    for(int k = 0; k<n_train; k++){
        reader_train.emplace_back(treeReader, var_MuonID_train_names.at(k));
    }
    TTreeReaderValue<double> reader_mueta(treeReader, "mu_eta");

    //prepare branch in output tree
    //Loop on input Tree
    Double_t score = 0;
    while (treeReader.Next()) {
        for(int k = 0; k<n_spec; k++){
            var_spec.at(k)= *reader_spec.at(k);
        }
        for(int k = 0; k<n_train; k++){
            var_train.at(k) = *reader_train.at(k);
        }
        
        //Evaluate method(s) and fill histogram or MiniTree
        if(abs(*reader_mueta) < 1.2) score = MVAreader_barrel->EvaluateMVA(method);
        else score  = MVAreader_endcap->EvaluateMVA(method);
        value_MuonID = score;
        tout->Fill();
    }
    delete MVAreader_barrel;
    delete MVAreader_endcap;
}


void evaluate_fillMuonID_tree(TString inputFile) 
{
    //open input file
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(inputFile);
    if (!f || !f->IsOpen()) f = new TFile(inputFile);
    std::cout<<"Opened input file: "<<inputFile<<std::endl;
    //Access Tree
    TTree *tin1 = (TTree*)f->Get("TreeMu1");
    TTree *tin2 = (TTree*)f->Get("TreeMu2");
    TTree *tin3 = (TTree*)f->Get("TreeMu3");

    Double_t value_MuonID = 0;
    //Book output Tree
    TString fout_path = inputFile.ReplaceAll(".root", "_MuonID.root");
    cout<<"Output file: "<<fout_path<<endl;
    TFile *fout = new TFile(fout_path, "recreate");
    fout->cd();
    TTree *tout1 = new TTree("MuonIDeval_Mu1","MuonIDeval_Mu1");
    tout1->Branch("MuonID",&value_MuonID,"MuonID/D");
    TTree *tout2 = new TTree("MuonIDeval_Mu2","MuonIDeval_Mu2");
    tout2->Branch("MuonID",&value_MuonID,"MuonID/D");
    TTree *tout3 = new TTree("MuonIDeval_Mu3","MuonIDeval_Mu3");
    tout3->Branch("MuonID",&value_MuonID,"MuonID/D");

    //Read branches, loop on tree, perform evaluate and fill histo
    filltree_MuonID(tin1, tout1, value_MuonID); 
    filltree_MuonID(tin2, tout2, value_MuonID); 
    filltree_MuonID(tin3, tout3, value_MuonID); 

    fout->Write();
    fout->Close();
    cout<<"+++++++++++++++++++++++++++\nWritten output file: "<<fout_path<<"\n\n"<<endl;

    return 0;
}

void run_fillMuonID_tree(){
    int n_data = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    for(int i=0; i<n_data; i++)
        evaluate_fillMuonID_tree(inputpath_datarun[i]);

    evaluate_fillMuonID_tree(inputpath_Ds);
    evaluate_fillMuonID_tree(inputpath_B0);
    evaluate_fillMuonID_tree(inputpath_Bp);

    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
}
