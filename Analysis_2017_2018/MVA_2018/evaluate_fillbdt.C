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

    //Variables
    std::vector<TString> var_train_name;
    std::vector<TString> var_train_def;
    std::vector<TString> var_spec_name;
    std::vector<TString> var_spec_def;
    TString BDTinVar;
    if(category.Contains("A")) BDTinVar = BDTinVar_A;
    if(category.Contains("B")) BDTinVar = BDTinVar_B;
    if(category.Contains("C")) BDTinVar = BDTinVar_C;

    readVarName(var_spec_name, var_spec_def, BDTspecVar);
    readVarName(var_train_name, var_train_def, BDTinVar);

    // Spectators declaration
    cout<<"Declaration of spectator variables - category "<<category<<" from file:"<<BDTspecVar<<endl;
    for(int k = 0; k<var_spec_name.size(); k++){
        cout<<k<<" - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
    }

    // Variables declaration
    cout<<"Declaration of variables for training - category "<<category<<" from file:"<<BDTinVar<<endl;
    for(int k = 0; k<var_train_name.size(); k++){
        cout<<k<<" - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
    }

    //number of spectators
    size_t n_spec = var_spec_name.size();
    std::vector<Float_t> var_spec;
    for(int j = 0; j<n_spec; j++) var_spec.push_back(0);
    // Spectators declaration
    for(int k = 0; k<n_spec; k++){
        MVAreader->TMVA::Reader::AddSpectator(var_spec_def.at(k), &var_spec.at(k));
        cout<<k<<" added to Reader - "<<var_spec_name.at(k)<<" - "<<var_spec_def.at(k)<<endl;
    }
    //number of variables used for training
    size_t n_train = var_train_name.size();
    std::vector<Float_t> var_train;
    for(int j = 0; j<n_train; j++) var_train.push_back(0);
    // Variables declaration
    for(int k = 0; k<n_train; k++){
        MVAreader->TMVA::Reader::AddVariable(var_train_def.at(k), &var_train.at(k));
        cout<<k<<" added to Reader - "<<var_train_name.at(k)<<" - "<<var_train_def.at(k)<<endl;
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
    //Branches for preselections
    TTreeReaderValue<double> reader_mu1_tLWM(treeReader, "TreeMu1.mu_trackerLayersWithMeasurement");
    TTreeReaderValue<double> reader_mu2_tLWM(treeReader, "TreeMu2.mu_trackerLayersWithMeasurement");
    TTreeReaderValue<double> reader_mu3_tLWM(treeReader, "TreeMu3.mu_trackerLayersWithMeasurement");
    TTreeReaderValue<double> reader_mu1_validFra(treeReader, "TreeMu1.mu_innerTrack_validFraction");
    TTreeReaderValue<double> reader_mu2_validFra(treeReader, "TreeMu2.mu_innerTrack_validFraction");
    TTreeReaderValue<double> reader_mu3_validFra(treeReader, "TreeMu3.mu_innerTrack_validFraction");
    TTreeReaderValue<double> reader_mu1_iTnC(treeReader, "TreeMu1.mu_innerTrack_normalizedChi2");
    TTreeReaderValue<double> reader_mu2_iTnC(treeReader, "TreeMu2.mu_innerTrack_normalizedChi2");
    TTreeReaderValue<double> reader_mu3_iTnC(treeReader, "TreeMu3.mu_innerTrack_normalizedChi2");
    TTreeReaderValue<double> reader_mu1_tKink(treeReader, "TreeMu1.mu_combinedQuality_trkKink");
    TTreeReaderValue<double> reader_mu2_tKink(treeReader, "TreeMu2.mu_combinedQuality_trkKink");
    TTreeReaderValue<double> reader_mu3_tKink(treeReader, "TreeMu3.mu_combinedQuality_trkKink");
    TTreeReaderValue<double> reader_fv_nC(treeReader, "fv_nC");

    TTreeReaderValue<double> reader_dxy3(treeReader, "dxy3");
    TTreeReaderValue<double> reader_dxyErr3(treeReader, "dxyErr3");
    TTreeReaderValue<double> reader_SVx(treeReader, "SVx");
    TTreeReaderValue<double> reader_SVy(treeReader, "SVy");
    TTreeReaderValue<double> reader_SVz(treeReader, "SVz");
    TTreeReaderValue<double> reader_PVx(treeReader, "RefVx1");
    TTreeReaderValue<double> reader_PVy(treeReader, "RefVy1");
    TTreeReaderValue<double> reader_PVz(treeReader, "RefVz1");


    //Read branches
    std::vector<TTreeReaderValue<double>> reader_spec;
    for(int k = 0; k<n_spec; k++){
        cout<<"   reading from TChain branch "<<var_spec_name.at(k)<<endl;
        reader_spec.emplace_back(treeReader, var_spec_name.at(k));
    }
    std::vector<TTreeReaderValue<double>> reader_train;
    for(int k = 0; k<n_train; k++){
        if(var_train_name.at(k)=="abs(dxy3/dxyErr3)" || var_train_def.at(k)=="abs(dxy3/dxyErr3)") continue;
        if(var_train_name.at(k)=="PS_SV_dxy" || var_train_def.at(k)=="sqrt((RefVx1-SVx)*(RefVx1-SVx)+(RefVy1-SVy)*(RefVy1-SVy))") continue;
        if(var_train_name.at(k)=="PS_SV_dz" || var_train_def.at(k)=="abs(RefVz1-SVz)" || var_train_def.at(k)=="abs(SVz-RefVz1)") continue;
        cout<<"   reading from TChain branch "<<var_train_name.at(k)<<endl;
        reader_train.emplace_back(treeReader, var_train_name.at(k));
    }
    TTreeReaderValue<double> reader_triplMassReso(treeReader, "tripletMassReso");
    TTreeReaderValue<double> reader_Ptmu1(treeReader, "Ptmu1");
    TTreeReaderValue<double> reader_triplMass(treeReader, "tripletMass");
    int count = 0;
    //Loop on input Tree
    while (treeReader.Next()) {
        count++;
        bool passedPresel = false;
        if (( *reader_fv_nC>0 && *reader_fv_nC<100) &&
            ( *reader_mu1_tLWM>8 ) &&
            ( *reader_mu2_tLWM>8 ) &&
            ( *reader_mu3_tLWM>8 ) &&
            ( *reader_mu1_validFra>0.5 ) &&
            ( *reader_mu2_validFra>0.5 ) &&
            ( *reader_mu3_validFra>0.5 ) &&
            ( *reader_mu1_iTnC<40 )&&
            ( *reader_mu2_iTnC<40 )&&
            ( *reader_mu3_iTnC<40 )&&
            ( *reader_mu1_tKink<900 ) &&
            ( *reader_mu2_tKink<900 ) &&
            ( *reader_mu3_tKink<900 ) ) passedPresel = true;
        if(!passedPresel) continue;

        for(int k = 0; k<n_spec; k++){
            var_spec.at(k) = *reader_spec.at(k);
        }
        for(int k = 0; k<n_train; k++){
            if(var_train_name.at(k)=="abs(dxy3/dxyErr3)" || var_train_def.at(k)=="abs(dxy3/dxyErr3)"){ 
                if((*reader_dxyErr3)!=0) var_train.at(k) = abs( (*reader_dxy3) / (*reader_dxyErr3) );
                else var_train.at(k) = 0;
            }
            else if(var_train_name.at(k)=="PS_SV_dxy" || var_train_def.at(k)=="sqrt((RefVx1-SVx)*(RefVx1-SVx)+(RefVy1-SVy)*(RefVy1-SVy))")
                var_train.at(k) = sqrt(pow(( (*reader_PVx) - (*reader_SVx) ),2) + pow(( (*reader_PVy) - (*reader_SVy) ),2));
            else if(var_train_name.at(k)=="PS_SV_dz" || var_train_def.at(k)=="abs(RefVz1-SVz)" || var_train_def.at(k)=="abs(SVz-RefVz1)")
                var_train.at(k) = abs( (*reader_PVz) - (*reader_SVz) );
            else var_train.at(k) = *reader_train.at(k);
        }
        //Evaluate method(s) and fill histogram or MiniTree
        auto BDTscore = MVAreader->EvaluateMVA(method);
        auto reso = *reader_triplMassReso;
        auto mu1pt = *reader_Ptmu1;
        auto mass = *reader_triplMass;
        bool isSB = true;
        if(!isMC){
            if( ( mass >= 1.62 && mass <= 1.75 ) || ( mass >= 1.80 && mass <= 2.0 ) ) isSB = true; 
            else isSB = false;
        }
        if(((category.Contains("A") && reso < 0.007 ) ||
            (category.Contains("B") && reso>= 0.007 && reso <= 0.0105 ) ||
            (category.Contains("C") && reso > 0.0105 )) && isSB ){
            if(((category.Contains("high") && mu1pt >= 7 ) ||
                (category.Contains("low")  && mu1pt <  7 ) ||
                (!(category.Contains("low")) && !(category.Contains("high"))))){
                hBDTdecision->Fill(BDTscore);
                //for(int k = 0; k<n_spec; k++){
                //    cout<<k<<"   spect. "<<var_spec_def.at(k)<<" "<<var_spec.at(k)<<endl;
                //}
                //for(int k = 0; k<n_train; k++){
                //    cout<<k<<"   train "<<var_train_def.at(k)<<" "<<var_train.at(k)<<endl;
                //}
                //cout<<"reso "<<reso<<" bdt "<<MVAreader->EvaluateMVA(method)<<endl;        
            }
        }
    }
    delete MVAreader;
}


void evaluate_fillbdt() 
{
    int ncat = sizeof(cat_label)/sizeof(cat_label[0]);

    int nrun = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
    //open input files
    std::vector<TTree*> sigTree, dataTree;

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

    //data
    TChain *tdata_mu1 = new TChain("TreeMu1");
    TChain *tdata_mu2 = new TChain("TreeMu2");
    TChain *tdata_mu3 = new TChain("TreeMu3");
    for(auto i=0; i<nrun; i++){
        tdata_mu1->Add(inputpath_datarun[i]); 
        tdata_mu2->Add(inputpath_datarun[i]); 
        tdata_mu3->Add(inputpath_datarun[i]); 
        std::cout<<"Opened input file: "<<inputpath_datarun[i]<<std::endl;
    }
    //MC Ds
    TChain *tmc1_mu1 = new TChain("TreeMu1");
    TChain *tmc1_mu2 = new TChain("TreeMu2");
    TChain *tmc1_mu3 = new TChain("TreeMu3");
    tmc1_mu1->Add(inputpath_Ds); 
    tmc1_mu2->Add(inputpath_Ds); 
    tmc1_mu3->Add(inputpath_Ds); 
    std::cout<<"Opened input file: "<<inputpath_Ds<<std::endl;
    //MC B0
    TChain *tmc2_mu1 = new TChain("TreeMu1");
    TChain *tmc2_mu2 = new TChain("TreeMu2");
    TChain *tmc2_mu3 = new TChain("TreeMu3");
    tmc2_mu1->Add(inputpath_B0); 
    tmc2_mu2->Add(inputpath_B0); 
    tmc2_mu3->Add(inputpath_B0); 
    std::cout<<"Opened input file: "<<inputpath_B0<<std::endl;
    //MC Bp
    TChain *tmc3_mu1 = new TChain("TreeMu1");
    TChain *tmc3_mu2 = new TChain("TreeMu2");
    TChain *tmc3_mu3 = new TChain("TreeMu3");
    tmc3_mu1->Add(inputpath_Bp); 
    tmc3_mu2->Add(inputpath_Bp); 
    tmc3_mu3->Add(inputpath_Bp); 
    std::cout<<"Opened input file: "<<inputpath_Bp<<std::endl;

    //data
    TChain *tdata_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tdata_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tdata_muid3 = new TChain("MuonIDeval_Mu3");
    for(auto i=0; i<nrun; i++){
        TString inputpath_data_muId = inputpath_datarun[i].ReplaceAll(".root", "_MuonID.root");
        tdata_muid1->Add(inputpath_datarun[i]); 
        tdata_muid2->Add(inputpath_datarun[i]); 
        tdata_muid3->Add(inputpath_datarun[i]); 
        std::cout<<"Opened input file: "<<inputpath_data_muId<<std::endl;
    }
    //MC Ds
    TChain *tmc1_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc1_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc1_muid3 = new TChain("MuonIDeval_Mu3");
    TString inputpath_Ds_muId = inputpath_Ds.ReplaceAll(".root", "_MuonID.root");
    tmc1_muid1->Add(inputpath_Ds_muId); 
    tmc1_muid2->Add(inputpath_Ds_muId); 
    tmc1_muid3->Add(inputpath_Ds_muId); 
    std::cout<<"Opened input file: "<<inputpath_Ds_muId<<std::endl;
    //MC B0
    TChain *tmc2_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc2_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc2_muid3 = new TChain("MuonIDeval_Mu3");
    TString inputpath_B0_muId = inputpath_B0.ReplaceAll(".root", "_MuonID.root");
    tmc2_muid1->Add(inputpath_B0_muId); 
    tmc2_muid2->Add(inputpath_B0_muId); 
    tmc2_muid3->Add(inputpath_B0_muId); 
    std::cout<<"Opened input file: "<<inputpath_B0_muId<<std::endl;
    //MC Bp
    TChain *tmc3_muid1 = new TChain("MuonIDeval_Mu1");
    TChain *tmc3_muid2 = new TChain("MuonIDeval_Mu2");
    TChain *tmc3_muid3 = new TChain("MuonIDeval_Mu3");
    TString inputpath_Bp_muId = inputpath_Bp.ReplaceAll(".root", "_MuonID.root");
    tmc3_muid1->Add(inputpath_Bp_muId); 
    tmc3_muid2->Add(inputpath_Bp_muId); 
    tmc3_muid3->Add(inputpath_Bp_muId); 
    std::cout<<"Opened input file: "<<inputpath_Bp_muId<<std::endl;

    tdata->AddFriend(tdata_mu1);
    tdata->AddFriend(tdata_mu2);
    tdata->AddFriend(tdata_mu3);
    tdata->AddFriend(tdata_muid1);
    tdata->AddFriend(tdata_muid2);
    tdata->AddFriend(tdata_muid3);

    tmc1->AddFriend(tmc1_mu1);
    tmc1->AddFriend(tmc1_mu2);
    tmc1->AddFriend(tmc1_mu3);
    tmc1->AddFriend(tmc1_muid1);
    tmc1->AddFriend(tmc1_muid2);
    tmc1->AddFriend(tmc1_muid3);

    tmc2->AddFriend(tmc2_mu1);
    tmc2->AddFriend(tmc2_mu2);
    tmc2->AddFriend(tmc2_mu3);
    tmc2->AddFriend(tmc2_muid1);
    tmc2->AddFriend(tmc2_muid2);
    tmc2->AddFriend(tmc2_muid3);

    tmc3->AddFriend(tmc3_mu1);
    tmc3->AddFriend(tmc3_mu2);
    tmc3->AddFriend(tmc3_mu3);
    tmc3->AddFriend(tmc3_muid1);
    tmc3->AddFriend(tmc3_muid2);
    tmc3->AddFriend(tmc3_muid3);

    //Loop on categories A, B, C
    for(auto i = 0; i<ncat; i++){
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

    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}
