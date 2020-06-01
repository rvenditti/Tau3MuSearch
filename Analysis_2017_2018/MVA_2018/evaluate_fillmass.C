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

    //Loop on input Tree
    while (treeReader.Next()) {
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
        auto reso = *reader_triplMassReso;
        auto mu1pt = *reader_Ptmu1;
        if(((category.Contains("A") && reso < 0.007 ) ||
            (category.Contains("B") && reso>= 0.007 && reso <= 0.0105 ) ||
            (category.Contains("C") && reso > 0.0105 )) ){
            if(((category.Contains("high") && mu1pt >= 7 ) ||
                (category.Contains("low")  && mu1pt <  7 ) ||
                (!(category.Contains("low")) && !(category.Contains("low"))))){
                auto BDTscore = MVAreader->EvaluateMVA(method);
                auto mass = *reader_triplMass; 
                if( BDTscore >= BDTcutvalues.a ) //category 1
                    hTriplMass1->Fill(mass);
                else if ( BDTscore < BDTcutvalues.a && BDTscore >= BDTcutvalues.b) //category 2
                    hTriplMass2->Fill(mass);
            }
        }
    }
    delete MVAreader;
}


void evaluate_fillmass(bool merge_high_low)
{
    int ncat = sizeof(cat_label)/sizeof(cat_label[0]);

    TString run[] = {"A", "B", "C", "D"};
    int nrun = sizeof(inputpath_datarun)/sizeof(inputpath_datarun[0]);
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
        TString fout_path = "";
        if(merge_high_low) fout_path = "datacardT3Mu_"+TMVA_inputpath+category+"_noptlabel.root";
        else fout_path = "datacardT3Mu_"+TMVA_inputpath+category+".root";
        TFile *fout = new TFile(fout_path,"recreate");
        fout->cd();
        TString cat_histo = category;
        if(merge_high_low) { cat_histo =cat_histo.ReplaceAll("lowpt",""); cat_histo = cat_histo.ReplaceAll("highpt",""); }
        TH1F * hTriplMass_data1 = new TH1F ("data_obs"+cat_histo+"1","Triplet mass "+cat_histo+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_data2 = new TH1F ("data_obs"+cat_histo+"2","Triplet mass "+cat_histo+"2",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Ds1   = new TH1F ("signalDs"+cat_histo+"1","Triplet mass "+cat_histo+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Ds2   = new TH1F ("signalDs"+cat_histo+"2","Triplet mass "+cat_histo+"2",42, 1.600000, 2.020000);
        TH1F * hTriplMass_B01   = new TH1F ("signalB0"+cat_histo+"1","Triplet mass "+cat_histo+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_B02   = new TH1F ("signalB0"+cat_histo+"2","Triplet mass "+cat_histo+"2",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Bp1   = new TH1F ("signalBp"+cat_histo+"1","Triplet mass "+cat_histo+"1",42, 1.600000, 2.020000);
        TH1F * hTriplMass_Bp2   = new TH1F ("signalBp"+cat_histo+"2","Triplet mass "+cat_histo+"2",42, 1.600000, 2.020000);
        
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
        TH1F * hSignal1 = new TH1F ("signal"+cat_histo+"1","Triplet mass "+cat_histo+"1",42, 1.600000, 2.020000);
        TH1F * hSignal2 = new TH1F ("signal"+cat_histo+"2","Triplet mass "+cat_histo+"2",42, 1.600000, 2.020000);        

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
        hTriplMass_data1->Write("background"+cat_histo+"1");
        hTriplMass_data2->Write("background"+cat_histo+"2");
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

    cout<<"Exiting ROOT"<<endl;
    gApplication->Terminate();
    return 0;
}
