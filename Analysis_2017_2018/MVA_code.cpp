#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

using namespace TMVA;

void MVA_code(TString categ){
    //Check on input argument
    if(strcmp(categ, "A") != 0 && strcmp(categ, "B") != 0 && strcmp(categ, "C") != 0){
        cout << "Please choose between 'A', 'B' and 'C'" << endl;
        return;
    }
    
    // Output file
    TFile *fout = new TFile("TMVA_file.root", "RECREATE");
    
    // Get the signal and background trees from TFile source(s);
    // SIGNAL tree
    // Sgn 1 : DsTau3Mu
    TFile *f_sig1 = (TFile*)gROOT->GetListOfFiles()->FindObject("../AnalysedTree/AnalysedTree_cat_3loose_MC_sgn_DsTau3mu.root");
    if (!f_sig1 || !f_sig1->IsOpen()) {
        f_sig1 = new TFile("../AnalysedTree/AnalysedTree_cat_3loose_MC_sgn_DsTau3mu.root");
    }
    TString tree_sgn_Name1 = "FinalTree"; tree_sgn_Name1 += categ; tree_sgn_Name1 += "_sgn";
    TTree* sigTree1  = (TTree*)f_sig1->Get(tree_sgn_Name1);
    // Sgn 2 : B0Tau3Mu
    TFile *f_sig2 = (TFile*)gROOT->GetListOfFiles()->FindObject("../AnalysedTree/AnalysedTree_cat_3loose_MC_sgn_B0Tau3mu.root");
    if (!f_sig2 || !f_sig2->IsOpen()) {
        f_sig2 = new TFile("../AnalysedTree/DsTau3Mu_updated/AnalysedTree_cat_3loose_MC_sgn_B0Tau3mu.root");
    }
    TString tree_sgn_Name2 = "FinalTree"; tree_sgn_Name2 += categ; tree_sgn_Name2 += "_sgn";
    TTree* sigTree2  = (TTree*)f_sig2->Get(tree_sgn_Name2);
    // Sgn 3 : BpTau3Mu
    TFile *f_sig3 = (TFile*)gROOT->GetListOfFiles()->FindObject("../AnalysedTree/AnalysedTree_cat_3loose_MC_sgn_BpTau3mu.root");
    if (!f_sig3 || !f_sig3->IsOpen()) {
        f_sig3 = new TFile("../AnalysedTree/DsTau3Mu_updated/AnalysedTree_cat_3loose_MC_sgn_BpTau3mu.root");
    }
    TString tree_sgn_Name3 = "FinalTree"; tree_sgn_Name3 += categ; tree_sgn_Name3 += "_sgn";
    TTree* sigTree3  = (TTree*)f_sig3->Get(tree_sgn_Name3);
    
    // BACKGROUND trees
    // Bkg 1 : 2017B
    TFile *f_bkg1 = (TFile*)gROOT->GetListOfFiles()->FindObject("../AnalysedTree/AnalysedTree_cat_3loose_Data2017B_Bkg.root");
    if (!f_bkg1 || !f_bkg1->IsOpen()) {
        f_bkg1 = new TFile("../AnalysedTree/AnalysedTree_cat_3loose_Data2017B_Bkg.root");
    }
    TString tree_bkg1_Name = "FinalTree"; tree_bkg1_Name += categ; tree_bkg1_Name += "_Bkg";
    TTree* bkgTree1 = (TTree*)f_bkg1->Get(tree_bkg1_Name);
    // Bkg 2 : 2017D
    TFile *f_bkg2 = (TFile*)gROOT->GetListOfFiles()->FindObject("../AnalysedTree/AnalysedTree_cat_3loose_Data2017D_Bkg.root");
    if (!f_bkg2 || !f_bkg2->IsOpen()) {
        f_bkg2 = new TFile("../AnalysedTree/AnalysedTree_cat_3loose_Data2017D_Bkg.root");
    }
    TString tree_bkg2_Name = "FinalTree"; tree_bkg2_Name += categ; tree_bkg2_Name += "_Bkg";
    TTree* bkgTree2 = (TTree*)f_bkg2->Get(tree_bkg2_Name);
    // Bkg 3 : 2017C
//    TFile *f_bkg3 = (TFile*)gROOT->GetListOfFiles()->FindObject("../AnalysedTree/DsTau3Mu_updated/AnalysedTree_cat_3loose_Data2017C_Bkg.root");
//    if (!f_bkg3 || !f_bkg3->IsOpen()) {
//        f_bkg3 = new TFile("../AnalysedTree/AnalysedTree_cat_3loose_Data2017C_Bkg.root");
//    }
//    TString tree_bkg3_Name = "FinalTree"; tree_bkg3_Name += categ; tree_bkg3_Name += "_Bkg";
//    TTree* bkgTree3 = (TTree*)f_bkg3->Get(tree_bkg3_Name);
    // Bkg 4 : 2017F
    TFile *f_bkg4 = (TFile*)gROOT->GetListOfFiles()->FindObject("../AnalysedTree/AnalysedTree_cat_3loose_Data2017F_Bkg.root");
    if (!f_bkg4 || !f_bkg4->IsOpen()) {
        f_bkg4 = new TFile("../AnalysedTree/AnalysedTree_cat_3loose_Data2017F_Bkg.root");
    }
    TString tree_bkg4_Name = "FinalTree"; tree_bkg4_Name += categ; tree_bkg4_Name += "_Bkg";
    TTree* bkgTree4 = (TTree*)f_bkg4->Get(tree_bkg4_Name);
    
    // Set the event weights per tree
    Double_t sigWeight1  = 1.0;
    Double_t sigWeight2  = 0.1;
    Double_t sigWeight3  = 0.1;
    Double_t bkgWeight1 = 1.0;
    Double_t bkgWeight2 = 1.0;
    Double_t bkgWeight3 = 1.0;
    Double_t bkgWeight4 = 1.0;
    
    Factory *factory = new Factory("TMVA_new", fout, "");
    DataLoader *dataloader = new DataLoader("dataset");
    
    // Register the trees
    dataloader->AddSignalTree(sigTree1, sigWeight1);
//    dataloader->AddSignalTree(sigTree2, sigWeight2);
//    dataloader->AddSignalTree(sigTree3, sigWeight3);
    dataloader->AddBackgroundTree(bkgTree1, bkgWeight1);
    dataloader->AddBackgroundTree(bkgTree2, bkgWeight2);
//    dataloader->AddBackgroundTree(bkgTree3, bkgWeight3);
    dataloader->AddBackgroundTree(bkgTree4, bkgWeight4);
    
    // Variables declaration
//    dataloader->AddVariable("Pmu3", 'F');
//    dataloader->AddVariable("cLP", 'F');
//    dataloader->AddVariable("tKink", 'F');
    dataloader->AddVariable("segmComp", 'F'); // segmComp MIN
    dataloader->AddVariable("fv_nC", 'F'); // Chi2 norm vertex
    dataloader->AddVariable("fv_dphi3D", 'F'); // Angle PVSV
    dataloader->AddVariable("fv_d3Dsig", 'F'); // Flight dist significance
//    dataloader->AddVariable("d0", 'F');
    dataloader->AddVariable("d0sig", 'F'); // Transverse IP significance
//    dataloader->AddVariable("mindca_iso", 'F'); // isolation MIN Distance of Closest Approach
    dataloader->AddVariable("caloComp", 'F'); // calo Comp MIN
//    dataloader->AddVariable("trkRel", 'F');
//    dataloader->AddVariable("Pmu1", 'F');
//    dataloader->AddVariable("Ptmu1", 'F');
//    dataloader->AddVariable("Etamu1", 'F');
//    dataloader->AddVariable("Pmu2", 'F');
//    dataloader->AddVariable("Ptmu2", 'F');
//    dataloader->AddVariable("Etamu2", 'F');
//    dataloader->AddVariable("Ptmu3", 'F');
//    dataloader->AddVariable("Etamu3", 'F');
//    dataloader->AddVariable("P_tripl", 'F');
//    dataloader->AddVariable("Pt_tripl", 'F');
//    dataloader->AddVariable("Eta_tripl", 'F');
//    dataloader->AddVariable("nStMu1", 'F');
//    dataloader->AddVariable("nStMu2", 'F');
//    dataloader->AddVariable("nStMu3", 'F');
//    dataloader->AddVariable("Iso03Mu1", 'F');
//    dataloader->AddVariable("Iso03Mu2", 'F');
//    dataloader->AddVariable("Iso03Mu3", 'F');
//    dataloader->AddVariable("Iso05Mu1", 'F');
//    dataloader->AddVariable("Iso05Mu2", 'F');
//    dataloader->AddVariable("Iso05Mu3", 'F');
//    dataloader->AddVariable("nMatchesMu1", 'F');
//    dataloader->AddVariable("nMatchesMu2", 'F');
//    dataloader->AddVariable("nMatchesMu3", 'F');
//    dataloader->AddVariable("timeAtIpInOut1", 'F');
//    dataloader->AddVariable("timeAtIpInOut2", 'F');
//    dataloader->AddVariable("timeAtIpInOut3", 'F');
//    dataloader->AddVariable("cQ_uS", 'F');
//    dataloader->AddVariable("cQ_tK", 'F');
//    dataloader->AddVariable("cQ_gK", 'F');
//    dataloader->AddVariable("cQ_tRChi2", 'F');
//    dataloader->AddVariable("cQ_sRChi2", 'F');
//    dataloader->AddVariable("cQ_Chi2LM", 'F');
//    dataloader->AddVariable("cQ_Chi2lD", 'F');
//    dataloader->AddVariable("cQ_gDEP", 'F');
//    dataloader->AddVariable("cQ_tM", 'F');
//    dataloader->AddVariable("cQ_gTP", 'F');
//    dataloader->AddVariable("calEn_emMu1", 'F');
//    dataloader->AddVariable("calEn_emMu2", 'F');
//    dataloader->AddVariable("calEn_emMu3", 'F');
//    dataloader->AddVariable("calEn_hadMu1", 'F');
//    dataloader->AddVariable("calEn_hadMu2", 'F');
//    dataloader->AddVariable("calEn_hadMu3", 'F');
//    dataloader->AddVariable("fliDistPVSV_Chi2", 'F');
//    dataloader->AddVariable("isGlb3", 'F');
//    dataloader->AddVariable("isTracker3", 'F');
//    dataloader->AddVariable("isLoose3", 'F');
//    dataloader->AddVariable("isSoft3", 'F');
//    dataloader->AddVariable("isPF3", 'F');
//    dataloader->AddVariable("isRPC3", 'F');
//    dataloader->AddVariable("isSA3", 'F');
//    dataloader->AddVariable("isCalo3", 'F');
//    dataloader->AddVariable("Vx1", 'F');
//    dataloader->AddVariable("Vx2", 'F');
//    dataloader->AddVariable("Vx3", 'F');
//    dataloader->AddVariable("Vy1", 'F');
//    dataloader->AddVariable("Vy2", 'F');
//    dataloader->AddVariable("Vy3", 'F');
//    dataloader->AddVariable("Vz1", 'F');
//    dataloader->AddVariable("Vz2", 'F');
//    dataloader->AddVariable("Vz3", 'F');
//    dataloader->AddVariable("RefVx1", 'F');
//    dataloader->AddVariable("RefVy1", 'F');
//    dataloader->AddVariable("RefVz1", 'F');
//    dataloader->AddVariable("SVx", 'F');
//    dataloader->AddVariable("SVy", 'F');
//    dataloader->AddVariable("SVz", 'F');
//    dataloader->AddVariable("had03", 'F');
//    dataloader->AddVariable("had05", 'F');
//    dataloader->AddVariable("nJets03", 'F');
//    dataloader->AddVariable("nJets05", 'F');
//    dataloader->AddVariable("nTracks03", 'F');
//    dataloader->AddVariable("nTracks05", 'F');
//    dataloader->AddVariable("sumPt03", 'F');
//    dataloader->AddVariable("sumPt05", 'F');
//    dataloader->AddVariable("hadVeto03", 'F');
//    dataloader->AddVariable("hadVeto05", 'F');
//    dataloader->AddVariable("emVeto03", 'F');
//    dataloader->AddVariable("emVeto05", 'F');
//    dataloader->AddVariable("trVeto03", 'F');
//    dataloader->AddVariable("trVeto05", 'F');
    

    TCut preselectionCut = "";
    dataloader->PrepareTrainingAndTestTree(preselectionCut, "SplitMode=Alternate");
    //    factory->PrepareTrainingAndTestTree(preselectionCut, "SplitMode=Random:NormMode=NumEvents");
    
    // Booking of MVA methods : BDT
    factory->BookMethod(dataloader, Types::kBDT, "BDT", "");
    
    // Training the MVA methods
    factory->TrainAllMethods();
    
    // Testing the MVA methods
    factory->TestAllMethods();
    
    // Evaluating the MVA methods
    factory->EvaluateAllMethods();
    
    // Save the output
    fout->Close();
    
    delete factory;
    delete dataloader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVAGui("TMVA_file.root");
    
}
