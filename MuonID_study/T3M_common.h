#include <iostream>
#include <vector>

//TMVA Training options
    TString TMVA_outputpath = "MuonMVA_6april_"; //name to give to TMVA output files
    //change it to perform 5-fold Cross Validation
    bool doCV = false;
    TString method = "BDT";
    TString TMVA_weightfilename = "/weights/TMVA_new_BDT.weights.xml"; //name given training BDT in "normal" way
    
   // if(doCV)
   //TString method = "BDTG";
   //TString TMVA_weightfilename = "/weights/TMVACrossValidation_BDTG.weights.xml"; //name given training BDT with crossvalidation
    
//TMVA Evaluating options
    TString TMVA_inputpath = "_";  //name to load TMVA results for evaluation

//data rootfiles
    TString inputpath_BdToPiPi =
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200407_1707/AnalysedTree_MC_2018BdToPiPi_muonid0.root";
    TString inputpath_BdToKK =
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200407_1648/AnalysedTree_MC_2018BdToKK_muonid0.root";
    TString inputpath_Ds = 
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200407_1712/AnalysedTree_MC_2018Ds_merged_muonid.root";
    TString inputpath_B0 = "";
    TString inputpath_Bp = "";

//Coefficients for signal normalisation


//TMVA settings
// Variables
    std::vector<TString> var_names = {
                                      //Muon reconstruction
                                      "mu_combinedQuality_chi2LocalMomentum",
                                      "mu_combinedQuality_chi2LocalPosition",
                                      "mu_combinedQuality_staRelChi2",
                                      "mu_combinedQuality_trkRelChi2",
                                      "mu_combinedQuality_globalDeltaEtaPhi",
                                      "mu_combinedQuality_trkKink",
                                      "mu_combinedQuality_glbKink",
                                      "mu_combinedQuality_glbTrackProbability",

                                      //collection of hits in the HitPattern
                                      //"", //Valid Tracker Hits
                                      "mu_Numberofvalidpixelhits",
                                      "mu_trackerLayersWithMeasurement",
                                      "mu_GLhitPattern_numberOfValidMuonHits",
                                      //"", //Hits in DT, CSC, RPC

                                      //muon track reconstruction
                                      "mu_numberOfMatchedStations",
                                      "mu_segmentCompatibility",
                                      "mu_timeAtIpInOut",
                                      "mu_timeAtIpInOutErr",

                                      //general track properties
                                      "mu_GLnormChi2",
                                      "mu_innerTrack_normalizedChi2",
                                      "mu_outerTrack_normalizedChi2",
                                      //"", Inner Valid Fraction
                                      "mu_QInnerOuter"
                                      //"", //dxyRef
                                      //"", //dzRef

                                      //custom variables track multiplicity
                                      };


    std::vector<TString> var_spec_names = {
                                      "mu_eta",
                                      "mu_pt",
                                      "mu_simPdgId",
                                      "mu_simMotherPdgId",
                                      "mu_simType",
                                      "mu_isGlobal"
                                      };
