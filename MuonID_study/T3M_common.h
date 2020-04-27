#include <iostream>
#include <vector>

//TMVA Training options
    TString TMVA_outputpath = "MuonMVA_27april_"; //name to give to TMVA output files
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
    TString inputpath_bkg[] = {
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200427_1050/AnalysedTree_MC_2018BdToKK_muonid0.root", //Bd KK
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200427_1048/AnalysedTree_MC_2018BdToPiPi_muonid0.root", //Bd PiPi
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200427_1046/AnalysedTree_MC_2018BdToKPi_muonid0.root", //Bd KPi
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200427_1039/AnalysedTree_MC_2018BdToKPi_2_muonid0.root", //Bd KPi_2
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200427_0927/AnalysedTree_MC_2018BsToKK_muonid0.root", //Bs KK
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200427_0856/AnalysedTree_MC_2018BsToKK_2_muonid0.root", //Bs KK_2
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200427_1056/AnalysedTree_MC_2018BsToPiPi_muonid0.root", //Bs PiPi
    };

    TString inputpath_Ds = 
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200422_2320/AnalysedTree_MC_2018Ds_muonid_22april.root";
    TString inputpath_B0 = 
    "/lustre/cms/store/user/fsimone/MuonID/Analysis/20200422_2320/AnalysedTree_MC_2018B0_muonid_22april.root";
    TString inputpath_Bp = "";

//Coefficients for signal normalisation


//TMVA settings
// Variables
    std::vector<TString> var_names = {
                                      //Muon reconstruction
                                      "mu_combinedQuality_chi2LocalMomentum>2000?2000:mu_combinedQuality_chi2LocalMomentum",
                                      "mu_combinedQuality_chi2LocalPosition>500?500:mu_combinedQuality_chi2LocalPosition",
                                      "mu_combinedQuality_staRelChi2>500?500:mu_combinedQuality_staRelChi2",
                                      "mu_combinedQuality_trkRelChi2",
                                      "mu_combinedQuality_globalDeltaEtaPhi",
                                      "mu_combinedQuality_trkKink>400?400:mu_combinedQuality_trkKink",
                                      "mu_combinedQuality_glbKink",
                                      "mu_combinedQuality_glbTrackProbability>2000?2000:mu_combinedQuality_glbTrackProbability",

                                      //collection of hits in the HitPattern
                                      "mu_Numberofvalidtrackerhits", //Valid Tracker Hits
                                      "mu_Numberofvalidpixelhits",
                                      "mu_trackerLayersWithMeasurement",
                                      "mu_GLhitPattern_numberOfValidMuonHits",
                                      "mu_validMuonHitComb", //Hits in DT, CSC, RPC

                                      //muon track reconstruction
                                      "mu_numberOfMatchedStations",
                                      "mu_segmentCompatibility",
                                      "mu_timeAtIpInOut",
                                      "mu_timeAtIpInOutErr",

                                      //general track properties
                                      "mu_GLnormChi2>200?200:mu_GLnormChi2",
                                      "mu_innerTrack_normalizedChi2>200?200:mu_innerTrack_normalizedChi2",
                                      "mu_outerTrack_normalizedChi2>200?200:mu_outerTrack_normalizedChi2",
                                      "mu_innerTrack_validFraction", //Inner Valid Fraction
                                      "mu_QInnerOuter"
                                      //"", //dxyRef
                                      //"", //dzRef

                                      //custom variables track multiplicity
                                      };


    std::vector<TString> var_spec_names = {
                                      "abs(mu_eta)",
                                      "mu_pt",
                                      "mu_phi",
                                      "mu_simPdgId",
                                      "mu_simMotherPdgId",
                                      "mu_simType",
                                      "mu_isGlobal",
                                      "mu_SoftMVA"
                                      };
