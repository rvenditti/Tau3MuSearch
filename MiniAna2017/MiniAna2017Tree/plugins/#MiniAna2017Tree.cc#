// -*- C++ -*-
//
// Package:    MiniAna2017/MiniAna2017Tree
// Class:      MiniAna2017Tree
//
/**\class MiniAna2017Tree MiniAna2017Tree.cc MiniAna2017/MiniAna2017Tree/plugins/MiniAna2017Tree.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  venditti
//         Created:  Tue, 18 Dec 2018 09:30:06 GMT
//
//


// system include files
#include <memory>
#include <algorithm>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TFile.h"
#include "TH1.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"


#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicTree.h"


#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"


#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1TGlobalPrescalesVetosRcd.h"
#include "CondFormats/L1TObjects/interface/L1TGlobalPrescalesVetos.h"
////
class MiniAna2017Tree : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
    explicit MiniAna2017Tree(const edm::ParameterSet&);
    ~MiniAna2017Tree();
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    float dR(float eta1, float eta2, float phi1, float phi2);
    float dRtriggerMatch(pat::Muon m, trigger::TriggerObjectCollection triggerObjects);
    void beginRun(edm::Run const &, edm::EventSetup const&, edm::Event const&);
    
    
private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    // virtual void beginRun(edm::Run const &, edm::EventSetup const&) override;
    virtual void endJob() override;
    edm::EDGetTokenT<edm::View<pat::Muon> > muons_;
    edm::EDGetTokenT<edm::View<reco::Vertex> > vertex_;
    edm::EDGetTokenT<edm::View<reco::Track> > trackToken_;
    edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > Cand3Mu_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticles_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puToken_ ;
    edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
    edm::EDGetTokenT<trigger::TriggerEvent> trigeventToken_;
    edm::EDGetToken algToken_;
    bool isMc;
    bool isAna;
    //  edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
    //  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjCollToken;
    const TransientTrackBuilder* theTransientTrackBuilder_;
    HLTConfigProvider hltConfig;
    //  TPMERegexp* _re;
    
    edm::Service<TFileService> fs;
    l1t::L1TGlobalUtil* gtUtil_;
  TH1F *hEvents;
  TH1F *hEventsAfterGoodCand;
  //      TH1F *hEventsAfterGoodCand;
  //    TH1F *hEventsAfterGoodCand;
   edm::InputTag algInputTag_;

    //tree
    TTree*      tree_;
    std::vector<float>  MuonPt, MuonEta, MuonPhi;
    std::vector<double> MuonEnergy,  MuonCharge;
    
    std::vector<int> GenParticle_PdgId, GenParticle_MotherPdgId;
    std::vector<double> GenParticle_Pt, GenParticle_Eta,    GenParticle_Phi;
    
    //Vtx position
    std::vector<double>  Muon_vx,  Muon_vy,  Muon_vz;
    
    //MuonID
    std::vector<double>  Muon_isGlobal,  Muon_isTracker,  Muon_isSoft,  Muon_isLoose, Muon_isTight,  Muon_isPF,  Muon_isRPCMuon,  Muon_isStandAloneMuon,  Muon_isTrackerMuon,  Muon_isCaloMuon,  Muon_isQualityValid,  Muon_isTimeValid,  Muon_isIsolationValid,  Muon_numberOfMatchedStations,  Muon_numberOfMatches;
    
    //MuonTime
    std::vector<double>  Muon_timeAtIpInOut,Muon_timeAtIpInOutErr;
    
    //Muon inner + outer track
    std::vector<double>  Muon_GLnormChi2, Muon_GLhitPattern_numberOfValidMuonHits,  Muon_trackerLayersWithMeasurement,  Muon_Numberofvalidpixelhits,  Muon_outerTrack_p,  Muon_outerTrack_eta,
    Muon_outerTrack_phi,  Muon_outerTrack_normalizedChi2,  Muon_outerTrack_muonStationsWithValidHits,  Muon_innerTrack_p,  Muon_innerTrack_eta,  Muon_innerTrack_phi,  Muon_innerTrack_normalizedChi2,  Muon_QInnerOuter;
    
    std::vector<double>   Muon_combinedQuality_updatedSta,  Muon_combinedQuality_trkKink,  Muon_combinedQuality_glbKink,  Muon_combinedQuality_trkRelChi2,  Muon_combinedQuality_staRelChi2,  Muon_combinedQuality_chi2LocalPosition,  Muon_combinedQuality_chi2LocalMomentum,  Muon_combinedQuality_localDistance,  Muon_combinedQuality_globalDeltaEtaPhi,  Muon_combinedQuality_tightMatch,  Muon_combinedQuality_glbTrackProbability,  Muon_calEnergy_em,  Muon_calEnergy_emS9,  Muon_calEnergy_emS25,  Muon_calEnergy_had,  Muon_calEnergy_hadS9,  Muon_segmentCompatibility,  Muon_caloCompatibility,  Muon_ptErrOverPt, Muon_BestTrackPt,  Muon_BestTrackPtErr, Muon_BestTrackEta,  Muon_BestTrackEtaErr,  Muon_BestTrackPhi,  Muon_BestTrackPhiErr;

    std::vector<int>  Muon_simPdgId, Muon_simMotherPdgId, Muon_simFlavour,  Muon_simType, Muon_simBX, Muon_simTpEvent, Muon_simMatchQuality;
    std::vector<double>  Mu1_Pt,  Mu1_Eta,  Mu1_Phi,  Mu2_Pt,  Mu2_Eta,  Mu2_Phi,  Mu3_Pt,  Mu3_Eta,  Mu3_Phi, GenMatchMu1_SimPt, GenMatchMu2_SimPt, GenMatchMu3_SimPt,GenMatchMu1_SimEta, GenMatchMu2_SimEta, GenMatchMu3_SimEta, GenMatchMu1_SimPhi, GenMatchMu2_SimPhi, GenMatchMu3_SimPhi,  GenMatchMu1_Pt,  GenMatchMu2_Pt,  GenMatchMu3_Pt,  GenMatchMu1_Eta,  GenMatchMu2_Eta,  GenMatchMu3_Eta,  GenMatchMu1_Phi,  GenMatchMu2_Phi,  GenMatchMu3_Phi;
    std::vector<float> Mu1_dRtriggerMatch, Mu2_dRtriggerMatch, Mu3_dRtriggerMatch;
 
    std::vector<double>     Muon_emEt03, Muon_hadEt03, Muon_nJets03, Muon_nTracks03, Muon_sumPt03, Muon_emEt05,    Muon_hadEt05, Muon_nJets05, Muon_nTracks05, Muon_sumPt05,
    Muon_hadVetoEt03,Muon_emVetoEt03,    Muon_trackerVetoPt03,    Muon_hadVetoEt05,    Muon_emVetoEt05,    Muon_trackerVetoPt05;
    //dd  Mu1_SimPt,  Mu1_SimEta,  Mu1_SimPhi,  Mu2_SimPt,  Mu2_SimEta,  Mu2_SimPhi, Mu3_SimPt,  Mu3_SimEta,  Mu3_SimPhi,
    
    std::vector<double>     Triplet_mindca_iso, Triplet_relativeiso;
 
    std::vector<int>  Mu1_TripletIndex,  Mu2_TripletIndex,  Mu3_TripletIndex;
    std::vector<int>  Mu1_NTracks03iso,  Mu2_NTracks03iso,  Mu3_NTracks03iso;
    
    int TripletCollectionSize, PVCollection_Size, MuonCollectionSize;
    std::vector<double>  TripletVtx_x,  TripletVtx_y,  TripletVtx_z,  TripletVtx_Chi2,  TripletVtx_NDOF,  Triplet_Mass,  Triplet_Pt,  Triplet_Eta,  Triplet_Phi, Triplet_Charge;
    
    std::vector<double> dxy_mu1, dxy_mu2, dxy_mu3, dxyErr_mu1, dxyErr_mu2, dxyErr_mu3; 
    
    std::vector<double>  RefittedPV_x;
    std::vector<double>  RefittedPV_y;
    std::vector<double>  RefittedPV_z;
    std::vector<double>  RefittedPV_NTracks;
    std::vector<int>     RefittedPV_isValid;
    
    //RefittedPV_Chi2.push_back(PVertex.);
    
    std::vector<double>  FlightDistPVSV;
    std::vector<double>  FlightDistPVSV_Err;
    std::vector<double>  FlightDistPVSV_Significance;
    std::vector<double>  FlightDistPVSV_chi2;
    
   std::vector<double> PV_x,  PV_y,  PV_z,  PV_NTracks;
  std::vector<int> NGoodTriplets;
    uint  evt, run, lumi, puN;
  std::vector<string>  Trigger_l1name;
  std::vector<int> Trigger_l1decision;
  std::vector<int> Trigger_l1prescale;

  std::vector<string>  Trigger_hltname;
  std::vector<int> Trigger_hltdecision;
    //SyncTree
    /*  TTree*      SyncTree_;
     std::vector<float>  allmuons_pt, leadmuon_pt, leadmuon_phi, leadmuon_eta;
     std::vector<float>  alltracks_pt, leadtrack_pt,  leadtrack_eta,  leadtrack_phi;
     uint nprimevtxs, nmuons, evt, run, lumi;
     */
    };
    
    
    
    MiniAna2017Tree::MiniAna2017Tree(const edm::ParameterSet& iConfig){
        edm::InputTag algInputTag_;
        isMc = iConfig.getUntrackedParameter<bool>("isMcLabel");
	isAna = iConfig.getUntrackedParameter<bool>("isAnaLabel");
        muons_ = consumes<edm::View<pat::Muon> >  (iConfig.getParameter<edm::InputTag>("muonLabel"));
        vertex_ = consumes<edm::View<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("VertexLabel"));
        trackToken_ = consumes<edm::View<reco::Track> > (edm::InputTag("generalTracks"));
        genParticles_ = consumes<edm::View<reco::GenParticle>  > (iConfig.getParameter<edm::InputTag>("genParticleLabel"));
        Cand3Mu_ = consumes<edm::View<reco::CompositeCandidate> > (iConfig.getParameter<edm::InputTag>("Cand3MuLabel"));
        puToken_ =   consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummary"));
	triggerToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
	trigeventToken_ = consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerSummary"));
	algToken_ = consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("AlgInputTag"));
	gtUtil_ = new l1t::L1TGlobalUtil(iConfig, consumesCollector(), *this, algInputTag_, algInputTag_);
        //  _hltInputTag(iConfig.getParameter<edm::InputTag>("hltInputTag")),
        //tauToken_(consumes(iConfig.getParameter("taus"))),
        //metToken_(consumes(iConfig.getParameter("mets")))
        //tree_(0);
        //MuonPt(0);
    }
    MiniAna2017Tree::~MiniAna2017Tree()
    {
        
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
    }
    
    
    float MiniAna2017Tree::dR(float eta1, float eta2, float phi1, float phi2){
        float dphi=(phi1-phi2);
        float deta=(eta1-eta2);
        float deltaR= TMath::Sqrt(dphi*dphi + deta*deta);
        return deltaR;
    }

    float MiniAna2017Tree::dRtriggerMatch(pat::Muon m, trigger::TriggerObjectCollection triggerObjects) {
        float dRmin = 1.;                                                                          
        for (unsigned int i = 0 ; i < triggerObjects.size() ; i++) {
            float deltaR = sqrt( reco::deltaR2(triggerObjects[i].eta(), triggerObjects[i].phi(), m.eta(), m.phi()));
            //float deltaR = sqrt( pow(triggerObjects[i].eta() - m.eta(), 2) + pow(acos(cos(triggerObjects[i].phi() - m.phi())), 2));
            if (deltaR < dRmin) dRmin = deltaR;                                                    
        }                                                                                          
        return dRmin;                                                                              
    }

    
    bool isGoodTrack(const reco::Track &track) {
        if(track.pt()>1){
            if(std::fabs(track.eta())<2.4){
                if(track.hitPattern().trackerLayersWithMeasurement()>5){
                    if(track.hitPattern().pixelLayersWithMeasurement()>1) return true;
                }
            }
        }
        return false;
    }
    
    typedef std::map<const reco::Track*, reco::TransientTrack> TransientTrackMap;
    // auxiliary function to exclude tracks associated to tau lepton decay "leg"
    // from primary event vertex refit
    bool tracksMatchByDeltaR(const reco::Track* trk1, const reco::Track* trk2)
    {
        if ( reco::deltaR(*trk1, *trk2) < 1.e-2 && trk1->charge() == trk2->charge() ) return true;
        else return false;
    }
    
    void removeTracks(TransientTrackMap& pvTracks_toRefit, const std::vector<reco::Track*> svTracks)
    {
        for ( std::vector<reco::Track*>::const_iterator svTrack = svTracks.begin(); svTrack != svTracks.end(); ++svTrack ){
            //--- remove track from list of tracks included in primary event vertex refit
            //    if track matches by reference or in eta-phi
            //    any of the tracks associated to tau lepton decay "leg"
            for ( TransientTrackMap::iterator pvTrack = pvTracks_toRefit.begin(); pvTrack != pvTracks_toRefit.end(); ++pvTrack ) {
                if ( tracksMatchByDeltaR(pvTrack->first, *svTrack) ) {
                    pvTracks_toRefit.erase(pvTrack);
                    break;
                }
            }
        }
    }
    
    
    
    void MiniAna2017Tree::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup, const edm::Event& iEvent) {
      /*  
      edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
        edm::InputTag trigResultsTag("TriggerResults"," ","HLT"); //make sure have correct process on MC
        iEvent.getByLabel(trigResultsTag,trigResults);
        
        bool changed = true;
        if (hltConfig.init(iRun, iSetup, trigResultsTag.process(), changed)) {
            // if init returns TRUE, initialisation has succeeded!
            std::cout << "HLT config with process name "<< trigResultsTag.process() << " successfully extracted"<<std::endl;
        }
        else {
            // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
            // with the file and/or code and needs to be investigated!
            std::cout << "Error! HLT config extraction with process name " << trigResultsTag.process() << " failed"<<std::endl;
            // In this case, all access methods will return empty values!
        }
      */
    }
    
    
    
    void
    MiniAna2017Tree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
    {
        using namespace edm;
        using namespace std;
        using namespace reco;
        using std::vector;
        
        
        
        edm::Handle< edm::View<reco::Vertex> >vertices;
        iEvent.getByToken(vertex_, vertices);
        if (vertices->empty()) return; // skip the event if no PV found
        const reco::Vertex &PV = vertices->front();
        
        edm::Handle< edm::View<pat::Muon> > muons;
        iEvent.getByToken(muons_, muons);
        
        
        edm::Handle<edm::View<reco::CompositeCandidate> > Cand3Mu;
        iEvent.getByToken(Cand3Mu_, Cand3Mu);
        
        
        edm::Handle< edm::View<reco::GenParticle> > genParticles;
        iEvent.getByToken(genParticles_, genParticles);
        
        
        edm::Handle<edm::View<reco::Track> > trackCollection;
        iEvent.getByToken(trackToken_, trackCollection);
        
	Handle<TriggerResults> triggerResults;
	iEvent.getByToken(triggerToken_, triggerResults);

	Handle<trigger::TriggerEvent> triggerSummary;
	iEvent.getByToken(trigeventToken_, triggerSummary);

	Handle<BXVector<GlobalAlgBlk>> alg;
	iEvent.getByToken(algToken_,alg);

	edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);
        theTransientTrackBuilder_ = theTransientTrackBuilder.product();

	//edm::Handle<SimTrackContainer> simTracks;
	//iEvent.getByLabel("g4SimHits",simTracks);

        hEvents->Fill(1);


	///////////////Fill Trigger Vars, L1 and HLT///////////////

	gtUtil_->retrieveL1(iEvent, iSetup, algToken_);
	const vector<pair<string, bool> > initialDecisions = gtUtil_->decisionsInitial();

	if (!iEvent.isRealData())
	  {
	    for (size_t i_l1t = 0; i_l1t < initialDecisions.size(); i_l1t++) 
	      {
		string l1tName = (initialDecisions.at(i_l1t)).first;
		if(l1tName.find("DoubleMu") != string::npos || l1tName.find("TripleMu") != string::npos){
		  Trigger_l1name.push_back( l1tName );
		  Trigger_l1decision.push_back( initialDecisions.at(i_l1t).second );
		  Trigger_l1prescale.push_back( 1 );
		}
	      }
	  }
	else
	  {
	    ESHandle<L1TGlobalPrescalesVetos> psAndVetos;
	    auto psRcd = iSetup.tryToGet<L1TGlobalPrescalesVetosRcd>();
	    if(psRcd) psRcd->get(psAndVetos);
	    int columnN= gtUtil_->prescaleColumn();
	    for (size_t i_l1t = 0; i_l1t < initialDecisions.size(); i_l1t++) {
	      string l1tName = (initialDecisions.at(i_l1t)).first;
	      if(l1tName.find("DoubleMu") != string::npos || l1tName.find("TripleMu") != string::npos){
		//cout<<"L1Seed="<<l1tName<<" decision="<<initialDecisions.at(i_l1t).second<<" prescale="<<(psAndVetos->prescale_table_)[columnN][i_l1t]<<endl;
		Trigger_l1name.push_back( l1tName );
		Trigger_l1decision.push_back( initialDecisions.at(i_l1t).second );
		Trigger_l1prescale.push_back( (psAndVetos->prescale_table_)[columnN][i_l1t]);
                }
            }
        }
        
            
	const TriggerNames &triggerNames = iEvent.triggerNames( *triggerResults );
	for (size_t i_hlt = 0; i_hlt != triggerResults->size(); ++i_hlt){
	    string hltName = triggerNames.triggerName(i_hlt);
	    if(hltName.find("HLT_DoubleMu") != string::npos){
	      //cout<<" HLTPath="<<hltName<<" isPassed="<<triggerResults->accept(i_hlt )<<endl;
	      Trigger_hltname.push_back(hltName);
	      Trigger_hltdecision.push_back(triggerResults->accept(i_hlt ));
	    }
	}
    
 
	///////////////Fill GenParticles///////////////
        if(isMc){
            uint j=0;
            uint ngenP=genParticles->size();
            std::vector<int> genPidx;
            //uint tauRaw=-999; //int DsRaw=-999;
            
	    //            cout<<"****************GenLevel Info Begin********************"<<endl;
            for(edm::View<reco::GenParticle>::const_iterator gp=genParticles->begin(); gp!=genParticles->end(), j<ngenP; ++gp , ++j){
                //if(fabs(gp->pdgId())==15) tauRaw = j;
                
	      if(fabs(gp->pdgId())==13 || fabs(gp->pdgId())==15  || fabs(gp->pdgId())==11 || fabs(gp->pdgId())==211 || fabs(gp->pdgId())==321 ||  fabs(gp->pdgId())==12  || fabs(gp->pdgId())==14 || fabs(gp->pdgId())==16 || fabs(gp->pdgId())==431 || fabs(gp->pdgId())==511 || fabs(gp->pdgId())==521) {
                    GenParticle_PdgId.push_back(gp->pdgId());
                    GenParticle_Pt.push_back(gp->pt());
                    GenParticle_Eta.push_back(gp->eta());
                    GenParticle_Phi.push_back(gp->phi());
                    if (gp->numberOfMothers()) {GenParticle_MotherPdgId.push_back(gp->mother(0)->pdgId());
                    }else{
                        GenParticle_MotherPdgId.push_back(-99);
                    }
                    for (uint i=0; i<gp->numberOfMothers();i++){
                        if(fabs(gp->mother(i)->pdgId())==15) {
			  //std::cout<<j<<"--genMu pt="<<gp->pt()<<" eta="<<gp->eta()<<" phi="<<gp->phi()<<" pdgID="<<gp->pdgId()<<" tau pt="<<gp->mother(i)->pt()<<" mu vtx_x="<<gp->vx()<<" mu vtx_y="<<gp->vy()<<" mu vtx_z="<<gp->vz()<<endl;
                            //cout<<tauRaw<<"--mother pdgID="<<gp->mother(i)->pdgId()<<" mother vtx_x="<<gp->mother(i)->vx()<<" vy="<<gp->mother(i)->vy()<<" vz="<<gp->mother(i)->vz()<<endl;
                            //cout<<"TauMother pdgId="<<gp->mother(i)->mother(0)->pdgId()<<" vx="<<gp->mother(i)->mother(0)->vx()<<" vy="<<gp->mother(i)->mother(0)->vy()<<" vz="<<gp->mother(i)->mother(0)->vz()<<endl;
                            genPidx.push_back(j);
                        }
                    }
                }
            }
            
	    //            cout<<"****************GenLevel Info End ********************"<<endl;
        }
        ///////////////Fill GenParticles///////////////
        
        //Primary Vtx
	//  const reco::Vertex* eventVertex;

        PVCollection_Size = vertices->size();
        //cout<<" PV size ="<<vertices->size()<<endl;

	//std::vector<TransientTrackMap> pvTrackMap_refitVec;	
	//for (reco::VertexCollection::const_iterator it = vertices.begin(); it != vertices.end() && VtxIt != vertices->size() ; ++it, ++VtxIt) {
	//for(uint VtxIt =0;VtxIt<vertices->size();VtxIt++ ){
	//std::vector<reco::TransientTrack> pvTracks_original2;
	//TransientTrackMap pvTrackMap_refit2;
	//for ( reco::Vertex::trackRef_iterator pvTrack =  (*vertices)[VtxIt].tracks_begin(); pvTrack != (*vertices)[VtxIt].tracks_end(); ++pvTrack ) {
	    //cout<<" pv track size"<<(*vertices)[VtxIt].tracks_size()<<endl;
	    //reco::TransientTrack pvTrack_transient =theTransientTrackBuilder_->build(pvTrack->get());
            //pvTracks_original2.push_back(pvTrack_transient);
            //pvTrackMap_refit2.insert(std::make_pair(pvTrack->get(), pvTrack_transient));
	    //pvTrackMap_refitVec.push_back(pvTrackMap_refit2);
	//}
	//cout<<"Vtx id="<<VtxIt<<" Number of tracks associated to the PV="<<pvTracks_original2.size()<<endl;
	//}
	
	
        //Triplets  Loop
        //cout<<"Number Of Triplets="<<Cand3Mu->size()<<endl;
	std::vector<int> NTripl;
	if(isAna){
        TripletCollectionSize = Cand3Mu->size() ;
        int TripletIndex =-99; uint trIn=0;
        for(edm::View<reco::CompositeCandidate>::const_iterator TauIt=Cand3Mu->begin(); TauIt!=Cand3Mu->end(), trIn<Cand3Mu->size(); ++TauIt, ++trIn){

	  const Candidate * c1 = TauIt->daughter(0)->masterClone().get();
	  const pat::Muon *mu1 = dynamic_cast<const pat::Muon *>(c1);
	  
	  const Candidate * c2 = TauIt->daughter(1)->masterClone().get();
	  const pat::Muon *mu2 = dynamic_cast<const pat::Muon *>(c2);
		
	  const Candidate * c3 = TauIt->daughter(2)->masterClone().get();
	  const pat::Muon *mu3 = dynamic_cast<const pat::Muon *>(c3);



            TrackRef trk1, trk2, trk3;
            if (mu1->isGlobalMuon()) { trk1 = mu1->get<TrackRef,reco::CombinedMuonTag>();}
            else { trk1 = mu1->get<TrackRef>();}
            if (mu2->isGlobalMuon()) { trk2 = mu2->get<TrackRef,reco::CombinedMuonTag>();}
            else{ trk2 = mu2->get<TrackRef>();}
            if (mu3->isGlobalMuon()) { trk3 = mu3->get<TrackRef,reco::CombinedMuonTag>();}
            else{  trk3 = mu3->get<TrackRef>();}
            //cout<<" trk1 id="<<trk1.id()<<" tr2:"<<trk2.id()<<" trk3="<<trk3.id()<<endl;
            const reco::TransientTrack transientTrack1=theTransientTrackBuilder_->build( trk1 );
            const reco::TransientTrack transientTrack2=theTransientTrackBuilder_->build( trk2 );
            const reco::TransientTrack transientTrack3=theTransientTrackBuilder_->build( trk3 );
            reco::Track Track1 =transientTrack1.track();
            reco::Track Track2 =transientTrack2.track();
            reco::Track Track3 =transientTrack3.track();
            reco::Track* TrackRef1=&Track1;
            reco::Track* TrackRef2=&Track2;
            reco::Track* TrackRef3=&Track3;
            vector<reco::Track*> SVTrackRef;
            SVTrackRef.push_back(TrackRef1);
            SVTrackRef.push_back(TrackRef2);
            SVTrackRef.push_back(TrackRef3);



            reco::Vertex TripletVtx = reco::Vertex(TauIt->vertex(), TauIt->vertexCovariance(), TauIt->vertexChi2(), TauIt->vertexNdof(), TauIt->numberOfDaughters() );
	    double dphi_pv = -1.0;
	    uint primaryvertex_index=0;
	    TLorentzVector ThreeCandidate;
	    ThreeCandidate.SetPtEtaPhiM(TauIt->pt(), TauIt->eta(), TauIt->phi(), TauIt->mass());


	    for(uint VtxIt =0;VtxIt<vertices->size();VtxIt++ ){
	      //cout<<"Vtx id="<<VtxIt<<" x="<<(*vertices)[VtxIt].x()<<endl;
	      TVector3 Dv3D_reco(TripletVtx.x() - (*vertices)[VtxIt].x(), TripletVtx.y() - (*vertices)[VtxIt].y(), TripletVtx.z() - (*vertices)[VtxIt].z());
	      double Cosdphi_3D = Dv3D_reco.Dot(ThreeCandidate.Vect())/(Dv3D_reco.Mag()*ThreeCandidate.Vect().Mag());
	      if(Cosdphi_3D>dphi_pv){
		dphi_pv = Cosdphi_3D;
		primaryvertex_index=VtxIt;
	      }
	    }


	    std::vector<reco::TransientTrack> pvTracks_original;
	    TransientTrackMap pvTrackMap_refit;

	    for ( reco::Vertex::trackRef_iterator pvTrack = (*vertices)[primaryvertex_index].tracks_begin(); pvTrack != (*vertices)[primaryvertex_index].tracks_end(); ++pvTrack ) {
	      reco::TransientTrack pvTrack_transient =theTransientTrackBuilder_->build(pvTrack->get());
	      pvTracks_original.push_back(pvTrack_transient);
	      pvTrackMap_refit.insert(std::make_pair(pvTrack->get(), pvTrack_transient));
	    }

	    removeTracks(pvTrackMap_refit,  SVTrackRef);
	    std::vector<reco::TransientTrack> pvTracks_refit;
            for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refit.begin();  pvTrack != pvTrackMap_refit.end(); ++pvTrack ) {
	      pvTracks_refit.push_back(pvTrack->second);}


	    cout<<" Closest PV index "<<primaryvertex_index<<" x="<<(*vertices)[primaryvertex_index].x()<<" y="<<(*vertices)[primaryvertex_index].y()<<" z="<<(*vertices)[primaryvertex_index].z()<<endl;

	    //std::vector<reco::transientTrack> pvtracks_original;
	    //Transienttrackmap Pvtrackmap_Refit;

	    //For ( reco::Vertex::trackRef_iterator pvTrack = (*vertices)[primaryvertex_index].tracks_begin(); pvTrack != (*vertices)[primaryvertex_index].tracks_end(); ++pvTrack ) {
	    // reco::TransientTrack pvTrack_transient =theTransientTrackBuilder_->build(pvTrack->get());
	    //  pvTracks_original.push_back(pvTrack_transient);
	    //  pvTrackMap_refit.insert(std::make_pair(pvTrack->get(), pvTrack_transient));
	    //}
	    RefittedPV_NTracks.push_back(pvTracks_refit.size());   

            if(pvTracks_refit.size() >1){
	      KalmanVertexFitter PV_fitter (true);
	      TransientVertex PVertex = PV_fitter.vertex(pvTracks_refit);
		
	      RefittedPV_isValid.push_back(PVertex.isValid());
                
	      //cout<<"Valid Vtx1="<<PVertex.isValid()<<endl;
	      if(PVertex.isValid()){

		NTripl.push_back(1);
            
		TripletIndex=trIn;
		//    if (!(TauIt->vertexChi2() < 20)) continue ;
		//    TripletsCounter.push_back(1);
            
		//Daughter Kinematics at reco+gen level
            
            Mu1_Pt.push_back(mu1->pt());
            Mu1_Eta.push_back(mu1->eta());
            Mu1_Phi.push_back(mu1->phi());
            Mu1_TripletIndex.push_back(TripletIndex);
            
            Mu2_Pt.push_back(mu2->pt());
            Mu2_Eta.push_back(mu2->eta());
            Mu2_Phi.push_back(mu2->phi());
            Mu2_TripletIndex.push_back(TripletIndex);
            
            Mu3_Pt.push_back(mu3->pt());
            Mu3_Eta.push_back(mu3->eta());
            Mu3_Phi.push_back(mu3->phi());
            Mu3_TripletIndex.push_back(TripletIndex);
            //cout<<"Reco mu1 pt="<<mu1->pt()<<" mu2 pt="<<mu2->pt()<<" mu3 pt="<<mu3->pt()<<endl;
            

            //TransientVertex TransientTripletVtx = reco::Vertex(TauIt->vertex(), TauIt->vertexCovariance(), TauIt->vertexChi2(), TauIt->vertexNdof(), TauIt->numberOfDaughters() );
            //cout<<" number of muons in triplet="<<TauIt->numberOfDaughters()<<endl;


	    ///////////////Check Trigger Matching///////////////
            float dR1 = 999., dR2 = 999., dR3 = 999.;
            std::vector<trigger::TriggerObject> trgobjs = triggerSummary->getObjects();
            trigger::TriggerObjectCollection MuonsObjects;
            edm::InputTag MuonFilterTag = edm::InputTag("hltTau3muTkVertexFilter", "", "HLT");
            size_t MuonFilterIndex = (*triggerSummary).filterIndex(MuonFilterTag); //find the index corresponding to the event
            if(MuonFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
            //save the trigger objetcs corresponding to muon leg
                const trigger::Keys &KEYS = (*triggerSummary).filterKeys(MuonFilterIndex);
                for (unsigned int ipart = 0; ipart < KEYS.size(); ipart++) {
                    trigger::TriggerObject foundObject = (trgobjs)[KEYS[ipart]];
                    MuonsObjects.push_back(foundObject);
                }
                dR1 = MiniAna2017Tree::dRtriggerMatch(*mu1, MuonsObjects);
                dR2 = MiniAna2017Tree::dRtriggerMatch(*mu2, MuonsObjects);
                dR3 = MiniAna2017Tree::dRtriggerMatch(*mu3, MuonsObjects);
            }
            Mu1_dRtriggerMatch.push_back(dR1);
            Mu2_dRtriggerMatch.push_back(dR2);
            Mu3_dRtriggerMatch.push_back(dR3);

            if(isMc){
                bool isMatch1=false; bool isMatch2=false; bool isMatch3=false;
                if( (mu1->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu1->simMotherPdgId()) == 15) ){
                    isMatch1=true;

                }
                if( (mu2->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu2->simMotherPdgId()) == 15) ){
                    isMatch2=true;
                }
                if( (mu3->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu3->simMotherPdgId()) == 15) ){
                    isMatch3=true;
                }
                
                //    cout<<TripletIndex<<"Triplet Mass:"<<TauIt->mass()<<" pt="<<TauIt->pt()<<" vtx.x="<<TauIt->vx()<<" vtx x="<<TripletVtx.x()<<" chi2="<<TauIt->vertexChi2()<<" ndof="<<TauIt->vertexNdof()<<endl;
                // cout<<TripletIndex<<"--Muon 1 pt="<<mu1->pt()<<" Muon2 pt="<<mu2->pt()<<" Mu3 pt="<<mu3->pt()<<" "<<endl;
                if( isMatch1 && isMatch2 && isMatch3) {
		  //cout<<" Matched Triplets mass="<<TauIt->mass()<<endl;
                    GenMatchMu1_SimPt.push_back(mu1->simPt());
                    GenMatchMu2_SimPt.push_back(mu2->simPt());
                    GenMatchMu3_SimPt.push_back(mu3->simPt());
                    
                    GenMatchMu1_SimEta.push_back(mu1->simEta());
                    GenMatchMu2_SimEta.push_back(mu2->simEta());
                    GenMatchMu3_SimEta.push_back(mu3->simEta());
                    
                    GenMatchMu1_SimPhi.push_back(mu1->simPhi());
                    GenMatchMu2_SimPhi.push_back(mu2->simPhi());
                    GenMatchMu3_SimPhi.push_back(mu3->simPhi());
                    
                    GenMatchMu1_Pt.push_back(mu1->pt());
                    GenMatchMu2_Pt.push_back(mu2->pt());
                    GenMatchMu3_Pt.push_back(mu3->pt());
                    
                    GenMatchMu1_Eta.push_back(mu1->eta());
                    GenMatchMu2_Eta.push_back(mu2->eta());
                    GenMatchMu3_Eta.push_back(mu3->eta());
                    
                    GenMatchMu1_Phi.push_back(mu1->phi());
                    GenMatchMu2_Phi.push_back(mu2->phi());
                    GenMatchMu3_Phi.push_back(mu3->phi());
                    
                    
                }
                
                
                //GenVtx vars to be added
            }
            
            //Triplets Vars
            
            
            TripletVtx_x.push_back(TauIt->vx());
            TripletVtx_y.push_back(TauIt->vy());
            TripletVtx_z.push_back(TauIt->vz());
            
            TripletVtx_Chi2.push_back(TauIt->vertexChi2());
            TripletVtx_NDOF.push_back(TauIt->vertexNdof());
            
            Triplet_Mass.push_back(TauIt->mass());
            Triplet_Pt.push_back(TauIt->pt());
            Triplet_Eta.push_back(TauIt->eta());
            Triplet_Phi.push_back(TauIt->phi());
            Triplet_Charge.push_back(TauIt->charge());
            //Matrix covariance to be added!!!!
            
            //Refitted Vars
            //vector < TransientTrack > ttrks = TripletVtx.refittedTracks();
            
            ////


	    /*
	    typedef std::vector<reco::TransientTrack> TTVect;
	    std::vector<TTVect> TTVectContainer;
	    for (uint k=0; k<pvTrackMap_refitVec.size(); k++){
	      removeTracks( pvTrackMap_refitVec.at(k),  SVTrackRef);
	      cout<<" k="<<k<<" size after TrackRemoval="<<pvTrackMap_refitVec.at(k)
	      //std::vector<reco::TransientTrack> pvTracks_refitTmp;
	      //for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refitVec.at(k).begin();  pvTrack != pvTrackMap_refitVec.at(k).end(); ++pvTrack ) {
	      //pvTracks_refitTmp.push_back(pvTrack->second);
	      //}  
	      //cout<<" pvTracks_refitTmp size"<<pvTracks_refitTmp.size()<<endl;
	      //TTVectContainer.push_back(pvTracks_refitTmp);
	    }
	    //cout<<"TTVectContainer size "<<TTVectContainer.size()<<endl;
	    */


	    //KalmanVertexFitter PV_fitter (true);
	    //TransientVertex PVertex = PV_fitter.vertex(pvTracks_refit);





	    /////////////PV Refit//////////////////////////////
	    //cout<<" pvTrackMap_refit size "<<pvTrackMap_refit.size()<<endl;

            /*for(uint i=0; i<pvTracks_refit.size(); i++){
             TrackRef tr = TrackRef(pvTracks_refit, i);
             //reco::Track pvTr=pvTracks_refit.at(i).track();
             //TrackRef pvTrRef = pvTr.get<TrackRef>();
             cout<<i<<"PV track ID="<<tr.id()<<endl;
             }*/
	    /////////////PV Refit//////////////////////////////

            //Defining ISO VAR related to the triplet
            math::XYZPoint SVertexPoint = math::XYZPoint(TripletVtx.x(), TripletVtx.y(), TripletVtx.z());
            TLorentzVector LV1=TLorentzVector( mu1->px(), mu1->py(), mu1->pz(), mu1->energy() );
            TLorentzVector LV2=TLorentzVector( mu2->px(), mu2->py(), mu2->pz(), mu2->energy() );
            TLorentzVector LV3=TLorentzVector( mu3->px(), mu3->py(), mu3->pz(), mu3->energy() );
            TLorentzVector LVTau = LV1 + LV2 + LV3;


            edm::View<reco::Track>::const_iterator trIt  = trackCollection->begin();
            edm::View<reco::Track>::const_iterator trEnd = trackCollection->end();

            int nTracks03_mu1=0, nTracks03_mu2=0, nTracks03_mu3=0;
            double mindist=9999;
            double sumPtTrack1=0, sumPtTrack2=0, sumPtTrack3=0, maxSumPtTracks=0;
            for (; trIt != trEnd; ++trIt) {
               const reco::Track track = (*trIt);
               if(  (track.pt()>1) && (fabs(track.eta())<2.4) && (track.hitPattern().trackerLayersWithMeasurement()>5) && (track.hitPattern().pixelLayersWithMeasurement()>1)  ){
                  double dR1 = dR(Track1.eta(), track.eta(), Track1.phi(), track.phi() );
                  double dR2 = dR(Track2.eta(), track.eta(), Track2.phi(), track.phi() );
                  double dR3 = dR(Track3.eta(), track.eta(), Track3.phi(), track.phi() );
                  //if (dR1 == 0 || dR2 == 0 || dR3 == 0) { 
		  //cout<<"Skip muon track"<<endl; 
		  if (dR1 < 0.01 || dR2 < 0.01 || dR3 < 0.01) { 
		    continue;}
                  double dz = abs(track.dz(SVertexPoint));
                  double dxy = abs(track.dxy(SVertexPoint));
                  double dca_fv = sqrt(dz*dz+dxy*dxy);
                  if(dca_fv<mindist && dca_fv>0) { 
		    
		    mindist = dca_fv;
		    //	    cout<<" MinDist="<<dca_fv<<endl; 
		  }
		  
                  //for eack track having pt>1, excluded the muon tracks,
                  //for each muon in the triplet, if deltaR<0.3 and the DCA is smaller than 1 mm
                  //the pt of the track is added -> I will take the largest total pt from the three muons
                  if (dca_fv < 0.1) {
                     if (dR1<0.3) {
                        sumPtTrack1+=track.pt();
                        nTracks03_mu1++;
                     }
                     if (dR2<0.3) {
                        sumPtTrack2+=track.pt();
                        nTracks03_mu2++;
                     }
                     if (dR3<0.3) {
                        sumPtTrack3+=track.pt();
                        nTracks03_mu3++;
                     }
                  }
               }
            }
            Triplet_mindca_iso.push_back(mindist);
            maxSumPtTracks = std::max(sumPtTrack1, std::max(sumPtTrack2,sumPtTrack3));
	    //cout<<TripletIndex<<" TauMass "<<TauIt->mass()<<" SumPt Tracks in cone="<<maxSumPtTracks<<" TauPt="<<TauIt->pt()<<endl;
            double relativeiso = maxSumPtTracks/LVTau.Pt();
            Triplet_relativeiso.push_back(relativeiso);
            
            Mu1_NTracks03iso.push_back(nTracks03_mu1);
            Mu2_NTracks03iso.push_back(nTracks03_mu2);
            Mu3_NTracks03iso.push_back(nTracks03_mu3);
 
                    //cout<<"Valid Vtx2="<<PVertex.isValid()<<endl;
                    //CachingVertex<5> fittedVertex = vertexFitter.vertex(tracksToVertex);
                    GlobalPoint PVertexPos  (PVertex.position());
                    GlobalPoint SVertexPos  (TripletVtx.x(), TripletVtx.y(), TripletVtx.z());
		    //cout<<" PV Coord after refit="<<PVertexPos.x()<<" y="<<PVertexPos.y()<<" z="<<PVertexPos.z()<<endl;
                    double FlightDist = TMath::Sqrt( pow(( PVertexPos.x() -SVertexPos.x()),2)+ pow(( PVertexPos.y() -SVertexPos.y()),2) + pow(( PVertexPos.z() -SVertexPos.z()),2));
                    
                    VertexDistance3D vertTool;
                    VertexState PVstate(PVertex.position(),PVertex.positionError());
                    //VertexState SVstate(SVertexPos,TripletVtx.position());
                    double distance = vertTool.distance(PVstate, TripletVtx).value();
                    double dist_err = vertTool.distance(PVstate, TripletVtx).error();
                    double dist_sign =vertTool.distance(PVstate, TripletVtx).significance();
                    double chi2 = vertTool.compatibility(PVstate, TripletVtx);
                    
		    //VertexDistance3D dist;
		    //double fv_d3D = dist.distance(Vertex(fv), pvv).value(); // = dv_reco.Mag() ??
		    //double fv_d3Dsig = dist.distance(Vertex(fv), pvv).significance();
		    //double fv_ppdl3D = fv_d3D*fv_cosdphi3D*m3mu_reco/vtau.P();

                    ////

		    PV_x.push_back( (*vertices)[primaryvertex_index].x());
		    PV_y.push_back( (*vertices)[primaryvertex_index].y());
		    PV_z.push_back( (*vertices)[primaryvertex_index].z());
		    PV_NTracks.push_back(pvTracks_original.size());

                    
                    RefittedPV_x.push_back(PVertexPos.x());
                    RefittedPV_y.push_back(PVertexPos.y());
                    RefittedPV_z.push_back(PVertexPos.z());
                    //RefittedPV_NTracks.push_back(pvTracks_refit.size());
                    //RefittedPV_Chi2.push_back(PVertex.);
                    
                    FlightDistPVSV.push_back(distance);
                    FlightDistPVSV_Err.push_back(dist_err);
                    FlightDistPVSV_Significance.push_back(dist_sign);
                    FlightDistPVSV_chi2.push_back(chi2);

		    //Add dxy info
		    GlobalVector dir1(mu1->px(), mu1->py(),  mu1->pz()); //To compute sign of IP
		    GlobalVector dir2(mu2->px(), mu2->py(),  mu2->pz()); //To compute sign of IP
		    GlobalVector dir3(mu3->px(), mu3->py(),  mu3->pz()); //To compute sign of IP
		    std::pair<bool,Measurement1D> signed_IP2D_mu1 = IPTools::signedTransverseImpactParameter(transientTrack1, dir1, PVertex);
		    std::pair<bool,Measurement1D> signed_IP2D_mu2 = IPTools::signedTransverseImpactParameter(transientTrack2, dir2, PVertex);
		    std::pair<bool,Measurement1D> signed_IP2D_mu3 = IPTools::signedTransverseImpactParameter(transientTrack3, dir3, PVertex);
		    //cout<<" IP mu1="<<signed_IP2D_mu1.second.value()<<" IP mu2="<<signed_IP2D_mu2.second.value()<<" IP mu3="<<signed_IP2D_mu3.second.value()<<endl;
		    dxy_mu1.push_back(signed_IP2D_mu1.second.value());
		    dxy_mu2.push_back(signed_IP2D_mu2.second.value());
		    dxy_mu3.push_back(signed_IP2D_mu3.second.value());

		    dxyErr_mu1.push_back(signed_IP2D_mu1.second.error());
		    dxyErr_mu2.push_back(signed_IP2D_mu2.second.error());
		    dxyErr_mu3.push_back(signed_IP2D_mu3.second.error());
		    //std::pair<bool,Measurement1D> signed_IP3D_mu1 = IPTools::signedImpactParameter3D(transientTrack1, dir1, PVertex);
		    //std::pair<bool,Measurement1D> signed_IP3D_mu2 = IPTools::signedImpactParameter3D(transientTrack2, dir2, PVertex);
		    //std::pair<bool,Measurement1D> signed_IP3D_mu3 = IPTools::signedImpactParameter3D(transientTrack3, dir3, PVertex);
		    //ip3d=signed_IP3D.second.value();
		    //ip3d_err=signed_IP3D.second.error();
		    //TransverseImpactPointExtrapolator extrapolator(transTrk.field());
		    //GlobalPoint pos  = extrapolator.extrapolate(transTrk.impactPointState(), RecoVertex::convertPos(PV->position())).globalPosition();
		    //poca=reco::Vertex::Point(pos.x(),pos.y(),pos.z());
		    //AnalyticalImpactPointExtrapolator extrapolator3D(transTrk.field());
		    //GlobalPoint pos3d = extrapolator3D.extrapolate(transTrk.impactPointState(),RecoVertex::convertPos(PV->position())).globalPosition();
		    //ip3d_poca=reco::Vertex::Point(pos3d.x(),pos3d.y(),pos3d.z());
      
	      }else{
		Mu1_Pt.push_back(-99);
		  Mu1_Eta.push_back(-99);
		  Mu1_Phi.push_back(-99);
		  Mu1_TripletIndex.push_back(-99);

		  Mu2_Pt.push_back(-99);
		  Mu2_Eta.push_back(-99);
		  Mu2_Phi.push_back(-99);
		  Mu2_TripletIndex.push_back(-99);

		  Mu3_Pt.push_back(-99);
		  Mu3_Eta.push_back(-99);
		  Mu3_Phi.push_back(-99);
		  Mu3_TripletIndex.push_back(-99);



		  TripletVtx_x.push_back(-99);
		  TripletVtx_y.push_back(-99);
		  TripletVtx_z.push_back(-99);
		    
		  TripletVtx_Chi2.push_back(-99);
		  TripletVtx_NDOF.push_back(-99);

		  Triplet_Mass.push_back(-99);
		  Triplet_Pt.push_back(-99);
		  Triplet_Eta.push_back(-99);
		  Triplet_Phi.push_back(-99);
		  Triplet_Charge.push_back(-99);
		  
		  Triplet_mindca_iso.push_back(-99);
		  Triplet_relativeiso.push_back(-99);

		  Mu1_NTracks03iso.push_back(-99);
		  Mu2_NTracks03iso.push_back(-99);
		  Mu3_NTracks03iso.push_back(-99);
		  RefittedPV_x.push_back(-99);
		  RefittedPV_y.push_back(-99);
		  RefittedPV_z.push_back(-99);
		  RefittedPV_NTracks.push_back(-99);
                    //RefittedPV_Chi2.push_back(PVertex.);
		  Mu1_dRtriggerMatch.push_back(-99);
		  Mu2_dRtriggerMatch.push_back(-99);
		  Mu3_dRtriggerMatch.push_back(-99);

                    FlightDistPVSV.push_back(-99);
                    FlightDistPVSV_Err.push_back(-99);
                    FlightDistPVSV_Significance.push_back(-99);
                    FlightDistPVSV_chi2.push_back(-99);
                    dxy_mu1.push_back(-99);
                    dxy_mu2.push_back(-99);
                    dxy_mu3.push_back(-99);

                    dxyErr_mu1.push_back(-99);
                    dxyErr_mu2.push_back(-99);
                    dxyErr_mu3.push_back(-99);
                  
		    GenMatchMu1_SimPt.push_back(-99);
                    GenMatchMu2_SimPt.push_back(-99);
                    GenMatchMu3_SimPt.push_back(-99);

                    GenMatchMu1_SimEta.push_back(-99);
                    GenMatchMu2_SimEta.push_back(-99);
                    GenMatchMu3_SimEta.push_back(-99);

                    GenMatchMu1_SimPhi.push_back(-99);
                    GenMatchMu2_SimPhi.push_back(-99);
                    GenMatchMu3_SimPhi.push_back(-99);

                    GenMatchMu1_Pt.push_back(-99);
                    GenMatchMu2_Pt.push_back(-99);
                    GenMatchMu3_Pt.push_back(-99);

                    GenMatchMu1_Eta.push_back(-99);
                    GenMatchMu2_Eta.push_back(-99);
                    GenMatchMu3_Eta.push_back(-99);

                    GenMatchMu1_Phi.push_back(-99);
                    GenMatchMu2_Phi.push_back(-99);
                    GenMatchMu3_Phi.push_back(-99);


  
                }
                
            }
            
        }
	}
        NGoodTriplets.push_back(NTripl.size());
	if(NTripl.size()>0) hEventsAfterGoodCand->Fill(1);
	//    cout<<"***Number of Muons="<<muons->size()<<endl; uint k=0;
        
        
        std::vector<int> MuFilter;
        vector<pat::Muon>    MyMu, MyMu2, SyncMu;
        double AllMuPt =0;
        MuonCollectionSize = muons->size();
	uint k=0;        
        for(edm::View<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(), k<muons->size(); ++mu, ++k){
            
            /*

             if (!(mu->track()->quality(reco::TrackBase::highPurity))) continue;
             if ((!(mu->track()->hitPattern().numberOfValidPixelHits()>0))  ) continue;
            
             */
            
            
            //mu.muonBestTrack()->dz(PV.position()), mu.isTightMuon(PV));
            //std::cout<<"RecoMu pt="<<mu->pt()<<" eta="<<mu->eta()<<" phi="<<mu->phi()<<" segm compatibility="<<mu->segmentCompatibility()<<" trkRelChi2="<<mu->combinedQuality().trkRelChi2<<"  vz="<<mu->vz()<<" mu charge="<<mu->charge()<<" mu ECAL energy="<<mu->calEnergy().em<<std::endl;
            //" SymType="<<mu->simExtType()<<" SimPt="<<mu->simPt()<<
            //cout<<" pixel hits="<<mu->track()->hitPattern().numberOfValidPixelHits()<<" pt="<<mu->track()->pt()<<" mu ECAL energy="<<mu->calEnergy().em<<" mu.vx="<<mu->vx()<<endl;
            
            
            MuFilter.push_back(1);
            MyMu.push_back(*mu);
            
            //Basic Kinematics
            MuonPt.push_back(mu->pt());
            MuonEta.push_back(mu->eta());
            MuonPhi.push_back(mu->phi());
            MuonEnergy.push_back(mu->energy());
            MuonCharge.push_back(mu->charge());
            
            Muon_simPdgId.push_back(mu->simPdgId());
            Muon_simMotherPdgId.push_back(mu->simMotherPdgId());
            Muon_simFlavour.push_back(mu->simFlavour());
	    Muon_simType.push_back(mu->simType());
	    Muon_simBX.push_back(mu->simBX());
	    //	    Muon_simTpEvent.push_back(mu->simTpEvent());
	    //	    Muon_simMatchQuality.push_back(mu->simMatchQuality());
            //Vtx position
            Muon_vx.push_back(mu->vx());
            Muon_vy.push_back(mu->vy());
            Muon_vz.push_back(mu->vz());
            
            //cout<<" Muon mother: "<<mu->simMotherPdgId()<<endl;
            
            //MuonID
            Muon_isGlobal.push_back(mu->isGlobalMuon());
            Muon_isSoft.push_back(mu->isSoftMuon(PV));
            Muon_isLoose.push_back(mu->isLooseMuon());
            Muon_isTight.push_back(mu->isTightMuon(PV));
            Muon_isPF.push_back(mu->isPFMuon());
            Muon_isRPCMuon.push_back(mu->isRPCMuon());
            Muon_isStandAloneMuon.push_back(mu->isStandAloneMuon());
            Muon_isTrackerMuon.push_back(mu->isTrackerMuon());
            Muon_isCaloMuon.push_back(mu->isCaloMuon());
            Muon_isQualityValid.push_back(mu->isQualityValid());
            Muon_isTimeValid.push_back(mu->isTimeValid());
            Muon_isIsolationValid.push_back(mu->isIsolationValid());
            Muon_numberOfMatchedStations.push_back(mu->numberOfMatchedStations());
            Muon_numberOfMatches.push_back(mu->numberOfMatches(reco::Muon::SegmentArbitration));
            
            Muon_timeAtIpInOut.push_back(mu->time().timeAtIpInOut);
            Muon_timeAtIpInOutErr.push_back(mu->time().timeAtIpInOutErr);
            
            if (mu->isGlobalMuon()) {
                Muon_GLnormChi2.push_back(mu->globalTrack()->normalizedChi2());
                Muon_GLhitPattern_numberOfValidMuonHits.push_back(mu->globalTrack()->hitPattern().numberOfValidMuonHits());
            }else
            {
                Muon_GLnormChi2.push_back(-999);
                Muon_GLhitPattern_numberOfValidMuonHits.push_back(-999);
            }
            
            if (mu->innerTrack().isNonnull()){
                Muon_trackerLayersWithMeasurement.push_back(mu->innerTrack()->hitPattern().trackerLayersWithMeasurement());
                Muon_Numberofvalidpixelhits.push_back(mu->innerTrack()->hitPattern().numberOfValidPixelHits());
                Muon_innerTrack_p.push_back(mu->innerTrack()->p());
                Muon_innerTrack_eta.push_back(mu->innerTrack()->eta());
                Muon_innerTrack_phi.push_back(mu->innerTrack()->phi());
                Muon_innerTrack_normalizedChi2.push_back(mu->innerTrack()->normalizedChi2());
            }else
            {
                Muon_trackerLayersWithMeasurement.push_back(-999);
                Muon_Numberofvalidpixelhits.push_back(-999);
                Muon_innerTrack_p.push_back(-999);
                Muon_innerTrack_eta.push_back(-999);
                Muon_innerTrack_phi.push_back(-999);
                Muon_innerTrack_normalizedChi2.push_back(-999);
            }
            if (mu->outerTrack().isNonnull()){
                Muon_outerTrack_p.push_back(mu->outerTrack()->p());
                Muon_outerTrack_eta.push_back(mu->outerTrack()->eta());
                Muon_outerTrack_phi.push_back(mu->outerTrack()->phi());
                Muon_outerTrack_normalizedChi2.push_back(mu->outerTrack()->normalizedChi2());
                Muon_outerTrack_muonStationsWithValidHits.push_back(mu->outerTrack()->hitPattern().muonStationsWithValidHits());
            }else{
                Muon_outerTrack_p.push_back(-999);
                Muon_outerTrack_eta.push_back(-999);
                Muon_outerTrack_phi.push_back(-999);
                Muon_outerTrack_normalizedChi2.push_back(-999);
                Muon_outerTrack_muonStationsWithValidHits.push_back(-999);
            }
            if (mu->innerTrack().isNonnull() && mu->outerTrack().isNonnull()){
                Muon_QInnerOuter.push_back(mu->outerTrack()->charge()*mu->innerTrack()->charge());
            }else{
                Muon_QInnerOuter.push_back(-999);
            }
            
            
            Muon_combinedQuality_updatedSta.push_back(mu->combinedQuality().updatedSta);
            Muon_combinedQuality_trkKink.push_back(mu->combinedQuality().trkKink);
            Muon_combinedQuality_glbKink.push_back(mu->combinedQuality().glbKink);
            Muon_combinedQuality_trkRelChi2.push_back(mu->combinedQuality().trkRelChi2);
            Muon_combinedQuality_staRelChi2.push_back(mu->combinedQuality().staRelChi2);
            Muon_combinedQuality_chi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition);
            Muon_combinedQuality_chi2LocalMomentum.push_back(mu->combinedQuality().chi2LocalMomentum);
            Muon_combinedQuality_localDistance.push_back(mu->combinedQuality().localDistance);
            Muon_combinedQuality_globalDeltaEtaPhi.push_back(mu->combinedQuality().globalDeltaEtaPhi);
            Muon_combinedQuality_tightMatch.push_back(mu->combinedQuality().tightMatch);
            Muon_combinedQuality_glbTrackProbability.push_back(mu->combinedQuality().glbTrackProbability);
            
            Muon_calEnergy_em.push_back(mu->calEnergy().em);
            Muon_calEnergy_emS9.push_back(mu->calEnergy().emS9);
            Muon_calEnergy_emS25.push_back(mu->calEnergy().emS25);
            Muon_calEnergy_had.push_back(mu->calEnergy().had);
            Muon_calEnergy_hadS9.push_back(mu->calEnergy().hadS9);
            
            Muon_segmentCompatibility.push_back(muon::segmentCompatibility(*mu));
            Muon_caloCompatibility.push_back(muon::caloCompatibility(*mu));
            
            Muon_ptErrOverPt.push_back( (mu->muonBestTrack()->ptError()/mu->muonBestTrack()->pt()) );
	    Muon_BestTrackPt.push_back( mu->muonBestTrack()->pt() );
	    Muon_BestTrackPtErr.push_back( mu->muonBestTrack()->ptError() );

	    Muon_BestTrackEta.push_back( mu->muonBestTrack()->eta() );
	    Muon_BestTrackEtaErr.push_back( mu->muonBestTrack()->etaError() );

	    Muon_BestTrackPhi.push_back( mu->muonBestTrack()->phi() );
	    Muon_BestTrackPhiErr.push_back( mu->muonBestTrack()->phiError() );
            
            const reco::MuonIsolation Iso03 = mu->isolationR03();
            const reco::MuonIsolation Iso05 = mu->isolationR05();
            if (mu->isIsolationValid()) {
                Muon_emEt03.push_back(Iso03.emEt);
                Muon_hadEt03.push_back(Iso03.hadEt);
                Muon_nJets03.push_back(Iso03.nJets);
                Muon_nTracks03.push_back(Iso03.nTracks);
                Muon_sumPt03.push_back(Iso03.sumPt);
                Muon_emVetoEt03.push_back(Iso03.emVetoEt);
                Muon_hadVetoEt03.push_back(Iso03.hadVetoEt);
                Muon_trackerVetoPt03.push_back(Iso03.trackerVetoPt);
                
                Muon_emEt05.push_back(Iso05.emEt);
                Muon_hadEt05.push_back(Iso05.hadEt);
                Muon_nJets05.push_back(Iso05.nJets);
                Muon_nTracks05.push_back(Iso05.nTracks);
                Muon_sumPt05.push_back(Iso05.sumPt);
                Muon_hadVetoEt05.push_back(Iso05.hadVetoEt);
                Muon_trackerVetoPt05.push_back(Iso05.trackerVetoPt);
                Muon_emVetoEt05.push_back(Iso05.emVetoEt);
                
            } else { // if isolation is not valid use -1 as default
                Muon_emEt03.push_back(-1);
                Muon_hadEt03.push_back(-1);
                Muon_nJets03.push_back(-1);
                Muon_nTracks03.push_back(-1);
                Muon_sumPt03.push_back(-1);
                
                Muon_hadVetoEt03.push_back(-1);
                Muon_emVetoEt03.push_back(-1);
                Muon_trackerVetoPt03.push_back(-1);
                
                Muon_emEt05.push_back(-1);
                Muon_hadEt05.push_back(-1);
                Muon_nJets05.push_back(-1);
                Muon_nTracks05.push_back(-1);
                Muon_sumPt05.push_back(-1);
                
                Muon_emVetoEt05.push_back(-1);
                Muon_trackerVetoPt05.push_back(-1);
                Muon_hadVetoEt05.push_back(-1);
            }
            
        }
        
        if (!iEvent.isRealData())
        {
            Handle<vector<PileupSummaryInfo> >  PupInfo;
            iEvent.getByToken(puToken_, PupInfo);
            puN = PupInfo->begin()->getTrueNumInteractions();
        }
        
        
        
        ////Synch Tree//////
        /*   double maxPt =0; double maxPhi=0, maxEta=0; vector<pat::Muon> SyncSortedMu ;
         double maxTrPt =0; double maxTrPhi=0, maxTrEta=0; vector<reco::Track> SyncSortedTr ;
         for(uint i=0; i<SyncMu.size();i++){
         if(SyncMu.at(i).pt() > maxPt){
         maxPt  = SyncMu.at(i).pt();
         maxPhi = SyncMu.at(i).phi();
         maxEta = SyncMu.at(i).eta();
         SyncSortedMu.push_back(SyncMu.at(i));
         }
         }
         
         allmuons_pt.push_back(AllMuPt);
         
         leadmuon_pt.push_back(maxPt);
         leadmuon_eta.push_back(maxEta);
         leadmuon_phi.push_back(maxPhi);
         nmuons = SyncMu.size();
         nprimevtxs =vertices->size();
         
         edm::View<reco::Track>::const_iterator trIt  = trackCollection->begin();
         edm::View<reco::Track>::const_iterator trEnd = trackCollection->end();
         
         double AllTrPt=0;
         for (; trIt != trEnd; ++trIt)
         {
         
         const reco::Track track = (*trIt);
         if(  (track.pt()>1) && (fabs(track.eta())<2.4) && (track.hitPattern().trackerLayersWithMeasurement()>5) && (track.hitPattern().pixelLayersWithMeasurement()>1)  ){
         AllTrPt +=trIt->pt();
         if(track.pt() > maxTrPt){
         maxTrPt = track.pt();
         maxTrEta = track.eta();
         maxTrPhi= track.phi();
         }
         }
         }
         
         alltracks_pt.push_back(AllTrPt);
         leadtrack_pt.push_back(maxTrPt);
         leadtrack_eta.push_back(maxTrEta);
         leadtrack_phi.push_back(maxTrPhi);
         
         evt   = iEvent.id().event();
         run = iEvent.id().run();
         lumi = iEvent.luminosityBlock();
         */
        ///////SyncTree
        
        
        
        evt   = iEvent.id().event();
        run = iEvent.id().run();
        lumi = iEvent.luminosityBlock();
        
        //  SyncTree_->Fill();
        tree_->Fill();
        
        /*
         allmuons_pt.clear();
         alltracks_pt.clear();
         leadmuon_pt.clear();
         leadmuon_phi.clear();
         leadmuon_eta.clear();
         allmuons_pt.clear();
         leadtrack_pt.clear();
         leadtrack_eta.clear();
         leadtrack_phi.clear();
         evt= -999;
         run= -999;
         lumi= -999;
         nmuons = -999;
         */
        run= -999;
        evt= -999;
        lumi= -999;
	puN= -999;
 
        
        GenParticle_PdgId.clear();
        GenParticle_Pt.clear();
        GenParticle_Eta.clear();
        GenParticle_Phi.clear();
        GenParticle_MotherPdgId.clear();
        
        MuonCollectionSize =0;
        MuonPt.clear();
        MuonEta.clear();
        MuonPhi.clear();
        
        Muon_simPdgId.clear();
        Muon_simMotherPdgId.clear();
        Muon_simFlavour.clear();
	Muon_simType.clear();
	Muon_simBX.clear();
	//	Muon_simTpEvent.clear();
	//	Muon_simMatchQuality.clear();
        MuonEnergy.clear();
        MuonCharge.clear();
        
        //Vtx position
        Muon_vx.clear();
        Muon_vy.clear();
        Muon_vz.clear();
        
        //MuonID
        Muon_isGlobal.clear();
        Muon_isSoft.clear();
        Muon_isLoose.clear();
        Muon_isTight.clear();
        Muon_isPF.clear();
        Muon_isRPCMuon.clear();
        Muon_isStandAloneMuon.clear();
        Muon_isTrackerMuon.clear();
        Muon_isCaloMuon.clear();
        Muon_isQualityValid.clear();
        Muon_isTimeValid.clear();
        Muon_isIsolationValid.clear();
        Muon_numberOfMatchedStations.clear();
        Muon_numberOfMatches.clear();
        
        //MuonTime
        Muon_timeAtIpInOut.clear();
        Muon_timeAtIpInOutErr.clear();
        
        //Muon inner + outer track
        Muon_GLnormChi2.clear();
        Muon_GLhitPattern_numberOfValidMuonHits.clear();
        
        Muon_trackerLayersWithMeasurement.clear();
        Muon_Numberofvalidpixelhits.clear();
        
        Muon_outerTrack_p.clear();
        Muon_outerTrack_eta.clear();
        Muon_outerTrack_phi.clear();
        Muon_outerTrack_normalizedChi2.clear();
        Muon_outerTrack_muonStationsWithValidHits.clear();
        Muon_innerTrack_p.clear();
        Muon_innerTrack_eta.clear();
        Muon_innerTrack_phi.clear();
        Muon_innerTrack_normalizedChi2.clear();
        Muon_QInnerOuter.clear();
        
        
        Muon_combinedQuality_updatedSta.clear();
        Muon_combinedQuality_trkKink.clear();
        Muon_combinedQuality_glbKink.clear();
        Muon_combinedQuality_trkRelChi2.clear();
        Muon_combinedQuality_staRelChi2.clear();
        Muon_combinedQuality_chi2LocalPosition.clear();
        Muon_combinedQuality_chi2LocalMomentum.clear();
        Muon_combinedQuality_localDistance.clear();
        Muon_combinedQuality_globalDeltaEtaPhi.clear();
        Muon_combinedQuality_tightMatch.clear();
        Muon_combinedQuality_glbTrackProbability.clear();
        
        Muon_calEnergy_em.clear();
        Muon_calEnergy_emS9.clear();
        Muon_calEnergy_emS25.clear();
        Muon_calEnergy_had.clear();
        Muon_calEnergy_hadS9.clear();
        
        Muon_segmentCompatibility.clear();
        Muon_caloCompatibility.clear();
        
        Muon_ptErrOverPt.clear();
    	Muon_BestTrackPt.clear();
	Muon_BestTrackPtErr.clear();

	Muon_BestTrackEta.clear();
	Muon_BestTrackEtaErr.clear();

	Muon_BestTrackPhi.clear();
	Muon_BestTrackPhiErr.clear();


        Muon_emEt03.clear();
        Muon_hadEt03.clear();
        Muon_nJets03.clear();
        Muon_nTracks03.clear();
        Muon_sumPt03.clear();
        Muon_hadVetoEt03.clear();
        Muon_emVetoEt03.clear();
        Muon_trackerVetoPt03.clear();
        
        Muon_emEt05.clear();
        Muon_hadEt05.clear();
        Muon_nJets05.clear();
        Muon_nTracks05.clear();
        Muon_sumPt05.clear();
        Muon_hadVetoEt05.clear();
        Muon_emVetoEt05.clear();
        Muon_trackerVetoPt05.clear();
        
        Triplet_mindca_iso.clear();
        Triplet_relativeiso.clear();

        PV_x.clear();
        PV_y.clear();
        PV_z.clear();
        PV_NTracks.clear();
        

        Mu1_Pt.clear();
        Mu1_Eta.clear();
        Mu1_Phi.clear();
        Mu1_NTracks03iso.clear();
        Mu1_dRtriggerMatch.clear();
        
        Mu2_Pt.clear();
        Mu2_Eta.clear();
        Mu2_Phi.clear();
        Mu2_NTracks03iso.clear();
        Mu2_dRtriggerMatch.clear();
        
        Mu3_Pt.clear();
        Mu3_Eta.clear();
        Mu3_Phi.clear();
        Mu3_NTracks03iso.clear();
        Mu3_dRtriggerMatch.clear();
        
        Mu1_TripletIndex.clear();
        Mu2_TripletIndex.clear();
        Mu3_TripletIndex.clear();
        /*
         Mu1_SimPt.clear();
         Mu1_SimEta.clear();
         Mu1_SimPhi.clear();
         
         Mu2_SimPt.clear();
         Mu2_SimEta.clear();
         Mu2_SimPhi.clear();
         
         Mu3_SimPt.clear();
         Mu3_SimEta.clear();
         Mu3_SimPhi.clear();
         */
        GenMatchMu1_SimPhi.clear();
        GenMatchMu2_SimPhi.clear();
        GenMatchMu3_SimPhi.clear();
        
        GenMatchMu1_SimPt.clear();
        GenMatchMu2_SimPt.clear();
        GenMatchMu3_SimPt.clear();
        
        GenMatchMu1_SimEta.clear();
        GenMatchMu2_SimEta.clear();
        GenMatchMu3_SimEta.clear();
        
        GenMatchMu1_Pt.clear();
        GenMatchMu2_Pt.clear();
        GenMatchMu3_Pt.clear();
        
        GenMatchMu1_Eta.clear();
        GenMatchMu2_Eta.clear();
        GenMatchMu3_Eta.clear();
        
        GenMatchMu1_Phi.clear();
        GenMatchMu2_Phi.clear();
        GenMatchMu3_Phi.clear();
        
        TripletCollectionSize = -99;
        
        TripletVtx_x.clear();
        TripletVtx_y.clear();
        TripletVtx_z.clear();
        
        TripletVtx_Chi2.clear();
        TripletVtx_NDOF.clear();
        
        Triplet_Mass.clear();
        Triplet_Pt.clear();
        Triplet_Eta.clear();
        Triplet_Phi.clear();
        Triplet_Charge.clear();
        
	dxy_mu1.clear();
	dxy_mu2.clear();
	dxy_mu3.clear();
	dxyErr_mu1.clear();
	dxyErr_mu2.clear();
	dxyErr_mu3.clear();

        RefittedPV_x.clear();
        RefittedPV_y.clear();
        RefittedPV_z.clear();
        RefittedPV_NTracks.clear();
        RefittedPV_isValid.clear();
        //RefittedPV_Chi2.push_back(PVertex.);
        
        FlightDistPVSV.clear();
        FlightDistPVSV_Err.clear();
        FlightDistPVSV_Significance.clear();
        FlightDistPVSV_chi2.clear();
        PVCollection_Size =0;

        Trigger_l1name.clear();
        Trigger_l1decision.clear();
        Trigger_l1prescale.clear();

	Trigger_hltname.clear();
	Trigger_hltdecision.clear();
	NGoodTriplets.clear();



    }
    
    
    // ------------ method called once each job just before starting event loop  ------------
    void
    MiniAna2017Tree::beginJob()
    {
        
      hEvents = fs->make<TH1F>("hEvents","hEvents",10,0,10);
      hEventsAfterGoodCand = fs->make<TH1F>("hEventsAfterGoodCand","hEventsAfterGoodCand",10,0,10);
        
        tree_ = fs->make<TTree>("ntuple","LFVTau ntuple");
        tree_->Branch("evt", &evt);
        tree_->Branch("run", &run);
        tree_->Branch("lumi", &lumi);
	tree_->Branch("NGoodTriplets", &NGoodTriplets);
        tree_->Branch("nPileUpInt", &puN);
        
	tree_->Branch("Trigger_l1name", &Trigger_l1name);
        tree_->Branch("Trigger_l1decision",&Trigger_l1decision);
	tree_->Branch("Trigger_l1prescale",&Trigger_l1prescale);

	tree_->Branch("Trigger_hltname",&Trigger_hltname);
	tree_->Branch("Trigger_hltdecision",&Trigger_hltdecision);

        tree_->Branch("GenParticle_PdgId", &GenParticle_PdgId);
        tree_->Branch("GenParticle_Pt", &GenParticle_Pt);
        tree_->Branch("GenParticle_Eta", &GenParticle_Eta);
        tree_->Branch("GenParticle_Phi", &GenParticle_Phi);
        tree_->Branch("GenParticle_MotherPdgId", &GenParticle_MotherPdgId);
        
        tree_->Branch("MuonCollectionSize",&MuonCollectionSize);
        tree_->Branch("MuonPt",&MuonPt);
        tree_->Branch("MuonEnergy", &MuonEnergy);
        tree_->Branch("MuonCharge", &MuonCharge);
        tree_->Branch("MuonEta",&MuonEta);
        tree_->Branch("MuonPhi",&MuonPhi);
        tree_->Branch("Muon_simPdgId", &Muon_simPdgId);
        tree_->Branch("Muon_simMotherPdgId", &Muon_simMotherPdgId);
        tree_->Branch("Muon_simFlavour", &Muon_simFlavour);
	tree_->Branch("Muon_simType", &Muon_simType);
	tree_->Branch("Muon_simBX", &Muon_simBX);
	//	tree_->Branch("Muon_simTpEvent", &Muon_simTpEvent);
	//	tree_->Branch("Muon_simMatchQuality", &Muon_simMatchQuality);

        
        //Vtx position
        tree_->Branch("Muon_vx", &Muon_vx);
        tree_->Branch("Muon_vy", &Muon_vy);
        tree_->Branch("Muon_vz", &Muon_vz);
        
        //MuonID
        tree_->Branch("Muon_isGlobal", &Muon_isGlobal);
        //tree_->Branch("Muon_isTracker", &Muon_isTracker);
        tree_->Branch("Muon_isSoft", &Muon_isSoft);
        tree_->Branch("Muon_isLoose", &Muon_isLoose);
        tree_->Branch("Muon_isTight", &Muon_isTight);
        tree_->Branch("Muon_isPF", &Muon_isPF);
        tree_->Branch("Muon_isRPCMuon", &Muon_isRPCMuon);
        tree_->Branch("Muon_isStandAloneMuon", &Muon_isStandAloneMuon);
        tree_->Branch("Muon_isTrackerMuon", &Muon_isTrackerMuon);
        tree_->Branch("Muon_isCaloMuon", &Muon_isCaloMuon);
        tree_->Branch("Muon_isQualityValid", &Muon_isQualityValid);
        tree_->Branch("Muon_isTimeValid", &Muon_isTimeValid);
        tree_->Branch("Muon_isIsolationValid", &Muon_isIsolationValid);
        tree_->Branch("Muon_numberOfMatchedStations", &Muon_numberOfMatchedStations);
        tree_->Branch("Muon_numberOfMatches", &Muon_numberOfMatches);
        
        
        tree_->Branch("Muon_timeAtIpInOut",&Muon_timeAtIpInOut);
        tree_->Branch("Muon_timeAtIpInOutErr",&Muon_timeAtIpInOutErr);
        //Muon inner + outer track                                                                                                                                                                   
        tree_->Branch("Muon_GLnormChi2", &Muon_GLnormChi2);
        tree_->Branch("Muon_GLhitPattern_numberOfValidMuonHits", &Muon_GLhitPattern_numberOfValidMuonHits);
        
        tree_->Branch("Muon_trackerLayersWithMeasurement", &Muon_trackerLayersWithMeasurement);
        tree_->Branch("Muon_Numberofvalidpixelhits", &Muon_Numberofvalidpixelhits);
        
        tree_->Branch("Muon_outerTrack_p", &Muon_outerTrack_p);
        tree_->Branch("Muon_outerTrack_eta", &Muon_outerTrack_eta);
        tree_->Branch("Muon_outerTrack_phi", &Muon_outerTrack_phi);
        tree_->Branch("Muon_outerTrack_normalizedChi2", &Muon_outerTrack_normalizedChi2);
        tree_->Branch("Muon_outerTrack_muonStationsWithValidHits", &Muon_outerTrack_muonStationsWithValidHits);
        tree_->Branch("Muon_innerTrack_p", &Muon_innerTrack_p);
        tree_->Branch("Muon_innerTrack_eta", &Muon_innerTrack_eta);
        
        tree_->Branch("Muon_innerTrack_phi", &Muon_innerTrack_phi);
        tree_->Branch("Muon_innerTrack_normalizedChi2", &Muon_innerTrack_normalizedChi2);
        tree_->Branch("Muon_QInnerOuter", &Muon_QInnerOuter);
        
        
        tree_->Branch("Muon_combinedQuality_updatedSta", &Muon_combinedQuality_updatedSta);
        tree_->Branch("Muon_combinedQuality_trkKink", &Muon_combinedQuality_trkKink);
        tree_->Branch("Muon_combinedQuality_glbKink", &Muon_combinedQuality_glbKink);
        tree_->Branch("Muon_combinedQuality_trkRelChi2", &Muon_combinedQuality_trkRelChi2);
        tree_->Branch("Muon_combinedQuality_staRelChi2", &Muon_combinedQuality_staRelChi2);
        tree_->Branch("Muon_combinedQuality_chi2LocalPosition", &Muon_combinedQuality_chi2LocalPosition);
        tree_->Branch("Muon_combinedQuality_chi2LocalMomentum", &Muon_combinedQuality_chi2LocalMomentum);
        tree_->Branch("Muon_combinedQuality_localDistance", &Muon_combinedQuality_localDistance);
        tree_->Branch("Muon_combinedQuality_globalDeltaEtaPhi", &Muon_combinedQuality_globalDeltaEtaPhi);
        tree_->Branch("Muon_combinedQuality_tightMatch", &Muon_combinedQuality_tightMatch); 
        tree_->Branch("Muon_combinedQuality_glbTrackProbability", &Muon_combinedQuality_glbTrackProbability);
        
        tree_->Branch("Muon_calEnergy_em", &Muon_calEnergy_em);
        tree_->Branch("Muon_calEnergy_emS9", &Muon_calEnergy_emS9);
        tree_->Branch("Muon_calEnergy_emS25", &Muon_calEnergy_emS25);
        tree_->Branch("Muon_calEnergy_had", &Muon_calEnergy_had);
        tree_->Branch("Muon_calEnergy_hadS9", &Muon_calEnergy_hadS9);
        
        tree_->Branch("Muon_segmentCompatibility", &Muon_segmentCompatibility);
        tree_->Branch("Muon_caloCompatibility", &Muon_caloCompatibility);
        
        tree_->Branch("Muon_ptErrOverPt", &Muon_ptErrOverPt);
	tree_->Branch("Muon_BestTrackPt", &Muon_BestTrackPt);
        tree_->Branch("Muon_BestTrackPtErr", &Muon_BestTrackPtErr);
        tree_->Branch("Muon_BestTrackEta", &Muon_BestTrackEta);
	tree_->Branch("Muon_BestTrackEtaErr", &Muon_BestTrackEtaErr);
	tree_->Branch("Muon_BestTrackPhi", &Muon_BestTrackPhi);
	tree_->Branch("Muon_BestTrackPhiErr", &Muon_BestTrackPhiErr);

        tree_->Branch("Muon_emEt03", &Muon_emEt03);
        tree_->Branch("Muon_hadEt03", &Muon_hadEt03);
        tree_->Branch("Muon_nJets03", &Muon_nJets03);
        tree_->Branch("Muon_nTracks03", &Muon_nTracks03);
        tree_->Branch("Muon_sumPt03", &Muon_sumPt03);
        tree_->Branch("Muon_hadVetoEt03", &Muon_hadVetoEt03);
        tree_->Branch("Muon_emVetoEt03", &Muon_emVetoEt03);
        tree_->Branch("Muon_trackerVetoPt03", &Muon_trackerVetoPt03);

        tree_->Branch("Muon_emEt05", &Muon_emEt05);
        tree_->Branch("Muon_hadEt05", &Muon_hadEt05);
        tree_->Branch("Muon_nJets05", &Muon_nJets05);
        tree_->Branch("Muon_nTracks05", &Muon_nTracks05);
        tree_->Branch("Muon_sumPt05", &Muon_sumPt05);
        tree_->Branch("Muon_hadVetoEt05", &Muon_hadVetoEt05);
        tree_->Branch("Muon_emVetoEt05", &Muon_emVetoEt05);
        tree_->Branch("Muon_trackerVetoPt05", &Muon_trackerVetoPt05);

        tree_->Branch("PVCollection_Size", &PVCollection_Size);
        tree_->Branch("PV_x", &PV_x);
        tree_->Branch("PV_y", &PV_y);
        tree_->Branch("PV_z", &PV_z);
        tree_->Branch("PV_NTracks", &PV_NTracks);
        
        tree_->Branch("TripletCollectionSize", &TripletCollectionSize);  
        tree_->Branch("Mu1_Pt",&Mu1_Pt);
        tree_->Branch("Mu1_Eta", &Mu1_Eta);
        tree_->Branch("Mu1_Phi", &Mu1_Phi);
        tree_->Branch("Mu1_NTracks03iso", &Mu1_NTracks03iso); 
        tree_->Branch("Mu1_dRtriggerMatch", &Mu1_dRtriggerMatch); 
        tree_->Branch("Mu1_TripletIndex", &Mu1_TripletIndex);
 
        tree_->Branch("Mu2_Pt", &Mu2_Pt);
        tree_->Branch("Mu2_Eta", &Mu2_Eta);
        tree_->Branch("Mu2_Phi", &Mu2_Phi);
        tree_->Branch("Mu2_NTracks03iso", &Mu2_NTracks03iso); 
        tree_->Branch("Mu2_dRtriggerMatch", &Mu2_dRtriggerMatch); 
        tree_->Branch("Mu2_TripletIndex", &Mu2_TripletIndex); 
        
        tree_->Branch("Mu3_Pt", &Mu3_Pt);
        tree_->Branch("Mu3_Eta",&Mu3_Eta);
        tree_->Branch("Mu3_Phi", &Mu3_Phi);
        tree_->Branch("Mu3_NTracks03iso", &Mu3_NTracks03iso); 
        tree_->Branch("Mu3_dRtriggerMatch", &Mu3_dRtriggerMatch); 
        tree_->Branch("Mu3_TripletIndex", &Mu3_TripletIndex); 

	tree_->Branch("dxy_mu1", &dxy_mu1);
	tree_->Branch("dxy_mu2", &dxy_mu2);
	tree_->Branch("dxy_mu3", &dxy_mu3);
	tree_->Branch("dxyErr_mu1", &dxyErr_mu1);
	tree_->Branch("dxyErr_mu2", &dxyErr_mu2);
	tree_->Branch("dxyErr_mu3", &dxyErr_mu3);


        /*  tree_->Branch("Mu1_SimPt", &Mu1_SimPt);
         tree_->Branch("Mu1_SimEta", &Mu1_SimEta);
         tree_->Branch("Mu1_SimPhi", &Mu1_SimPhi);
         tree_->Branch("Mu2_SimPt", &Mu2_SimPt);
         tree_->Branch("Mu2_SimEta", &Mu2_SimEta);
         tree_->Branch("Mu2_SimPhi", &Mu2_SimPhi);
         tree_->Branch("Mu3_SimPt", &Mu3_SimPt);
         tree_->Branch("Mu3_SimEta", &Mu3_SimEta);
         tree_->Branch("Mu3_SimPhi", &Mu3_SimPhi);
         */
        tree_->Branch("GenMatchMu1_SimPt", &GenMatchMu1_SimPt);
        tree_->Branch("GenMatchMu2_SimPt", &GenMatchMu2_SimPt);
        tree_->Branch("GenMatchMu3_SimPt", &GenMatchMu3_SimPt);
        
        tree_->Branch("GenMatchMu1_SimEta", &GenMatchMu1_SimEta);
        tree_->Branch("GenMatchMu2_SimEta", &GenMatchMu2_SimEta);
        tree_->Branch("GenMatchMu3_SimEta", &GenMatchMu3_SimEta);
        
        tree_->Branch("GenMatchMu1_SimPhi", &GenMatchMu1_SimPhi);
        tree_->Branch("GenMatchMu2_SimPhi", &GenMatchMu2_SimPhi);
        tree_->Branch("GenMatchMu3_SimPhi", &GenMatchMu3_SimPhi);
        
        tree_->Branch("GenMatchMu1_Pt", &GenMatchMu1_Pt);
        tree_->Branch("GenMatchMu2_Pt", &GenMatchMu2_Pt);
        tree_->Branch("GenMatchMu3_Pt", &GenMatchMu3_Pt);
        
        tree_->Branch("GenMatchMu1_Eta", &GenMatchMu1_Eta);
        tree_->Branch("GenMatchMu2_Eta", &GenMatchMu2_Eta);
        tree_->Branch("GenMatchMu3_Eta", &GenMatchMu3_Eta);
        
        tree_->Branch("GenMatchMu1_Phi", &GenMatchMu1_Phi);
        tree_->Branch("GenMatchMu2_Phi", &GenMatchMu2_Phi);
        tree_->Branch("GenMatchMu3_Phi", &GenMatchMu3_Phi);
        
        tree_->Branch("TripletVtx_x", &TripletVtx_x);
        tree_->Branch("TripletVtx_y", &TripletVtx_y);
        tree_->Branch("TripletVtx_z", &TripletVtx_z);
        
        tree_->Branch("TripletVtx_Chi2", &TripletVtx_Chi2);
        tree_->Branch("TripletVtx_NDOF", &TripletVtx_NDOF);
        
        tree_->Branch("Triplet_Mass", &Triplet_Mass);
        tree_->Branch("Triplet_Pt", &Triplet_Pt);
        tree_->Branch("Triplet_Eta", &Triplet_Eta);
        tree_->Branch("Triplet_Phi", &Triplet_Phi);
        tree_->Branch("Triplet_Charge", &Triplet_Charge);
        
        tree_->Branch("Triplet_mindca_iso", &Triplet_mindca_iso);
        tree_->Branch("Triplet_relativeiso", &Triplet_relativeiso);
        
        tree_->Branch("RefittedPV_x", &RefittedPV_x);
        tree_->Branch("RefittedPV_y", &RefittedPV_y);
        tree_->Branch("RefittedPV_z", &RefittedPV_z);
        tree_->Branch("RefittedPV_NTracks", &RefittedPV_NTracks);
        tree_->Branch("RefittedPV_isValid", &RefittedPV_isValid);
        //RefittedPV_Chi2.push_back(PVertex.);                                                                                                                                                     
        
        tree_->Branch("FlightDistPVSV", &FlightDistPVSV);
        tree_->Branch("FlightDistPVSV_Err", &FlightDistPVSV_Err);
        tree_->Branch("FlightDistPVSV_Significance", &FlightDistPVSV_Significance);
        tree_->Branch("FlightDistPVSV_chi2", &FlightDistPVSV_chi2);
        
        
        
        /*  SyncTree_ = fs->make<TTree>("t","Sync ntuple");
         SyncTree_ ->Branch("allmuons_pt",&allmuons_pt);
         SyncTree_->Branch("leadmuon_pt",&leadmuon_pt);
         SyncTree_->Branch("leadmuon_phi",&leadmuon_phi);
         SyncTree_->Branch("leadmuon_eta",&leadmuon_eta);
         SyncTree_->Branch("nmuons",&nmuons);
         SyncTree_->Branch("nprimevtxs",&nprimevtxs); 
         
         SyncTree_->Branch("leadtrack_pt", &leadtrack_pt);
         SyncTree_->Branch("leadtrack_eta", &leadtrack_eta);
         SyncTree_->Branch("leadtrack_phi", &leadtrack_phi);
         SyncTree_->Branch("alltracks_pt", &alltracks_pt);
         SyncTree_->Branch("evt", &evt);
         SyncTree_->Branch("run", &run);
         SyncTree_->Branch("lumi", &lumi);
         */
        
        
        
    }
    
    
    // ------------ method called once each job just after ending the event loop  ------------
    void 
    MiniAna2017Tree::endJob() 
    {
        //  tree_->GetDirectory()->cd();
        tree_->Write();
        
        //  SyncTree_->GetDirectory()->cd();
        //  SyncTree_->Write();
        
    }
    
    // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
    
    
    void MiniAna2017Tree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
        //The following says we do not know what parameters are allowed so do no validation
        // Please change this to state exactly what you do use, even if it is no parameters
        edm::ParameterSetDescription desc;
        desc.setUnknown();
        descriptions.addDefault(desc);
    }
    //define this as a plug-in
    DEFINE_FWK_MODULE(MiniAna2017Tree);
