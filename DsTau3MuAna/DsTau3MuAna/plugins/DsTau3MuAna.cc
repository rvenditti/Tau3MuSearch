// -*- C++ -*-
//
// Package:    DsTau3MuAna/DsTau3MuAna
// Class:      DsTau3MuAna
// 
/**\class DsTau3MuAna DsTau3MuAna.cc DsTau3MuAna/DsTau3MuAna/plugins/DsTau3MuAna.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  venditti
//         Created:  Wed, 30 Jan 2019 10:10:11 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DsTau3MuAna : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DsTau3MuAna(const edm::ParameterSet&);
      ~DsTau3MuAna();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  edm::EDGetTokenT<edm::View<pat::Muon> > muons_;
  edm::EDGetTokenT<edm::View<reco::Vertex> > vertex_;
  edm::EDGetTokenT<edm::View<reco::CompositeCandidate> > Cand3Mu_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticles_;
  //typedef std::vector<SimVertex> SimVertexContainer;
  //edm::EDGetTokenT<SimVertexContainer> simvertsToken_;

  const TransientTrackBuilder* theTransientTrackBuilder_;
  edm::Service<TFileService> fs;
  TH1F *hEvt;
  TH1F *hSize3MuCand, *hValidKinFitVtx;
  TH1F *hVtxMass_KinVsCand, *hVtxChi2_Kin,  *hVtxNDof_Kin,  *hProb_Kin, *hVtxMass_KalVsGen;
  TH1F *hDx, *hDy, *hDz; 

  TH1F *  hVtxMass_KalVsCand;
  TH1F *  hVtxChi2_Kal;
  TH1F *  hVtxNDof_Kal;
  TH1F *  hProb_Kal, *hValidVtx_size_Kalman;
 TH1F *   hDxGenReco;
 TH1F *   hDyGenReco;
 TH1F *   hDzGenReco;

 TH1F *   hDxGenRefit;
 TH1F *   hDyGenRefit;
 TH1F *   hDzGenRefit;

 TH1F *   hVtxMass_KinVsGen;
 TH1F *   hVtxMass_RecoVsGen,  *hDPx_GenKal, *hDPy_GenKal, *hDPz_GenKal,  *hDPx_GenReco, *hDPx_GenKin,  *hDPy_GenKin, *hDPy_GenReco, *hDPz_GenKin, *hDPz_GenReco;



  TH1F *  hDist_PVSV_Kin, *  hDist_PVSV_Kal;
  TH1F *  hDistErr_PVSV_Kin, * hDistErr_PVSV_Kal;
  TH1F *  hDistComp_PVSV_Kin, *  hDistComp_PVSV_Kal, *hDistSign;

TH1F *   hGenSV_KinSV ;
TH1F *   hDeltaGenSVKinSV;

TH1F *   hGenSV_KalSV;
TH1F *   hDeltaGenSVKalSV;

TH1F *   hDxGenKal;
TH1F *   hDyGenKal;
TH1F *   hDzGenKal;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DsTau3MuAna::DsTau3MuAna(const edm::ParameterSet& iConfig)

{
  muons_ = consumes<edm::View<pat::Muon> >  (iConfig.getParameter<edm::InputTag>("muonLabel"));
  vertex_ = consumes<edm::View<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("VertexLabel"));
  Cand3Mu_ = consumes<edm::View<reco::CompositeCandidate> > (iConfig.getParameter<edm::InputTag>("Cand3MuLabel"));
  genParticles_ = consumes<edm::View<reco::GenParticle>  > (iConfig.getParameter<edm::InputTag>("genParticleLabel"));
  //simvertsToken_ = consumes<SimVertexContainer> (edm::InputTag("SimVertexs"));
  usesResource("TFileService");

}


DsTau3MuAna::~DsTau3MuAna()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
std::vector<float> RefitMu(auto kinPart){
  std::vector<float> Mu;
  for(int i=0; i<7; i++){
    float ref=kinPart->currentState().kinematicParameters().vector().At(i);
    Mu.push_back(ref);
  }
  return Mu;
}

std::vector<float> KalmanRefMu(auto TrTack){
  std::vector<float> MyKalRefit;
  double  muon_mass = 0.10565837; //pdg mass 
    MyKalRefit.push_back(TrTack.impactPointState().globalPosition().x()); //track coordinates at SV
    MyKalRefit.push_back(TrTack.impactPointState().globalPosition().y());
    MyKalRefit.push_back(TrTack.impactPointState().globalPosition().z());
    MyKalRefit.push_back(TrTack.impactPointState().globalMomentum().x());//track momentum at SV
    MyKalRefit.push_back(TrTack.impactPointState().globalMomentum().y());//track momentum at SV
    MyKalRefit.push_back(TrTack.impactPointState().globalMomentum().z());//track momentum at SV
    double E=TMath::Sqrt(MyKalRefit.at(3)*MyKalRefit.at(3) + MyKalRefit.at(4)*MyKalRefit.at(4) + MyKalRefit.at(5)*MyKalRefit.at(5)+muon_mass*muon_mass);
    MyKalRefit.push_back(E);
      
  return MyKalRefit;
}

bool tracksMatchByDeltaR(const reco::Track *trk1, const reco::Track *trk2)
{
  //std::cout<<" trk 1 eta:"<<trk1->eta()<<" trk2 eta: "<<trk2->eta()<<std::endl;

  if ( reco::deltaR(*trk1, *trk2) < 3.e-1 && trk1->charge() == trk2->charge() ) {
    return true;
  }
  else return false;
}

// auxiliary function to exclude tracks associated to tau lepton decay "leg"                                                             
// from primary event vertex refit                                                                                               
 
typedef std::map<const reco::Track*, reco::TransientTrack> TransientTrackMap;
void removeTracks(TransientTrackMap& pvTracks_toRefit, const std::vector<reco::Track*> svTracks)
{
  for ( std::vector<reco::Track*>::const_iterator svTrack = svTracks.begin(); svTrack != svTracks.end(); ++svTrack ) {

    for ( TransientTrackMap::iterator pvTrack = pvTracks_toRefit.begin(); pvTrack != pvTracks_toRefit.end(); ++pvTrack ) {
      
      if ( tracksMatchByDeltaR(pvTrack->first, *svTrack) ) {
	//std::cout<<"pvTrack eta"<<pvTrack->first->eta()<<" sv track "<<svTrack->eta()<<std::endl;
	pvTracks_toRefit.erase(pvTrack);
	//break;
      }
    }
  }
}

// ------------ method called for each event  ------------
void
DsTau3MuAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

   //   edm::Handle<SimVertexContainer> simverts;
   //   iEvent.getByToken(simvertsToken_, simverts);

   edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);
   theTransientTrackBuilder_ = theTransientTrackBuilder.product();

   hEvt->Fill(1);
   hSize3MuCand->Fill(Cand3Mu->size());

   std::vector<reco::TransientTrack> pvTracks_original;
   TransientTrackMap pvTrackMap_refit;
   const reco::Vertex* eventVertex;

   for ( reco::Vertex::trackRef_iterator pvTrack = (*vertices)[0].tracks_begin(); pvTrack != (*vertices)[0].tracks_end(); ++pvTrack ) {
     reco::TransientTrack pvTrack_transient =theTransientTrackBuilder_->build(pvTrack->get());
     pvTracks_original.push_back(pvTrack_transient);
     pvTrackMap_refit.insert(std::make_pair(pvTrack->get(), pvTrack_transient));
  }

   cout<<" Number of tracks associated to the PV="<<pvTracks_original.size()<<endl;

   uint j=0;
   uint ngenP=genParticles->size();
   std::vector<int> genPidx;
   vector<reco::GenParticle> GenMu;
   //   GlobalPoint DsVtx (0.,0.,0.);
   float gvx, gvy, gvz;
   cout<<"****************GenLevel Info Begin********************"<<endl;
   for(edm::View<reco::GenParticle>::const_iterator gp=genParticles->begin(); gp!=genParticles->end(), j<ngenP; ++gp , ++j){
     if(fabs(gp->pdgId())==13) {
       for (uint i=0; i<gp->numberOfMothers();i++){
	 if(fabs(gp->mother(i)->pdgId())==15) {
	   std::cout<<j<<"--genMu pt="<<gp->pt()<<" eta="<<gp->eta()<<" phi="<<gp->phi()<<" pdgID="<<gp->pdgId()<<" tau pt="<<gp->mother(i)->pt()<<" mu vtx_x="<<gp->vx()<<" mu vtx_y="<<gp->vy()<<" mu vtx_z="<<gp->vz()<<endl;
	   genPidx.push_back(j);
	   GenMu.push_back(*gp);
	   if(fabs(gp->mother(i)->mother(0)->pdgId())==431) {
	     //GlobalPoint DsVtx (gp->mother(i)->mother(0)->vx(), gp->mother(i)->mother(0)->vy() , gp->mother(i)->mother(0)->vz());
	     gvx=gp->mother(i)->mother(0)->vx();
	     gvy=gp->mother(i)->mother(0)->vy();
	     gvz=gp->mother(i)->mother(0)->vz();
	   }
	 }
       }
     }
   }
   GlobalPoint DsVtx (gvx, gvy, gvz);
   cout<<"# Gen Muons: "<<genPidx.size()<<endl;
   cout<<"****************GenLevel Info End********************"<<endl;

   cout<<"Number of triplets="<<Cand3Mu->size()<<" number of muons="<<muons->size()<<endl;
   vector<int> ValidVtx;
   vector<RefCountedKinematicVertex>  VertexContainer_Valid, VertexContainer_ValidChi2;
   vector<TransientVertex> KVContainer, KV_ValidContainer;
   vector<CompositeCandidate> MatchedTripl, MatchedTripl2;

   for(edm::View<reco::CompositeCandidate>::const_iterator TauIt=Cand3Mu->begin(); TauIt!=Cand3Mu->end(); ++TauIt){

   vector<reco::Track*>KinTrackRef, KalTrackRef;     
   /*
     float mu1_pt=TauIt->daughter(0)->pt();
     float mu2_pt=TauIt->daughter(1)->pt();
     float mu3_pt=TauIt->daughter(2)->pt();
   */

   const Candidate * c1 = TauIt->daughter(0)->masterClone().get();
   const pat::Muon *mu1 = dynamic_cast<const pat::Muon *>(c1);

   const Candidate * c2 = TauIt->daughter(1)->masterClone().get();
   const pat::Muon *mu2 = dynamic_cast<const pat::Muon *>(c2);

   const Candidate * c3 = TauIt->daughter(2)->masterClone().get();
   const pat::Muon *mu3 = dynamic_cast<const pat::Muon *>(c3);


   bool isMatch1=false; bool isMatch2=false; bool isMatch3=false;
   if( (mu1->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu1->simMotherPdgId()) == 15) )     isMatch1=true;
   if( (mu2->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu2->simMotherPdgId()) == 15) )     isMatch2=true;
   if( (mu3->simType() == reco::MatchedMuonFromHeavyFlavour) && (fabs(mu3->simMotherPdgId()) == 15) )     isMatch3=true;


     bool isTauMatch = (isMatch1 && isMatch2 && isMatch3)  ;
     if (!isTauMatch) continue;
     /*
     cout<<" bestdr1="<<bestdr1<<" k="<<bestidx1<<" mu1match "<<mu1match<<"  GenMu_x="<<GenMu.at(bestidx1).vx()<<endl;
     cout<<" bestdr2="<<bestdr2<<" k="<<bestidx2<<" mu2match "<<mu2match<<"  GenMu_x="<<GenMu.at(bestidx2).vx()<<endl; 
     cout<<" bestdr3="<<bestdr3<<" k="<<bestidx3<<" mu3match "<<mu3match<<"  GenMu_x="<<GenMu.at(bestidx3).vx()<<endl;
     */
     TLorentzVector GenMu1, GenMu2, GenMu3, GenTau;
     GenMu1.SetPtEtaPhiM(mu1->simPt(), mu1->simEta(), mu1->simPhi(), 0.10565837);
     GenMu2.SetPtEtaPhiM(mu2->simPt(), mu2->simEta(), mu2->simPhi(), 0.10565837);
     GenMu3.SetPtEtaPhiM(mu3->simPt(), mu3->simEta(), mu3->simPhi(), 0.10565837);
     GenTau= GenMu1+GenMu2+GenMu3;
     //cout<<" gentau mass="<<GenTau.M()<<endl;

     //cout<<"Sim Ass Mu Vtx.x="<<mu1->vertex.x()<<endl;
     int bestidx1=0; int bestidx2=1; int  bestidx3=2;
     GlobalPoint GenMuVtx  (GenMu.at(bestidx1).vx(),  GenMu.at(bestidx1).vy(),  GenMu.at(bestidx1).vz());
     GlobalPoint GenMuVtx2  (GenMu.at(bestidx2).vx(),  GenMu.at(bestidx2).vy(),  GenMu.at(bestidx2).vz());
     GlobalPoint GenMuVtx3  (GenMu.at(bestidx3).vx(),  GenMu.at(bestidx3).vy(),  GenMu.at(bestidx3).vz());
     cout<<"GenMu1 vtx "<<GenMuVtx.x()<<" y="<<GenMuVtx.y()<<" z="<<GenMuVtx.z()<<endl;
     //cout<<" GenMu2 vtx "<<GenMuVtx2.x()<<" y="<<GenMuVtx2.y()<<" z="<<GenMuVtx2.z()<<endl;
     //cout<<" GenMu3 vtx "<<GenMuVtx3.x()<<" y="<<GenMuVtx3.y()<<" z="<<GenMuVtx3.z()<<endl;
     MatchedTripl.push_back(*TauIt);
     const Candidate * dau1 = TauIt->daughter( 0 );
     const Candidate * dau2 = TauIt->daughter( 1 );
     const Candidate * dau3 = TauIt->daughter( 2 );

     TrackRef trk1, trk2, trk3;
     if (dau1->isGlobalMuon()) { trk1 = dau1->get<TrackRef,reco::CombinedMuonTag>();} 
     else { trk1 = dau1->get<TrackRef>();}
     if (dau2->isGlobalMuon()) { trk2 = dau2->get<TrackRef,reco::CombinedMuonTag>();}
     else{ trk2 = dau2->get<TrackRef>();}
     if (dau3->isGlobalMuon()) { trk3 = dau3->get<TrackRef,reco::CombinedMuonTag>();}
     else{  trk3 = dau3->get<TrackRef>();}
     //cout<<" trk1 id="<<trk1.id()<<" tr2:"<<trk2.id()<<" trk3="<<trk3.id()<<endl;
     const reco::TransientTrack transientTrack1=theTransientTrackBuilder_->build( trk1 ); 
     const reco::TransientTrack transientTrack2=theTransientTrackBuilder_->build( trk2 ); 
     const reco::TransientTrack transientTrack3=theTransientTrackBuilder_->build( trk3 );
     reco::Track KinTrack1 =transientTrack1.track();
     reco::Track KinTrack2 =transientTrack2.track();
     reco::Track KinTrack3 =transientTrack3.track();
     reco::Track* KinTrackRef1=&KinTrack1;
     reco::Track* KinTrackRef2=&KinTrack2;
     reco::Track* KinTrackRef3=&KinTrack3;
     KinTrackRef.push_back(KinTrackRef1);     
     KinTrackRef.push_back(KinTrackRef2);     
     KinTrackRef.push_back(KinTrackRef3);     
     //cout<<" track eta="<<KinTrTracks.at(0).track().eta()<<endl;
     KinematicParticleFactoryFromTransientTrack pFactory;
     vector<RefCountedKinematicParticle> muonParticles;
     //initial chi2 and ndf before kinematic fits.                                                                                                                        
     ParticleMass muon_mass = 0.10565837; //pdg mass                                                                                                                                
     float muon_sigma = muon_mass*1.e-6;                                                                                                                                   

     
     float chi = 0.;
     float ndf = 0.;
     muonParticles.push_back( pFactory.particle(transientTrack1,muon_mass,chi,ndf,muon_sigma));
     muonParticles.push_back( pFactory.particle(transientTrack2,muon_mass,chi,ndf,muon_sigma));
     muonParticles.push_back( pFactory.particle(transientTrack3,muon_mass,chi,ndf,muon_sigma));          

     KinematicParticleVertexFitter fitter;
     RefCountedKinematicTree TauVertexFitTree;
     TauVertexFitTree = fitter.fit(muonParticles);


     if (TauVertexFitTree->isValid()) {

       if(bestidx1!=bestidx2 && bestidx2!=bestidx3 && bestidx1!=bestidx3){

	MatchedTripl2.push_back(*TauIt); 
       // VertexContainer_Valid
       ValidVtx.push_back(1);
       TauVertexFitTree->movePointerToTheTop();
 
       RefCountedKinematicParticle Tau_vFit_noMC = TauVertexFitTree->currentParticle();
       RefCountedKinematicVertex Tau_vFit_vertex_noMC = TauVertexFitTree->currentDecayVertex();
       
       //hValidKinFitVtx->Fill(Tau_vFit_vertex_noMC->isValid());

       float Chi2_KinFit = Tau_vFit_vertex_noMC->chiSquared();
       float NDof_KinFit = Tau_vFit_vertex_noMC->degreesOfFreedom();
       float svprob_KinFit = TMath::Prob(Chi2_KinFit, int(NDof_KinFit));

       //cout<<"3muon vtx: Kin Fit mass="<<Tau_vFit_noMC->currentState().mass()<<" Cand mass="<<TauIt->mass()<<" vtx_chi2="<<Tau_vFit_vertex_noMC->chiSquared()<<" Ndof="<<NDof_KinFit<<endl;
       


       hVtxChi2_Kin->Fill(Chi2_KinFit);
       hVtxNDof_Kin->Fill(NDof_KinFit );
       hProb_Kin->Fill(svprob_KinFit);

       TauVertexFitTree->movePointerToTheFirstChild();
       auto kmu1 = TauVertexFitTree->currentParticle();
       vector<float> RefMu1=RefitMu(kmu1);

       TauVertexFitTree->movePointerToTheNextChild();
       auto kmu2 = TauVertexFitTree->currentParticle();
       vector<float> RefMu2=RefitMu(kmu2);

       TauVertexFitTree->movePointerToTheNextChild();
       auto kmu3 = TauVertexFitTree->currentParticle();
       vector<float> RefMu3=RefitMu(kmu3);
       
       TLorentzVector RefKinMu1, RefKinMu2,RefKinMu3, RefTau;
       RefKinMu1.SetPxPyPzE(RefMu1.at(3), RefMu1.at(4), RefMu1.at(5), TMath::Sqrt(RefMu1.at(3)*RefMu1.at(3)+ RefMu1.at(4)*RefMu1.at(4) + RefMu1.at(5)*RefMu1.at(5) + RefMu1.at(6)*RefMu1.at(6) ));
       RefKinMu2.SetPxPyPzE(RefMu2.at(3), RefMu2.at(4), RefMu2.at(5), TMath::Sqrt(RefMu2.at(3)*RefMu2.at(3)+ RefMu2.at(4)*RefMu2.at(4) + RefMu2.at(5)*RefMu2.at(5) + RefMu2.at(6)*RefMu2.at(6) ));
       RefKinMu3.SetPxPyPzE(RefMu3.at(3), RefMu3.at(4), RefMu3.at(5), TMath::Sqrt(RefMu3.at(3)*RefMu3.at(3)+ RefMu3.at(4)*RefMu3.at(4) + RefMu3.at(5)*RefMu3.at(5) + RefMu3.at(6)*RefMu3.at(6) ));
       RefTau = RefKinMu1+RefKinMu2+RefKinMu3;



        
       cout<<"mu1: px="<<TauIt->daughter(0)->px()<<" ref px="<<RefMu1.at(3)<<" GenMuPx="<<GenMu1.Px()<<endl;
       cout<<"mu2: px="<<TauIt->daughter(1)->px()<<" ref px="<<RefMu2.at(3)<<" GenMuPx="<<GenMu2.Px()<<endl;
       cout<<"mu3: px="<<TauIt->daughter(2)->px()<<" ref px="<<RefMu3.at(3)<<" GenMuPx="<<GenMu3.Px()<<endl;
     


	 removeTracks(pvTrackMap_refit,  KinTrackRef);
	 std::vector<reco::TransientTrack> pvTracks_refit;
	 for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refit.begin();  pvTrack != pvTrackMap_refit.end(); ++pvTrack ) {
	   pvTracks_refit.push_back(pvTrack->second);}
	 cout<<" PV Tracks after refit="<<pvTracks_refit.size()<<endl;
	 /*for(uint i=0; i<pvTracks_refit.size(); i++){
	   TrackRef tr = TrackRef(pvTracks_refit, i);
	   //reco::Track pvTr=pvTracks_refit.at(i).track();
	   //TrackRef pvTrRef = pvTr.get<TrackRef>();
	   cout<<i<<"PV track ID="<<tr.id()<<endl;
	   }*/

	 KalmanVertexFitter PV_fitter (true);
	 TransientVertex PVertex = PV_fitter.vertex(pvTracks_refit);
	 //CachingVertex<5> fittedVertex = vertexFitter.vertex(tracksToVertex);                                                 
	 GlobalPoint PVertexPos  (PVertex.position());
	 GlobalPoint SVertexPos  (Tau_vFit_vertex_noMC->position());
	 double FlightDist = TMath::Sqrt( pow(( PVertexPos.x() -SVertexPos.x()),2)+ pow(( PVertexPos.y() -SVertexPos.y()),2) + pow(( PVertexPos.z() -SVertexPos.z()),2));
	 
	 VertexDistance3D vertTool;
	 VertexState PVstate(PVertex.position(),PVertex.positionError());
	 VertexState SVstate(Tau_vFit_vertex_noMC->position(),Tau_vFit_vertex_noMC->error());
	 double distance = vertTool.distance(PVstate, SVstate).value();
	 double dist_err = vertTool.distance(PVstate, SVstate).error();
	 double dist_sign =vertTool.distance(PVstate, SVstate).significance();
	 double chi2 = vertTool.compatibility(PVstate, SVstate);

	 hDist_PVSV_Kin->Fill(distance);
	 hDistErr_PVSV_Kin->Fill(dist_err);
	 hDistSign->Fill(dist_sign);
	 hDistComp_PVSV_Kin->Fill(chi2);

	 double GenFlightDist = TMath::Sqrt( pow((DsVtx.x()-GenMu.at(bestidx1).vx()),2)+pow((DsVtx.y()-GenMu.at(bestidx1).vy()),2) +pow((DsVtx.z()-GenMu.at(bestidx1).vz()),2));
	 double  GenMuSV =  TMath::Sqrt(GenMu.at(bestidx1).vx()*GenMu.at(bestidx1).vx()+GenMu.at(bestidx1).vy()*GenMu.at(bestidx1).vy()+GenMu.at(bestidx1).vz()*GenMu.at(bestidx1).vz());
	 double Ref_DistPVSV = TMath::Sqrt( pow(RefMu1.at(0)-PVertexPos.x(),2)+ pow(RefMu1.at(1)-PVertexPos.y(),2)+ pow(RefMu1.at(2)-PVertexPos.z(),2)  );
	 //const edm::SimVertexContainer simVC = *(simverts.product());
	 cout<<"GenPV : x="<<DsVtx.x()<<" y="<<DsVtx.y()<<" z="<<DsVtx.z()<<endl;
	 //cout<<"SimPV : x="<<simVC.begin()->position().x()<<endl;
	 cout<<"RecoPV: x="<<PVertexPos.x()<<" y="<<PVertexPos.y()<<" z="<<PVertexPos.z()<<endl;

	 cout<<"GenSV : x="<<GenMu.at(bestidx1).vx()<<" y="<<GenMu.at(bestidx1).vy()<<" z="<<GenMu.at(bestidx1).vz()<<endl;
	 cout<<"RecosV: x="<<SVertexPos.x()<<" y="<<SVertexPos.y()<<" z="<<SVertexPos.z()<<endl;

	 cout<<" Dist="<<distance<<" distErr="<<dist_err<<" GenDist="<<GenMuSV<<" MyFlightDist="<<FlightDist<<" GenFlightDist"<<GenFlightDist<<" refitted="<<Ref_DistPVSV<<endl;

	 hGenSV_KinSV->Fill(GenMuSV);
	 hDeltaGenSVKinSV->Fill(GenMuSV-distance);
	 hVtxMass_KinVsGen->Fill(RefTau.M()-GenTau.M());
	 hVtxMass_KinVsCand->Fill(Tau_vFit_noMC->currentState().mass()-TauIt->mass());
	 hVtxMass_RecoVsGen->Fill(TauIt->mass()-GenTau.M());

	 hDxGenReco->Fill(TauIt->daughter(0)->vx()-GenMu.at(bestidx1).vx());
	 hDxGenReco->Fill(TauIt->daughter(1)->vx()-GenMu.at(bestidx2).vx());
	 hDxGenReco->Fill(TauIt->daughter(2)->vx()-GenMu.at(bestidx3).vx());

	 hDyGenReco->Fill(TauIt->daughter(0)->vy()-GenMu.at(bestidx1).vy());
	 hDyGenReco->Fill(TauIt->daughter(1)->vy()-GenMu.at(bestidx2).vy());
	 hDyGenReco->Fill(TauIt->daughter(2)->vy()-GenMu.at(bestidx3).vy());

       hDzGenReco->Fill(TauIt->daughter(0)->vz()-GenMu.at(bestidx1).vz());
       hDzGenReco->Fill(TauIt->daughter(1)->vz()-GenMu.at(bestidx2).vz());
       hDzGenReco->Fill(TauIt->daughter(2)->vz()-GenMu.at(bestidx3).vz());
	 

       hDxGenRefit->Fill(RefMu1.at(0)-GenMu.at(bestidx1).vx());
       hDxGenRefit->Fill(RefMu2.at(0)-GenMu.at(bestidx2).vx());
       hDxGenRefit->Fill(RefMu3.at(0)-GenMu.at(bestidx3).vx());

       hDyGenRefit->Fill(RefMu1.at(1)-GenMu.at(bestidx1).vy());
       hDyGenRefit->Fill(RefMu2.at(1)-GenMu.at(bestidx2).vy());
       hDyGenRefit->Fill(RefMu3.at(1)-GenMu.at(bestidx3).vy());

       hDzGenRefit->Fill(RefMu1.at(2)-GenMu.at(bestidx1).vz());
       hDzGenRefit->Fill(RefMu2.at(2)-GenMu.at(bestidx2).vz());
       hDzGenRefit->Fill(RefMu3.at(2)-GenMu.at(bestidx3).vz());
	 
       
     hDx->Fill(TauIt->daughter(0)->vx()-RefMu1.at(0));
     hDx->Fill(TauIt->daughter(1)->vx()-RefMu2.at(0));
     hDx->Fill(TauIt->daughter(2)->vx()-RefMu3.at(0));

     hDy->Fill(TauIt->daughter(0)->vy()-RefMu1.at(1));
     hDy->Fill(TauIt->daughter(1)->vy()-RefMu2.at(1));
     hDy->Fill(TauIt->daughter(2)->vy()-RefMu3.at(1));

     hDz->Fill(TauIt->daughter(0)->vz()-RefMu1.at(2));
     hDz->Fill(TauIt->daughter(1)->vz()-RefMu2.at(2));
     hDz->Fill(TauIt->daughter(2)->vz()-RefMu3.at(2));

     hDPx_GenKin->Fill(GenMu1.Px() - RefMu1.at(3));
     hDPx_GenKin->Fill(GenMu2.Px() - RefMu2.at(3));
     hDPx_GenKin->Fill(GenMu3.Px() - RefMu3.at(3));

     hDPy_GenKin->Fill(GenMu1.Py() - RefMu1.at(4));
     hDPy_GenKin->Fill(GenMu2.Py() - RefMu2.at(4));
     hDPy_GenKin->Fill(GenMu3.Py() - RefMu3.at(4));

     hDPz_GenKin->Fill(GenMu1.Pz() - RefMu1.at(5));
     hDPz_GenKin->Fill(GenMu2.Pz() - RefMu2.at(5));
     hDPz_GenKin->Fill(GenMu3.Pz() - RefMu3.at(5));

       }
     }//Valid Kin Fit

     //Kalman VertexFit
     std::vector<reco::TransientTrack> Ttracks;
     Ttracks.push_back(transientTrack1);
     Ttracks.push_back(transientTrack2);
     Ttracks.push_back(transientTrack3);
     
     KalmanVertexFitter KVfitter (true);
     TransientVertex KVertex = KVfitter.vertex(Ttracks);
     //CachingVertex<5> fittedVertex = vertexFitter.vertex(tracksToVertex);
     KVContainer.push_back(KVertex);

     if(KVertex.hasRefittedTracks() && isTauMatch ){
     if (KVertex.isValid()) KV_ValidContainer.push_back(KVertex);
     float Chi2_KalFit = KVertex.totalChiSquared();
     float NDof_KalFit = KVertex.degreesOfFreedom();
     float svprob_KalFit = TMath::Prob(Chi2_KalFit, int(NDof_KalFit));

     //hVtxMass_KalVsCand->Fill(KVertex.mass() - TauIt->mass());
     hVtxChi2_Kal->Fill(Chi2_KalFit);
     hVtxNDof_Kal->Fill(NDof_KalFit );
     hProb_Kal->Fill(svprob_KalFit);

     vector < TransientTrack > ttrks = KVertex.refittedTracks();      
     TLorentzVector KMu1, KMu2, KMu3, KTau;
     KMu1.SetPxPyPzE( KalmanRefMu(ttrks.at(0)).at(3), KalmanRefMu(ttrks.at(0)).at(4), KalmanRefMu(ttrks.at(0)).at(5), KalmanRefMu(ttrks.at(0)).at(6));
     KMu2.SetPxPyPzE( KalmanRefMu(ttrks.at(1)).at(3), KalmanRefMu(ttrks.at(1)).at(4), KalmanRefMu(ttrks.at(1)).at(5), KalmanRefMu(ttrks.at(1)).at(6)); 
     KMu3.SetPxPyPzE( KalmanRefMu(ttrks.at(2)).at(3), KalmanRefMu(ttrks.at(2)).at(4), KalmanRefMu(ttrks.at(2)).at(5), KalmanRefMu(ttrks.at(2)).at(6));
     KTau = KMu1+KMu2+KMu3;
     hVtxMass_KalVsCand->Fill(KTau.M() - TauIt->mass()); 
     hVtxMass_KalVsGen->Fill(KTau.M() - GenTau.M());
     cout<<"Mu1: GenPx="<<GenMu1.Px()<<" Kalman px="<<KalmanRefMu(ttrks.at(0)).at(3)<<" reco px="<<dau1->px()<<endl;
     cout<<"Mu2: GenPx="<<GenMu2.Px()<<" Kalman px="<<KalmanRefMu(ttrks.at(1)).at(3)<<" reco px="<<dau2->px()<<endl;
     cout<<"Mu3: GenPx="<<GenMu3.Px()<<" Kalman px="<<KalmanRefMu(ttrks.at(2)).at(3)<<" reco px="<<dau3->px()<<endl;

       reco::Track KalTrack1 =ttrks.at(0).track();
       reco::Track KalTrack2 =ttrks.at(1).track();
       reco::Track KalTrack3 =ttrks.at(2).track();
       reco::Track* KalTrackRef1=&KalTrack1;
       reco::Track* KalTrackRef2=&KalTrack2;
       reco::Track* KalTrackRef3=&KalTrack3;
       KalTrackRef.push_back(KalTrackRef1);
       KalTrackRef.push_back(KalTrackRef2);
       KalTrackRef.push_back(KalTrackRef3);
       removeTracks(pvTrackMap_refit,  KalTrackRef);
       std::vector<reco::TransientTrack> pvTracks_refit2;
       for ( TransientTrackMap::iterator pvTrack = pvTrackMap_refit.begin();  pvTrack != pvTrackMap_refit.end(); ++pvTrack ) {
	 pvTracks_refit2.push_back(pvTrack->second);}
       cout<<" Kalman PV Tracks after refit="<<pvTracks_refit2.size()<<endl;

       KalmanVertexFitter PV_fitter2 (true);
       TransientVertex PVertex2 = PV_fitter2.vertex(pvTracks_refit2);

       VertexDistance3D vertTool2;
       VertexState PVstate2(PVertex2.position(),PVertex2.positionError());
       VertexState SVstate2(KVertex.position(),KVertex.positionError());
       double distanceK = vertTool2.distance(PVstate2, SVstate2).value();
       double dist_errK = vertTool2.distance(PVstate2, SVstate2).error();
       double chi2K = vertTool2.compatibility(PVstate2, SVstate2);

       hDist_PVSV_Kal->Fill(distanceK);
       hDistErr_PVSV_Kal->Fill(dist_errK);
       hDistComp_PVSV_Kal->Fill(chi2K);

       double GenMuSV2 =  TMath::Sqrt(GenMu.at(bestidx1).vx()*GenMu.at(bestidx1).vx()+GenMu.at(bestidx1).vy()*GenMu.at(bestidx1).vy()+GenMu.at(bestidx1).vz()*GenMu.at(bestidx1).vz());
       hGenSV_KalSV->Fill(GenMuSV2);
       hDeltaGenSVKalSV->Fill(GenMuSV2-distanceK);

       hDxGenKal->Fill(KVertex.position().x()-GenMu.at(bestidx1).vx());
       hDyGenKal->Fill(KVertex.position().y()-GenMu.at(bestidx1).vy());
       hDzGenKal->Fill(KVertex.position().z()-GenMu.at(bestidx1).vz());

       //cout<<"Kalman Dist="<<distanceK<<" distErr="<<dist_errK<<" GenDist="<<GenMuSV2<<endl;
       hDPx_GenKal->Fill(GenMu1.Px()-KMu1.Px());
       hDPx_GenKal->Fill(GenMu2.Px()-KMu2.Px());
       hDPx_GenKal->Fill(GenMu3.Px()-KMu3.Px());

       hDPy_GenKal->Fill(GenMu1.Py()-KMu1.Py());
       hDPy_GenKal->Fill(GenMu2.Py()-KMu2.Py());
       hDPy_GenKal->Fill(GenMu3.Py()-KMu3.Py());

       hDPz_GenKal->Fill(GenMu1.Pz()-KMu1.Pz());
       hDPz_GenKal->Fill(GenMu2.Pz()-KMu2.Pz());
       hDPz_GenKal->Fill(GenMu3.Pz()-KMu3.Pz());

       hDPx_GenReco->Fill(GenMu1.Px()-dau1->px());
       hDPx_GenReco->Fill(GenMu2.Px()-dau2->px());
       hDPx_GenReco->Fill(GenMu3.Px()-dau3->px());

       hDPy_GenReco->Fill(GenMu1.Py()-dau1->py());
       hDPy_GenReco->Fill(GenMu2.Py()-dau2->py());
       hDPy_GenReco->Fill(GenMu3.Py()-dau3->py());

       hDPz_GenReco->Fill(GenMu1.Pz()-dau1->pz());
       hDPz_GenReco->Fill(GenMu2.Pz()-dau2->pz());
       hDPz_GenReco->Fill(GenMu3.Pz()-dau3->pz());

       //TrajectoryStateClosestToPoint TSCP1 = ttrks[0].trajectoryStateClosestToPoint(KVertex.position());
       //TrajectoryStateClosestToPoint two_TSCP = tracksToVertex[1].trajectoryStateClosestToPoint(fittedVertex.position());
       //GlobalVector one_momentum = TSCP1.momentum();
       //GlobalVector two_momentum = two_TSCP.momentum();

       // calculate mass
       //double total_energy = sqrt(one_momentum.mag2() + 0.14*0.14) + sqrt(two_momentum.mag2() + 0.14*0.14);
       //double total_px = one_momentum.x();
       //double total_py = one_momentum.y() + two_momentum.y();
       //double total_pz = one_momentum.z() + two_momentum.z();
       //double mass = sqrt(pow(total_energy, 2) - pow(total_px, 2) - pow(total_py, 2) - pow(total_pz, 2));
       
     }
   }//loop over comp cand
 
   cout<<" Matched triplets: "<<MatchedTripl.size()<<"  Arbitration="<<MatchedTripl2.size()<<endl;
   hValidVtx_size_Kalman->Fill(KV_ValidContainer.size());
}


// ------------ method called once each job just before starting event loop  ------------
void 
DsTau3MuAna::beginJob()
{
  hEvt= fs->make<TH1F>("hEvt","# Events",3,0,3);
  hSize3MuCand = fs->make<TH1F>("hSize3MuCand","Number of Triplets",200, 0, 200);
  hVtxMass_KinVsCand =  fs->make<TH1F>("hVtxMass_KinVsCand","Three muons Mass_{KinFit}-Mass_{inv}",200,-0.1,0.1);
  hVtxMass_KinVsGen=  fs->make<TH1F>("hVtxMass_KinVsGen","Three muons Mass_{KinFit}-Mass_{gen}",200,-0.1,0.1);
  hVtxMass_RecoVsGen=  fs->make<TH1F>("hVtxMass_RecoVsGen","Three muons Mass_{reco}-Mass_{gen}",200,-0.1,0.1);

  hValidKinFitVtx = fs->make<TH1F>("hValidKinFitVtx","Valid KinFit Vtx",2,0,1);
  hVtxChi2_Kin = fs->make<TH1F>("hVtxChi2_Kin","Vtx Chi2 Kin Fit",200,0,200);
  hVtxNDof_Kin = fs->make<TH1F>("hVtxNDof_Kin","Vtx # DOF Kin Fit",5,0,5);
  hProb_Kin = fs->make<TH1F>("hProb_Kin","Vtx Prob Kin Fit",1000,0,1);
  hDx = fs->make<TH1F>("hDx","x_{mu}-x_{muRef}",100,-5,5);
  hDy = fs->make<TH1F>("hDy","y_{mu}-y_{muRef}",100,-5,5);
  hDz = fs->make<TH1F>("hDz","z_{mu}-z_{muRef}",100,-5,5);

  hDxGenReco = fs->make<TH1F>("hDxGenReco","x_{gen}-x_{reco}",100,-5,5);
  hDyGenReco = fs->make<TH1F>("hDyGenReco","y_{gen}-y_{reco}",100,-5,5);
  hDzGenReco = fs->make<TH1F>("hDzGenReco","z_{gen}-z_{reco}",100,-5,5);

  hDxGenRefit = fs->make<TH1F>("hDxGenRefit","x_{Gen}-x_{Kin}",100,-5,5);
  hDyGenRefit = fs->make<TH1F>("hDyGenRefit","y_{Gen}-y_{Kin}",100,-5,5);
  hDzGenRefit = fs->make<TH1F>("hDzGenRefit","z_{Gen}-z_{Kin}",100,-5,5);

  hDxGenKal = fs->make<TH1F>("hDxGenKal","x_{Gen}-x_{Kal}",100,-5,5);
  hDyGenKal = fs->make<TH1F>("hDyGenKal","y_{Gen}-y_{Kal}",100,-5,5);
  hDzGenKal = fs->make<TH1F>("hDzGenKal","z_{Gen}-z_{Kal}",100,-5,5);

  hDPx_GenKal = fs->make<TH1F>("hDPx_GenKal", "px residuals(gen, kalman refit)", 200,-0.5,0.5);
  hDPx_GenKin = fs->make<TH1F>("hDPx_GenKin",  "px residuals(gen, kin refit)", 200,-0.5,0.5);
  hDPx_GenReco=fs->make<TH1F>("hDPx_GenReco", "px residuals(gen, reco)",200,-0.5,0.5);

  hDPy_GenKal  = fs->make<TH1F>("hDPy_GenKal", "py residuals(gen, kalman refit)", 200,-0.5,0.5);
  hDPy_GenKin = fs->make<TH1F>("hDPy_GenKin",  "py residuals(gen, kin refit)", 200,-0.5,0.5);
  hDPy_GenReco=fs->make<TH1F>("hDPy_GenReco", "py residuals(gen, reco)",200,-0.5,0.5);

  hDPz_GenKal = fs->make<TH1F>("hDPz_GenKal", "pz residuals(gen, kalman refit)", 200,-0.5,0.5);
  hDPz_GenKin = fs->make<TH1F>("hDPz_GenKin",  "pz residuals(gen, kin refit)", 200,-0.5,0.5);
  hDPz_GenReco=fs->make<TH1F>("hDPz_GenReco", "pz residuals(gen, reco)",200,-0.5,0.5);

  hDist_PVSV_Kin= fs->make<TH1F>("hDist_PVSV_Kin", "hDist_PVSV_Kin", 250, 0, 5);
  hDistErr_PVSV_Kin= fs->make<TH1F>("hDistErr_PVSV_Kin", "hDistErr_PVSV_Kin", 100, 0, 0.1);
  hDistComp_PVSV_Kin=fs->make<TH1F>("hDistComp_PVSV_Kin", "hDistComp_PVSV_Kin", 100, 0, 100);
  hDistSign=fs->make<TH1F>("hDistSign", "hDistSign", 100, -1, 1);

  hDist_PVSV_Kal= fs->make<TH1F>("hDist_PVSV_Kal", "hDist_PVSV_Kal", 250, 0, 5);
  hDistErr_PVSV_Kal= fs->make<TH1F>("hDistErr_PVSV_Kal", "hDistErr_PVSV_Kal", 100, 0, 0.1);
  hDistComp_PVSV_Kal=fs->make<TH1F>("hDistComp_PVSV_Kal", "hDistComp_PVSV_Kal", 100, 0, 100);

  hVtxMass_KalVsCand = fs->make<TH1F>("hVtxMass_KalVsCand","Three muons Mass_{KinFit}-Mass_{inv}",200,-0.5,0.5);
  hVtxMass_KalVsGen= fs->make<TH1F>("hVtxMass_KalVsGen", "Three muons Mass_{KalFit}-Mass_{Gen}",200,-0.5,0.5);
  hVtxChi2_Kal = fs->make<TH1F>("hVtxChi2_Kal","Vtx Chi2 Kin Fit",200,0,200);
  hVtxNDof_Kal = fs->make<TH1F>("hVtxNDof_Kal","Vtx # DOF Kin Fit",5,0,5);
  hProb_Kal = fs->make<TH1F>("hProb_Kal","Vtx Prob Kin Fit",1000,0,1);
  hValidVtx_size_Kalman = fs->make<TH1F>("hValidVtx_size_Kalman","hValidVtx_size_Kalman",100,0,100);
  hGenSV_KinSV = fs->make<TH1F>("hGenSV_KinSV", "hGenSV_KinSV", 250, 0, 5);
  hDeltaGenSVKinSV = fs->make<TH1F>("hDeltaGenSVKinSV", "hDeltaGenSVKinSV", 200, -1, 1);

  hGenSV_KalSV = fs->make<TH1F>("hGenSV_KalSV","hGenSV_KalSV",250, 0, 5);
  hDeltaGenSVKalSV = fs->make<TH1F>("hDeltaGenSVKalSV","hDeltaGenSVKalSV", 200, -1, 1);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
DsTau3MuAna::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DsTau3MuAna::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DsTau3MuAna);
