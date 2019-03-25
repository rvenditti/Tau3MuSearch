import FWCore.ParameterSet.Config as cms

import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *


PatMuons = patMuons.clone(
    src = cms.InputTag("muons"),
    useParticleFlow = cms.bool( False ),
    #embedHighLevelSelection = cms.bool(True),
    computeMiniIso = cms.bool(False), 
    computeMuonMVA= cms.bool(False),
    computeSoftMuonMVA = cms.bool(True),
    addTriggerMatching = cms.bool(False),
    addGenMatch   = cms.bool(False),
    embedGenMatch = cms.bool(True),
)



looseMuons = cms.EDFilter("PATMuonSelector",
                          src = cms.InputTag("PatMuons"),
                           cut = cms.string('pt>0.5 && abs(eta)<2.4 && (innerTrack().isNonnull) && (charge!=0) && (innerTrack().hitPattern().numberOfValidPixelHits()>0)'), 
                          filter = cms.bool(False)                                
)

ThreeMuonsCand = cms.EDProducer("CandViewShallowCloneCombiner",
                         checkCharge = cms.bool(False),
                         cut = cms.string(''),
                         decay = cms.string("looseMuons looseMuons looseMuons")
)

ThreeMuonsVtxKalmanFit = cms.EDProducer("KalmanVertexFitCompositeCandProducer",
                                        src = cms.InputTag("ThreeMuonsCand"),
                                        #cut = cms.string('mass <5'), 
                                        #cut = cms.string('(vertexChi2 < 40) && (vertexNdof == 3) && (mass <5)'),
                                        ##filter = cms.bool(True)
                                        )


########################Define Histograms########################
InitialPlots = cms.EDAnalyzer('SimpleEventCounter',
                                   muonsInputTag = cms.InputTag("muons"),
                                   )

PlotsAfterLooseMu = cms.EDAnalyzer('RecoMuonAnalyzer', muonsInputTag = cms.InputTag("looseMuons"),)
PlotsAfter3Mu     = cms.EDAnalyzer('RecoMuonAnalyzer', muonsInputTag = cms.InputTag("looseMuons"),)
PlotsAfter3MuVtx  = cms.EDAnalyzer('RecoMuonAnalyzer', muonsInputTag = cms.InputTag("looseMuons"),)



PatMuonSequence = cms.Sequence(InitialPlots *
                               PatMuons *
                               looseMuons *
                               PlotsAfterLooseMu *
                               ThreeMuonsCand *
                               PlotsAfter3Mu *
                               ThreeMuonsVtxKalmanFit *
                               PlotsAfter3MuVtx
                               )








