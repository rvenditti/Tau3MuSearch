import FWCore.ParameterSet.Config as cms

process = cms.Process('DsPhiPiSkim')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('DsPhiPiTreeMaker.DsPhiPiTreeMaker.DsPhiPiSkimAOD_cff')
#from DsPhiPiTreeMaker.DsPhiPiTreeMaker.DsPhiPiSkimAOD_cff import *

#Tau3MuSkimAODForSync_cff.py
#process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6' #data2017
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6' #mc2016
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_RealisticBS_25ns_13TeV2016_v1_mc' #mc2016 gives proper mass distr
process.GlobalTag.globaltag = '94X_dataRun2_v11'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/data/Run2017B/DoubleMuonLowMass/AOD/23Jun2017-v1/90000/FC4BB0C3-E358-E711-9EFE-0025904B739A.root',
        'root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/AOD/09May2018-v1/80000/FEDF5D97-BEB0-E811-95BF-0CC47AD98B94.root',
        'root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/AOD/09May2018-v1/80000/AECC4C56-BAB0-E811-B92A-008CFA1979AC.root'
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_10.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_11.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_12.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_13.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_14.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_15.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_16.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_17.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_18.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_19.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_20.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_21.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_22.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_23.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_24.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_25.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_26.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_27.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_28.root',
             #'file:/lustre/cms/store/user/fsimone/DsPhiPi/crab_crab_DsPhiPi_13TeV_RECO/190225_141742/0000/DsPhiPi_13TeV_RECO_29.root',
    )
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Tree_PhiPi.root"))


process.Tree3Mu = cms.EDAnalyzer("DsPhiPiTreeMaker",
                                 isMcLabel = cms.untracked.bool(False),
                              #is3MuLabel = cms.untracked.bool(False),
                              muonLabel=cms.InputTag("looseMuons"),
                              VertexLabel=cms.InputTag("offlinePrimaryVerticesWithBS"),
                              TracksLabel=cms.InputTag("LooseTrack"),
                              genParticleLabel=cms.InputTag("genParticles"),
                              #Cand3MuLabel=cms.InputTag("ThreeMuonsVtxKalmanFit"),
                              Cand2Mu1TrackLabel=cms.InputTag("TwoMuonsOneTrackKalmanVtxFit"),
                              DiMuonLabel=cms.InputTag("DiMuonsVtxFit"),
                              pileupSummary = cms.InputTag("addPileupInfo"),
                              triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                              triggerSummary = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                              AlgInputTag = cms.InputTag( "gtStage2Digis" )
)



process.DsPhiPiSkim = cms.Path( process.TwoMuOneTrackSelSeq
                              * process.Tree3Mu)

"""
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("file_AODSIM_test_ForTreeMaker.root"),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('Tau3MuSkim')),
                               outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_looseMuons_*_Tau3MuSkim',
        'keep recoVertexCompositeCandidates_*_*_Tau3MuSkim',
        'keep recoCompositeCandidates_*_*_Tau3MuSkim',
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_generator_*_*',
        'keep *_offlineBeamSpot_*_*',
        'keep recoMuon_muons_*_*',
        'keep *_TriggerResults_*_*',
        'keep *_gtStage2Digis_*_*',
        'keep *_gmtStage2Digis_Muon_*',
        'keep *_scalersRawToDigi_*_*',
        'keep *_Trigger_*_*',
        'keep *_addPileupInfo_*_*',
        'keep *_genParticles_*_*',
        'keep recoVertex_offlinePrimaryVertices_*_*',
        'keep *_generalTracks_*_RECO',
        #'keep TrackExtra_generalTracks_*_RECO',
        'keep *_globalMuons_*_RECO',
        'keep *_standAloneMuons_*_RECO',
        'keep PSimHits_g4SimHits_MuonCSCHits_SIM',
        'keep PSimHits_g4SimHits_MuonDTHits_SIM',
        'keep PSimHits_g4SimHits_MuonRPCHits_SIM',
        'keep SimVertexs_g4SimHits__SIM',
        'keep *_csc2DRecHits_*_*',
        'keep *_dt1DRecHits_*_*',
        'keep *_rpcRecHits_*_*',
        'keep SimVertexs_g4SimHits_*_*',
        'drop floats_generalTracks_MVAValues_RECO',        
        'drop recoTrack_globalMuons_*_RECO',
        'drop recoTrack_standAloneMuons_*_RECO',
        #'keep recoTracks_generalTracks_*_*',

        )
)
"""
#type_label_instance_process
#process.outpath = cms.EndPath(process.out) 



"""
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
from PhysicsTools.PatAlgos.tools.trigTools import *
patTriggerEvent.patTriggerMatches  = cms.VInputTag( "muonTriggerMatchHLTMuons" )
switchOnTrigger( process )
switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ] )
"""
