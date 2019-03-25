import FWCore.ParameterSet.Config as cms

process = cms.Process('Tau3MuSkim')

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
process.load("MiniAna2017.MiniAna2017Tree.Tau3MuSkimAOD_cff")
#Tau3MuSkimAODForSync_cff.py
#process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6' #data2017
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6' #mc2016
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_RealisticBS_25ns_13TeV2016_v1_mc' #mc2016 gives proper mass distr
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_530.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_531.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_532.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_533.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_534.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_535.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_536.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_537.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_538.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_539.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_540.root',
         'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_541.root',

    )
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histoSkim.root"))


process.Tau3MuSkim = cms.Path(process.ThreeMuonSelSeq 
                     )


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
        'keep *_PileupInfo_*_*',
        'keep *_genParticles_*_*',
        'keep recoVertex_offlinePrimaryVertices_*_*',
        'keep *_generalTracks_*_RECO',
        #'keep *_g4SimHits_*_*',
        'keep PSimHits_g4SimHits_MuonCSCHits_SIM',
        'keep PSimHits_g4SimHits_MuonDTHits_SIM',
        'keep PSimHits_g4SimHits_MuonRPCHits_SIM',
        'keep SimVertexs_g4SimHits__SIM',
        'keep *_csc2DRecHits_*_*',
        'keep *_dt1DRecHits_*_*',
        'keep *_rpcRecHits_*_*',
        'keep SimVertexs_g4SimHits_*_*',
        'drop floats_generalTracks_MVAValues_RECO'        
        #'keep recoTracks_generalTracks_*_*',

        )
)

#type_label_instance_process
process.outpath = cms.EndPath(process.out) 



"""
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
from PhysicsTools.PatAlgos.tools.trigTools import *
patTriggerEvent.patTriggerMatches  = cms.VInputTag( "muonTriggerMatchHLTMuons" )
switchOnTrigger( process )
switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ] )
"""
