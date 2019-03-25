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
process.load("SkimTools.SkimTau3Mu.Tau3MuSkimAOD_cff")

#Tau3MuSkimAODForSync_cff.py
#process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6' #data2017
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6' #mc2016
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_RealisticBS_25ns_13TeV2016_v1_mc' #mc2016 gives proper mass distr
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1010.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1008.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1007.root',  
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1006.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1005.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1004.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1003.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1002.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1001.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1000.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_999.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_998.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_997.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_996.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_995.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_994.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_993.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_992.root',
         'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_991.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_1.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_10.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_100.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_101.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_102.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_103.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_104.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_105.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_106.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_107.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_108.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_109.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_11.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_110.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_111.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_112.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_113.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_114.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_115.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_116.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_117.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_118.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_119.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_12.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_120.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_121.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_122.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_123.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_124.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_125.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_126.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_127.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_128.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_129.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_13.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_130.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_131.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_132.root',
'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0000/custom_DsTau3Mu_13TeV_RECO_crab350_133.root',


        #'file:./B0A2E5AE-60AF-E811-8BF5-0CC47A7E6A5C.root' #Run2017F AOD
        #'/store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/ECAA9E63-C5FD-E711-9C1E-008CFAE451DC.root'
        #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/D699CEF5-E1FA-E711-A7EC-02163E013399.root'
        #'root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/MINIAOD/17Nov2017-v1/00000/129CB4BC-17FD-E711-8D49-FA163E838299.root'

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

#vector<PileupSummaryInfo>             "addPileupInfo"             ""                "HLT"
#vector<reco::TrackExtra>              "generalTracks"             ""                "RECO"
#vector<reco::TrackExtra>              "globalMuons"               ""                "RECO"
#vector<reco::TrackExtra>              "refittedStandAloneMuons"   ""                "RECO"
#vector<reco::TrackExtra>              "standAloneMuons"           ""                "RECO"
#        vector<reco::TrackExtra>              "generalTracks"             ""                "RECO"
        #'keep *_g4SimHits_*_*',


#type_label_instance_process
process.outpath = cms.EndPath(process.out) 



"""
from PhysicsTools.PatAlgos.triggerLayer1.triggerEventProducer_cfi import *
from PhysicsTools.PatAlgos.tools.trigTools import *
patTriggerEvent.patTriggerMatches  = cms.VInputTag( "muonTriggerMatchHLTMuons" )
switchOnTrigger( process )
switchOnTriggerMatching( process, triggerMatchers = [ 'muonTriggerMatchHLTMuons' ] )
"""
