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
#process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.GlobalTag.globaltag = '94X_mc2017_realistic_v17'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'root://xrootd-cms.infn.it//store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/110000/42304418-1347-E911-828D-003048F5B2F1.root',
         #'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_1.root', 
         #'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_10.root', 
         #'file:/lustre/cms/store/user/fsimone/B0Tau3Mu/crab_crab_B0Tau3Mu_13TeV_RECO/190225_140205/0000/B0Tau3Mu_13TeV_RECO_100.root'

         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/BE4475A5-1848-E911-BFAA-801844E561B0.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/BC850D08-A347-E911-A461-0CC47ADAF60A.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/BC2C5249-1B48-E911-8813-008CFAC94184.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/B6A9486A-7C4E-E911-8773-001E67DDC0FB.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/B2469434-AF4E-E911-8D81-90B11CBCFF9C.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/B08EF594-934E-E911-9FBD-F8C288DA8479.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/A89C54F2-2847-E911-BA06-1866DAEA6D08.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/A4E4757A-E448-E911-838A-D067E5F91FF8.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/98D51CE8-E147-E911-B51A-002590490020.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/8E0C741C-3848-E911-8ABA-0242AC130002.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/8C7C3FD5-AF47-E911-8BA2-002590491B1E.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/86DF70E9-6148-E911-B58D-200009DEFE80.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/8638DD7A-964E-E911-8FDA-3417EBE47FCA.root',
         '/store/mc/RunIIFall17DRPremix/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/120000/82338802-5651-E911-8F2C-44A842CF058B.root',


        #'/store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/ECAA9E63-C5FD-E711-9C1E-008CFAE451DC.root'
        #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/D699CEF5-E1FA-E711-A7EC-02163E013399.root'
        #'root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/MINIAOD/17Nov2017-v1/00000/129CB4BC-17FD-E711-8D49-FA163E838299.root'

    )
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Tree.root"))


process.Tree3Mu = cms.EDAnalyzer("MiniAna2017Tree",
                              isMcLabel = cms.untracked.bool(True),
                              isAnaLabel = cms.untracked.bool(True),
                              muonLabel=cms.InputTag("looseMuons"),
                              VertexLabel=cms.InputTag("offlinePrimaryVertices"),
                              genParticleLabel=cms.InputTag("genParticles"),
                              Cand3MuLabel=cms.InputTag("ThreeMuonsVtxKalmanFit"),
                              pileupSummary = cms.InputTag("addPileupInfo"),
                              triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                              triggerSummary = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                              AlgInputTag = cms.InputTag( "gtStage2Digis" )   
)



process.Tau3MuSkim = cms.Path(process.ThreeMuonSelSeq *
                              process.Tree3Mu
                     )

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
