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

#process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.GlobalTag.globaltag = '94X_dataRun2_v11' #data2017 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://xrootd-cms.infn.it//store/data/Run2017B/DoubleMuonLowMass/AOD/23Jun2017-v1/90000/FC4BB0C3-E358-E711-9EFE-0025904B739A.root',
        'root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/AOD/09May2018-v1/80000/FEDF5D97-BEB0-E811-95BF-0CC47AD98B94.root'
        #root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/AOD/09May2018-v1/80000/AECC4C56-BAB0-E811-B92A-008CFA1979AC.root'
        #"file:/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_4/src/CrabSubmission/MC/PiGun_RECO.root"
        #"file:/lustre/cms/store/user/rosma/PionGun_Pt0to30GeV/PiGun_13TeV_MC2017_RECO/190313_143541/0000/PiGun_RECO_979.root"
        #'file:/lustre/cms/store/user/fsimone/DsTau3Mu/crab_crab_DsTau3Mu__13TeV_MC2016_RECO/190120_140919/0001/custom_DsTau3Mu_13TeV_RECO_crab350_1010.root',


    )
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("TreeData.root"))


process.TreeMakerBkg = cms.EDAnalyzer("MiniAna2017Tree",
                              isMcLabel = cms.untracked.bool(False),
                              isAnaLabel = cms.untracked.bool(True),
                              muonLabel=cms.InputTag("looseMuons"),
                              VertexLabel=cms.InputTag("offlinePrimaryVertices"),
                              genParticleLabel=cms.InputTag("genParticles"),
                              pileupSummary = cms.InputTag("addPileupInfo"),
                              Cand3MuLabel=cms.InputTag("ThreeMuonsVtxKalmanFit"),
                              triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
                              triggerSummary = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
                              AlgInputTag = cms.InputTag( "gtStage2Digis" )
)




process.Tau3MuSkim = cms.Path(process.ThreeMuonSelSeq*
                              process.TreeMakerBkg
                     )


#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("file_AODSIM_test_ForTreeMaker.root"),
#                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('PatMuonSequence')),
#                               outputCommands = cms.untracked.vstring(
#        'drop *',
#        'keep *_offlinePrimaryVertices_*_*',
#        'keep *_generator_*_*',
#        'keep *_offlineBeamSpot_*_*',
#        'keep recoMuon_muons_*_*',
#        'keep *_TriggerResults_*_*',
#        'keep *_gtStage2Digis_*_*',
#        'keep *_gmtStage2Digis_Muon_*',
#        'keep *_scalersRawToDigi_*_*',
#        'keep *_Trigger_*_*',
#        'keep *_addPileupInfo_*_*',
#        'keep *_genParticles_*_*',
#         'keep recoVertex_offlinePrimaryVertices_*_*',
#       'keep *_generalTracks_*_RECO',
#       'keep TrackExtra_generalTracks_*_RECO',
#         'keep *_globalMuons_*_RECO',
 #        'keep *_standAloneMuons_*_RECO',
  #       'keep PSimHits_g4SimHits_MuonCSCHits_SIM',
   #      'keep PSimHits_g4SimHits_MuonDTHits_SIM',
    #     'keep PSimHits_g4SimHits_MuonRPCHits_SIM',
     #    'keep SimVertexs_g4SimHits__SIM',
      #   'keep *_csc2DRecHits_*_*',
       #  'keep *_dt1DRecHits_*_*',
       #  'keep *_rpcRecHits_*_*',
       #  'keep SimVertexs_g4SimHits_*_*',
       #  'drop floats_generalTracks_MVAValues_RECO',        
       #  'drop recoTrack_globalMuons_*_RECO',
       #  'drop recoTrack_standAloneMuons_*_RECO',
        #'keep recoTracks_generalTracks_*_*',

        #)
#)

#type_label_instance_process
#process.outpath = cms.EndPath(process.out) 




