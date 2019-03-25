import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


#process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v6'
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_RealisticBS_25ns_13TeV2016_v1_mc'
#process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

#       'file:/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_5/src/MiniAna2017/MiniAna2017Tree/file_AODSIM.root'
        #'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017/190129_213256/0000/file_AODSIM_1.root',
        #'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017/190129_213256/0000/file_AODSIM_2.root'
        '/store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/ECAA9E63-C5FD-E711-9C1E-008CFAE451DC.root'
        #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/D699CEF5-E1FA-E711-A7EC-02163E013399.root'
        #'root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/MINIAOD/17Nov2017-v1/00000/129CB4BC-17FD-E711-8D49-FA163E838299.root'

    )
)

process.demo = cms.EDAnalyzer("DsTau3MuAna",
                              #isMcLabel = cms.untracked.bool(True),
                              muonLabel=cms.InputTag("slimmedMuons"),
                              VertexLabel=cms.InputTag("slimmedOfflinePrimaryVertices"),
                              genParticleLabel=cms.InputTag("genParticles"), 
                              Cand3MuLabel=cms.InputTag("ThreeMuonsCand")
)


#process.options = cms.untracked.PSet(
#  SkipEvent = cms.untracked.vstring( "Error: uninitialized ProxyBase used" ),
  #IgnoreCompletely = cms.untracked.vstring( "ProductNotFound" )
#)

process.TFileService = cms.Service("TFileService",
fileName = cms.string("histo.root")
)

process.p = cms.Path(process.demo)
