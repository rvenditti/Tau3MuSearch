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
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_RealisticBS_25ns_13TeV2016_v1_mc'
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_69.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_7.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_70.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_71.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_72.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_73.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_74.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v3/190213_150844/0000/file_AODSIM_75.root',

#'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_1.root',  
#'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_10.root',                            
#'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_11.root',                             
#'file:/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_5/src/file_AODSIM_SimVertex.root',

#       'file:/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_5/src/MiniAna2017/MiniAna2017Tree/file_AODSIM.root'
        #'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017/190129_213256/0000/file_AODSIM_1.root',
        #'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017/190129_213256/0000/file_AODSIM_2.root'
        #'/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
        #'/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/90000/FEFE887F-E944-E811-ADCB-A0369FE2C14A.root'
        #'/store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/ECAA9E63-C5FD-E711-9C1E-008CFAE451DC.root'
        #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/60000/D699CEF5-E1FA-E711-A7EC-02163E013399.root'
        #'root://xrootd-cms.infn.it//store/data/Run2017F/DoubleMuonLowMass/MINIAOD/17Nov2017-v1/00000/129CB4BC-17FD-E711-8D49-FA163E838299.root'

    )
)

process.demo = cms.EDAnalyzer("DsTau3MuAna",
                              #isMcLabel = cms.untracked.bool(True),
                              muonLabel=cms.InputTag("looseMuons"),
                              VertexLabel=cms.InputTag("offlinePrimaryVertices"),
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

"""
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_1.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_10.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_11.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_12.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_13.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_14.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_15.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_16.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_17.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_18.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_19.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_2.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_20.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_21.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_22.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_23.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_24.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_25.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_26.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_27.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_28.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_29.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_3.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_30.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_31.root',
'file:/lustre/cms/store/user/rosma/DsTau3Mu/SkimTau3Mu_DsSignal_2017_v2/190130_142447/0000/file_AODSIM_32.root',
"""
