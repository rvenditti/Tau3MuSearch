from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TreeMaker_DsSignal_2017_v4'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_4/src/SkimTools/SkimTau3Mu/test/run_Tau3MuSkimAndTree_cfg.py'
#config.Data.inputDataset = '/DsTau3Mu/fsimone-crab_crab_DsTau3Mu__13TeV_MC2016_RECO-c3763c515b5d7d94a8137c090655e1bb/USER'
config.Data.inputDataset = '/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'TreeMaker_DsSignal_2017_v4'

config.Site.storageSite = 'T2_IT_Bari'
