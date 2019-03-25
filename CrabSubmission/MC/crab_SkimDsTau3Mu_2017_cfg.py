from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'SkimTau3Mu_DsSignal_2017_v4'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_4/src/SkimTools/SkimTau3Mu/test/run_Tau3MuSkim_AODSIM_cfg.py'
config.Data.inputDataset = '/DsTau3Mu/fsimone-crab_crab_DsTau3Mu__13TeV_MC2016_RECO-c3763c515b5d7d94a8137c090655e1bb/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimTau3Mu_DsSignal_2017_v4'
#config.Data.ignoreLocality = True
config.Site.storageSite = 'T2_IT_Bari'
