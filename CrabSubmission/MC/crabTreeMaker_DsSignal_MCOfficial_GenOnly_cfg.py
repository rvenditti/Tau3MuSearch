from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TreeMaker_DsSignalOfficial_2017_v4_GenOnly'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/SkimTools/SkimTau3Mu/test/run_Tau3MuSkimAndTree_NoSelections_cfg.py'
config.Data.inputDataset = '/DsToTau_To3Mu_MuFilter_TuneCUEP8M1_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'TreeMaker_DsSignalOfficial_2017_v4_GenOnly'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
