from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DsToPhiPi_ToMuMu_MuFilter_SkimPhiPi_MC_2017_v7'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/Tau3Mu2017_dev/CMSSW_9_4_4/src/DsPhiPiTreeMaker/DsPhiPiTreeMaker/test/run_DsPhiPiSkimAndTree_cfg.py'

config.Data.inputDataset = '/DsToPhiPi_ToMuMu_MuFilter_TuneCUEP8M1_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'DsToPhiPi_ToMuMu_MuFilter_SkimPhiPi_MC_2017_v7'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
#config.Site.ignoreGlobalBlacklist =True
#config.Site.whitelist = ['T2_IT_Bari']

