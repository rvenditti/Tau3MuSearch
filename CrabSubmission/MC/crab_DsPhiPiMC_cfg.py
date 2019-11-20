from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'SkimPhiPi_MC_Nim'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/DsPhiPiTreeMaker/DsPhiPiTreeMaker/test/run_DsPhiPiSkimAndTree_cfg.py'


config.Data.inputDataset = '/DsTau3Mu/nperera-CRAB3_MC2016_DsTau3Mu_13TeV_RECO2-4c560ee1ffc2a2467def179a8eb13ba9/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimPhiPi_MC_Nim'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist =True
config.Site.whitelist = ['T2_IT_Bari']

