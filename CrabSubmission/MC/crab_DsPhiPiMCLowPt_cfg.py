from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'SkimPhiPi_MC_LowPt_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/DsPhiPiTreeMaker/DsPhiPiTreeMaker/test/run_DsPhiPiSkimAndTreeLowPt_cfg.py'


config.Data.inputDataset = '/DsPhiMMPi13TeV_Pt0p5To2p7_GENSIM/rosma-DsPhiMMPi13TeVPt0p5To2p7_MC2017_RECO-475fbd8e0e9ef3c18fb1b10f3da8b578/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimPhiPi_MC_LowPt_v2'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist =True
config.Site.whitelist = ['T2_IT_Bari']

