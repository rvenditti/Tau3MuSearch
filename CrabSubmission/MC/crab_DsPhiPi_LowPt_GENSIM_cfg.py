from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DsPhiPi13TeV_Pt0To2p7_GENSIM'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/CrabSubmission/MC/DsPhiPiLowPt_GEN_SIM_cfg.py'

config.Data.outputPrimaryDataset = 'DsPhiMMPi13TeV_Pt0p5To2p7_GENSIM'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1000
NJOBS = 10000 #total number of events = 500000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.inputDBS = 'phys03'

config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'DsPhiPi13TeV_Pt0To2p7_GENSIM'
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
