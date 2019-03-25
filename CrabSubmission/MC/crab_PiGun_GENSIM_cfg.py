from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'KaonGun_13TeV_MC2017_GENSIM'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = '/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_4/src/SingleKPt30Eta2p5_cfi_GEN_SIM.py'

config.Data.outputPrimaryDataset = 'KaonGun_Pt0to30GeV'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 40 #we expect 50 events per file
NJOBS = 5000 #total number of events = 500000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.inputDBS = 'phys03'

config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'KaonGun_Pt0to30GeV'

config.Site.storageSite = 'T2_IT_Bari'
