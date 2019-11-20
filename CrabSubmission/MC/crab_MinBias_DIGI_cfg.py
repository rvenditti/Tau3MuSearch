from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'MinBias_13TeV_MC2017_kPiFilter_DIGI_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/CrabSubmission/MC/PiGun_13TeV_PU0_DIGI.py'


config.JobType.pluginName = 'Analysis'
config.Data.inputDataset ='/MinBias_TuneCP5_inelasticON_13TeV-pythia8/rosma-MinBiasPiKaonFilter_GENSIM_v1-b3ef5126d7ab537584b69881b81f2394/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 5000 #total number of events = 500000
#config.Data.totalUnits =
config.Data.inputDBS = 'phys03'

#config.Data.totalUnits = config.Data.unitsPerJo
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'MinBias_13TeV_MC2017_kPiFilter_DIGI_v1'
config.JobType.maxMemoryMB=2500
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist =True
config.Site.whitelist = ['T2_IT_Bari']
