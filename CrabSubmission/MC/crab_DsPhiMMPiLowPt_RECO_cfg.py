from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DsPhiMMPi13TeVPt0p5To2p7_MC2017_RECO'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/CrabSubmission/MC/PiGun_13TeV_RECO.py'



config.JobType.pluginName = 'Analysis'
config.Data.inputDataset ='/DsPhiMMPi13TeV_Pt0p5To2p7_GENSIM/rosma-DsPhiMMPi13TeV_MC2017_LowMuPt_DIGI-071c106241af289b3c5038b86ba29bf4/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
#NJOBS = 5000 #total number of events = 500000
#config.Data.totalUnits =
config.Data.inputDBS = 'phys03'

#config.Data.totalUnits = config.Data.unitsPerJo
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'DsPhiMMPi13TeVPt0p5To2p7_MC2017_RECO'
config.JobType.maxMemoryMB=2500

config.Site.storageSite = 'T2_IT_Bari'
config.JobType.allowUndistributedCMSSW = True
config.Site.ignoreGlobalBlacklist =True
config.Site.whitelist = ['T2_IT_Bari']
