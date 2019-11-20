from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DsPhiMMPi13TeV_MC2017_LowMuPt_DIGI'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/CrabSubmission/MC/PiGun_13TeV_DIGI.py'


config.JobType.pluginName = 'Analysis'
config.Data.inputDataset ='/DsPhiMMPi13TeV_Pt0p5To2p7_GENSIM/rosma-DsPhiPi13TeV_Pt0To2p7_GENSIM-5a4f83e2e84c2b975c29622d0650561b/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 5000 #total number of events = 500000
#config.Data.totalUnits =
config.Data.inputDBS = 'phys03'

#config.Data.totalUnits = config.Data.unitsPerJo
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'DsPhiMMPi13TeV_MC2017_LowMuPt_DIGI'
config.JobType.maxMemoryMB=2500
config.JobType.allowUndistributedCMSSW = True
config.Site.storageSite = 'T2_IT_Bari'
config.Site.ignoreGlobalBlacklist =True
config.Site.whitelist = ['T2_IT_Bari']
