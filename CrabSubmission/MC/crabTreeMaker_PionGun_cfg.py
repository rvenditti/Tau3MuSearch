from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'TreeMaker_PionGun_2017_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_4/src/SkimTools/SkimTau3Mu/test/run_BkgStudy_PatAndTree_cfg.py'
config.Data.inputDataset = '/PionGun_Pt0to30GeV/rosma-PiGun_13TeV_MC2017_RECO-62ce2dc6bda0cbaae9cceab15421a128/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 24
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'TreeMaker_PionGun_2017_v1'

config.Site.storageSite = 'T2_IT_Bari'
