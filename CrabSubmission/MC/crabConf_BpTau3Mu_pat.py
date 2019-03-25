from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.General.requestName = 'BpTau3Mu_MC_pat'
config.General.workArea = 'crabProjects'
config.General.transferOutputs = True



config.JobType.pluginName = 'Analysis'
#config.JobType.psetName ='/lustrehome/fsimone/CMSSW_9_4_5/src/Analysis/MiniAna2017Tree/run_BpTau3MuSkim_AODSIM_cfg.py'
config.JobType.psetName ='/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_4/src/SkimTools/SkimTau3Mu/test/run_BpTau3MuSkim_AODSIM_cfg.py'

config.Data.inputDataset='/BpTau3Mu/fsimone-crab_crab_BpTau3Mu_13TeV_RECO-c3763c515b5d7d94a8137c090655e1bb/USER'


config.Data.splitting = 'FileBased'
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True

config.Data.ignoreLocality = True
#config.Site.whitelist   = ['T2_US_UCSD', 'T2_US_Vanderbilt', 'T2_CH_CSCS','T2_DE_DESY','T2_ES_IFCA']
config.Site.storageSite = 'T2_IT_Bari'
config.Site.whitelist   = ['T2_IT_Bari']
config.JobType.maxMemoryMB=2500
