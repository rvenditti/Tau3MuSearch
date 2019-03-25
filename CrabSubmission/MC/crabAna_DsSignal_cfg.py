from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'VtxAna_DsSignal_2017_v4'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/TestMiniAOD2017/CMSSW_9_4_4/src/DsTau3MuAna/DsTau3MuAna/python/ConfFile_cfg.py'
config.Data.inputDataset = '/DsTau3Mu/rosma-SkimTau3Mu_DsSignal_2017_v3-60e6045cf7c9f7395b636516d3893df1/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'VtxAna_DsSignal_2017_v4'

config.Site.storageSite = 'T2_IT_Bari'
