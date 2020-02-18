from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'SkimTau3Mu_DoubleMuonLowMass_Run2016G_missLumi_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = '/lustrehome/venditti/Tau3Mu2016/CMSSW_8_0_21/src/SkimTools/SkimTau3Mu/test/run_Data2016_PatAndTree_cfg.py'

config.Data.inputDataset = '/DoubleMuonLowMass/Run2016G-23Sep2016-v1/AOD'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'LumiBased'
config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 35
config.Data.lumiMask = '/lustrehome/venditti/Tau3Mu2016/CMSSW_8_0_21/src/CrabSubmission/Data/SkimTau3Mu/crab_projects/crab_SkimTau3Mu_DoubleMuonLowMass_Run2016G_v1/results/notFinishedLumis.json'
#config.Data.runRange = '193093-193999' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimTau3Mu_DoubleMuonLowMass_Run2016G_AOD_missLumi__v1'
#config.JobType.allowUndistributedCMSSW = True 
config.Site.storageSite = 'T2_IT_Bari'


config.JobType.allowUndistributedCMSSW = True
