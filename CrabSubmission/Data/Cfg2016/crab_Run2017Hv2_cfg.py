from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'SkimTau3Mu_DoubleMuonLowMass_Run2016H_v2_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'

config.JobType.psetName = '/lustrehome/venditti/Tau3Mu2016/CMSSW_8_0_21/src/SkimTools/SkimTau3Mu/test/run_Data2016_PatAndTree_cfg.py'

config.Data.inputDataset = '/DoubleMuonLowMass/Run2016H-PromptReco-v2/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 100
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Era/Prompt/Cert_280919-284044_13TeV_PromptReco_Collisions16_JSON_eraH.txt'
#config.Data.runRange = '193093-193999' # '193093-194075'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'SkimTau3Mu_DoubleMuonLowMass_Run2016Hv2_v1'
config.JobType.allowUndistributedCMSSW = True 
config.Site.storageSite = 'T2_IT_Bari'


