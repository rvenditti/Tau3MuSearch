from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'MinBiasPiKaonFilter_GENSIM_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/lustrehome/venditti/Tau3MU_29072019/CMSSW_9_4_4/src/MinimumBias_Filter.py'



config.Data.inputDataset = '/MinBias_TuneCP5_inelasticON_13TeV-pythia8/RunIIFall17GS-93X_mc2017_realistic_v3-v1/GEN-SIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'MinBiasPiKaonFilter_GENSIM_v1'
#config.Data.ignoreLocality = True
config.Site.storageSite = 'T2_IT_Bari'
config.JobType.allowUndistributedCMSSW = True
