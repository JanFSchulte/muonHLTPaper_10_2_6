import sys


from CRABClient.UserUtilities import config

config = config()

config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple.root']#, 'DQMIO.root']muonNtupleDataIOIter2.root

config.Data.unitsPerJob     = 50
#config.Data.totalUnits      = 500000
config.Data.splitting       = 'LumiBased'
#config.Data.splitting       = 'Automatic'

#config.Data.useParent       = True #!!!!
config.Data.useParent       = False #!!!!

config.Site.storageSite     = 'T3_US_FNALLPC' #'T3_US_FNALLPC'
config.JobType.numCores     = 1
config.JobType.maxMemoryMB  = 2500
config.JobType.allowUndistributedCMSSW = True

#config.Data.lumiMask ='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

#config.Data.lumiMask ='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
config.Data.lumiMask ='/afs/cern.ch/user/k/kplee/public/MuonTriggerPaper/JSON/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_Last5invfb.txt'


from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

tag = "muHLTPaper_IterL32018_RAWRECO_NoIDv2"

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
config.General.workArea   = tag
config.Data.outLFNDirBase = '/store/user/jschulte/' + tag

#config.Data.runRange = '320500-325175'

config.JobType.psetName    = 'HLTCfg2018Data_Mu.py' # 'hlt_2018.py'#'hlt_2018_CkfTrajectoryBuilder.py'#'hlt_2018.py'
#config.JobType.psetName    = 'HLTCfg2018Data_CascadeTkMu.py' # 'hlt_2018.py'#'hlt_2018_CkfTrajectoryBuilder.py'#'hlt_2018.py'
config.General.requestName = tag
#config.Data.inputDataset = '/SingleMuon/Run2018A-PromptReco-v1/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018A-PromptReco-v2/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018A-PromptReco-v3/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018B-PromptReco-v1/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018B-PromptReco-v2/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018C-PromptReco-v1/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018C-PromptReco-v2/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018C-PromptReco-v3/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2018D-PromptReco-v2/AOD'
config.Data.inputDataset = '/SingleMuon/Run2018D-ZMu-PromptReco-v2/RAW-RECO'
                           
config.Data.outputDatasetTag   = tag
