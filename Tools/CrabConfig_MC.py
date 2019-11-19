import sys


from CRABClient.UserUtilities import config

config = config()

config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple.root']#, 'DQMIO.root']muonNtupleDataIOIter2.root

config.Data.unitsPerJob     = 10000
config.Data.totalUnits      = -1
config.Data.splitting       = 'EventAwareLumiBased'

config.Data.useParent       = True #!!!!
#config.Data.useParent       = False #!!!!

config.Site.storageSite     = 'T3_US_FNALLPC'
config.JobType.numCores     = 1
config.JobType.maxMemoryMB  = 2500
config.JobType.allowUndistributedCMSSW = True
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

tag = "muonHLTpaper_DYJets_Fall102XFlat_NoID"
#tag = "muonHLTpaper_WJets_Fall102XFlat_CascadeTkMu_v5"
#tag = "muonHLTpaper_DYJets_Fall102XFlat_CascadeTkMu_v5"

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
config.General.workArea   = tag
config.Data.outLFNDirBase = '/store/user/jschulte/' + tag



#config.JobType.psetName    = 'HLTCfg2018MC_CascadeTkMu.py' # 
config.JobType.psetName    = 'HLTCfg2018MC_Mu.py' # 
config.General.requestName = tag
config.Data.inputDataset ='/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18DR-FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/AODSIM'
#config.Data.inputDataset ='/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISpring18DR-NZSPU40to70_100X_upgrade2018_realistic_v10-v1/AODSIM'

config.Data.outputDatasetTag   = tag
