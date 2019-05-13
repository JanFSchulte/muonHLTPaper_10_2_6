import sys


from CRABClient.UserUtilities import config

config = config()

config.JobType.pluginName   = 'Analysis'
config.JobType.outputFiles  = ['muonNtuple.root']#, 'DQMIO.root']muonNtupleDataIOIter2.root

config.Data.unitsPerJob     = 5000
config.Data.totalUnits      = 1000000
config.Data.splitting       = 'EventAwareLumiBased'

config.Data.useParent       = True #!!!!
#config.Data.useParent       = False #!!!!

config.Site.storageSite     = 'T3_US_FNALLPC' #'T3_US_FNALLPC'
config.JobType.numCores     = 1
config.JobType.maxMemoryMB  = 2500

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

tag = "muonHLTpaper_DYJets_Fall102XFlat"

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
config.General.workArea   = tag
config.Data.outLFNDirBase = '/store/user/bmahakud/' + tag



config.JobType.psetName    = 'HLTCfg2018MC_Mu.py' # 
config.General.requestName = tag
config.Data.inputDataset ='/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18DR-FlatPU28to62NZS_102X_upgrade2018_realistic_v15-v1/AODSIM'

config.Data.outputDatasetTag   = tag
