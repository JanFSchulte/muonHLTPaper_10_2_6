# Simplified setup for muonHLTPaper for plots using CMSSW_10_2_6 (2018 data and MC)

This setup is tested in cmslpc. Should work in lxplus too. 

```
cmsrel CMSSW_10_2_6
cd src
cmsenv
git cms-addpkg HLTrigger/Configuration
git clone https://github.com/bmahakud/muonHLTPaper_10_2_6 MuonHLTNtuples
scram b -j 8
voms-proxy-init -voms cms
cd MuonHLTNtuples/Tools
cmsRun HLTCfg2018Data_Mu.py
```
The above config file is generated using the command 

```
hltGetConfiguration orcoff:/cdaq/physics/Run2018/2e34/v3.6.1/HLT/V2 \
--globaltag 101X_dataRun2_HLT_v7 \
--path HLTriggerFirstPath,\
HLT_IsoMu24_v*,\
HLT_IsoMu27_v*,\
HLT_Mu50_v*,\
HLT_OldMu100_v*,\
HLT_TkMu100_v*,\
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*,\
HLTriggerFinalPath,\
HLTAnalyzerEndpath \
--input /store/data/Run2018D/SingleMuon/RAW/v1/000/325/170/00000/13A7C0D8-97DF-324A-8741-F9BA4BF0C3B8.root \
--process MYHLT --full --offline \
--l1-emulator uGT \
--l1 L1Menu_Collisions2018_v2_1_0-d1_xml \
--prescale none --max-events 100 --output none > HLTCfg2018Data_Mu.py
```



If you endup with a  succesful cmsRun crabjobs could be submitted using CrabConfig_DataD_v2.py 












