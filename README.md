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
If you endup with a  succesful cmsRun crabjobs could be submitted using CrabConfig_DataD_v2.py 









