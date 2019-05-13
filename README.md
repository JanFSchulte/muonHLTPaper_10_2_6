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


# Ntuple location 
Some ntuples are already produced (using 2018 data D) and the location is here in lpceos
```
/eos/uscms/store/user/bmahakud/muHLTPaper_SingleMu_PromptReco2018D_v2/SingleMuon/muHLTPaper_SingleMu_PromptReco2018D_v2/190511_164115/0000/
```

# Getting the L3/L1 efficiency stack plots for paper
```
cd src/MuonHLTNtuples/Tools/Macros/IterL3
```

Now run the macro to produce the efficiency plots 
```
root -l -b -q 'readNtuplesPrefilter_IterL3ForStack.C("/eos/uscms/store/user/bmahakud/muHLTPaper_SingleMu_PromptReco2018D_v2/SingleMuon/muHLTPaper_SingleMu_PromptReco2018D_v2/190511_164115/0000/muonNtuple_*.root","PromptRecoD")'
```
This will end up produing a file PromptRecoD_IterL3preFilter.root that will contain all the required efficiency histograms.

#Produce the plot
```
cd /src/MuonHLTNtuples/Tools/Macros/IterL3/Plots
```
Run the following macro to produce the stack plots

```
root -l PlotStack.C
```
It should produce plots  similar to what you see in the following location

```
http://bmahakud.web.cern.ch/bmahakud/MuonHLT/muonHLTPaperPlots/IterL3/
```














