import sys,os

#path="/eos/uscms/store/user/bmahakud/ProductionHLTAN_LPC_IterL3HighStat/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/ProductionHLTAN_LPC_IterL3HighStat/181130_193653/0000/muonNtupleData_*.root"

#path="/eos/uscms/store/user/bmahakud/SingleMu_PromptReco2018A_v1_Crv2/SingleMuon/SingleMu_PromptReco2018A_v1_Crv2/190203_010755/0000/muonNtuple_*.root"
#path="/eos/uscms/store/user/bmahakud/SingleMu_PromptReco2018B_v1/SingleMuon/SingleMu_PromptReco2018B_v1/190203_010528/0000/muonNtuple_*.root"

#path="/uscms_data/d3/bmahakud/PostDocWork/MuonHLT/SetupForHLTNote/GlobalSetupForPaper/CMSSW_10_2_6/src/MuonHLTNtuples/Tools/SingleMuon_debug.root"


path="/eos/uscms/store/user/bmahakud/muHLTPaper_SingleMu_PromptReco2018D_v2/SingleMuon/muHLTPaper_SingleMu_PromptReco2018D_v2/190511_164115/0000/muonNtuple_*.root"

sample="PromptRecoD"

command4 ="root -l -b -q 'readNtuplesPrefilter_IterL3ForStack.C(\"%s\",\"%s\")'"%(str(path),str(sample))
os.system(command4)
















