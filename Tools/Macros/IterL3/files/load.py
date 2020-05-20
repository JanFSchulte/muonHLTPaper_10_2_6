import os, subprocess


command = "gfal-ls gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/muHLTPaper_SingleMu_PromptReco2018D_v2/SingleMuon/muHLTPaper_SingleMu_PromptReco2018D_v2/190516_181455/0000"

child = os.popen(command)
data = child.read()
err = child.close()

files = data.splitlines()
print files


for file in files:
	if 'DQM' in file: continue
	print file
	print ['gfal-copy','-f','gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/muHLTPaper_SingleMu_PromptReco2018D_v2/SingleMuon/muHLTPaper_SingleMu_PromptReco2018D_v2/190516_181455/0000/%s'%file, "file:///%s"%file]
	subprocess.call(['gfal-copy -f gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/muHLTPaper_SingleMu_PromptReco2018D_v2/SingleMuon/muHLTPaper_SingleMu_PromptReco2018D_v2/190516_181455/0000/%s file:///%s'%(file,file)])
	#subprocess.call(['gfal-copy -f gsiftp://cms-gridftp.rcac.purdue.edu/store/user/jschulte/muHLTPaper_SingleMu_PromptReco2018D_v2/SingleMuon/muHLTPaper_SingleMu_PromptReco2018D_v2/190516_181455/0000/%s file:///%s'%(file,file)],shell=True,stdout=open(os.devnull, 'wb'))
