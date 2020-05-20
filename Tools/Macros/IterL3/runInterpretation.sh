#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
### for case 1. EOS have the following line, otherwise remove this line in case 2.
xrdcp -s root://cmseos.fnal.gov//store/user/jschulte/CMSSW10213.tgz .
tar -xf CMSSW10213.tgz
rm CMSSW10213.tgz
export SCRAM_ARCH=slc6_amd64_gcc530
cd CMSSW_10_2_13/src/ZPrimeCombine
scramv1 b ProjectRename
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
python runInterpretation.py -c Run2CI -t CMSWEEK -r --CI -m ${1} -L ${2} -n ${3} # runs the actual calculations
for filename in results_Run2CI_CMSWEEK/*.root; do
	fileBase=$(basename $filename)
	xrdcp -f results_Run2CI_CMSWEEK/${fileBase} root://cmseos.fnal.gov//store/user/jschulte/limits/results_Run2CI_CMSWEEK/${fileBase}
done
### remove the output file if you don't want it automatically transferred when the job ends
rm -r dataCards_*
rm -r results_*
cd ${_CONDOR_SCRATCH_DIR}
rm -rf CMSSW_10_2_13
