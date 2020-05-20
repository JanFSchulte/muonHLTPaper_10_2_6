jdlTemplate='''
universe = vanilla
Executable = runPtRes.sh
Should_Transfer_Files = YES
whenToTransferOutput = ON_EXIT
Transfer_Input_Files = runInterpretation.sh
Output = zPrimeCombine_$(Cluster)_$(Process).stdout
Error = zPrimeCombine_$(Cluster)_$(Process).stderr
Log = zPrimeCombine_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
'''

scriptTemplate='''
#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a bash script, use .sh instead of .csh
### for case 1. EOS have the following line, otherwise remove this line in case 2.
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
root
for filename in results_Run2CI_CMSWEEK/*.root; do
        fileBase=$(basename $filename)
        xrdcp -f results_Run2CI_CMSWEEK/${fileBase} root://cmseos.fnal.gov//store/user/jschulte/limits/results_Run2CI_CMSWEEK/${fileBase}
done
### remove the output file if you don't want it automatically transferred when the job ends
cd ${_CONDOR_SCRATCH_DIR}
rm *

'''
