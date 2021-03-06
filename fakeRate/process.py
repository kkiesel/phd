#!/usr/bin/env python2
import argparse
import os
import multiprocessing
import glob
import subprocess

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# own libray
import run

periods = [
        "SingleElectron_Run2016B-23Sep2016-v3_nTuple.root",
        "SingleElectron_Run2016C-23Sep2016-v1_nTuple.root",
        "SingleElectron_Run2016D-23Sep2016-v1_nTuple.root",
        "SingleElectron_Run2016E-23Sep2016-v1_nTuple.root",
        "SingleElectron_Run2016F-23Sep2016-v1_nTuple.root",
        "SingleElectron_Run2016G-23Sep2016-v1_nTuple.root",
        "SingleElectron_Run2016H-PromptReco-v2_nTuple.root",
        "SingleElectron_Run2016H-PromptReco-v3_nTuple.root",
        "DYJetsToLL_M-50_nTuple.root",
        "DYJetsToLL_M-50_ext_nTuple.root",
        "DYJetsToLL_M-50-herwigpp_30M_nTuple.root",
        ]
dir = "/net/scratch_cms1b1/cms/user/kiesel/v18/"


#############################################
# Select datasets to process
#############################################

# compile only
run.run()

parser = argparse.ArgumentParser()
parser.add_argument('--condor', action="store_true")
args = parser.parse_args()

if args.condor:
    for x in periods:
        with open("submit","w") as f:
            f.write("""
Universe   = vanilla
Executable = run.sh
Arguments  = /net/scratch_cms1b1/cms/user/kiesel/v18/{0}
Log        = logs/{0}.log
Output     = logs/{0}.out
Error      = logs/{0}.error
Queue
""".format(x))
        subprocess.call(["condor_submit", "submit"])

else: # local processing

    files = [dir+x for x in periods]
    files.sort(key=os.path.getsize, reverse=True)

    p = multiprocessing.Pool()
    p.map(run.run, files)
