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

ds={
    "sp": [
        "SinglePhoton_Run2016B-PromptReco-v2_nTuple.root",
        "SinglePhoton_Run2016C-PromptReco-v2_nTuple.root",
        "SinglePhoton_Run2016D-PromptReco-v2_nTuple.root",
        "SinglePhoton_Run2016E-PromptReco-v2_nTuple.root",
        "SinglePhoton_Run2016F-PromptReco-v1_nTuple.root",
        "SinglePhoton_Run2016G-PromptReco-v1_nTuple.root",
        ],
    "jh": [
        "JetHT_Run2016B-PromptReco-v2_nTuple.root",
        "JetHT_Run2016C-PromptReco-v2_nTuple.root",
        "JetHT_Run2016D-PromptReco-v2_nTuple.root",
        "JetHT_Run2016E-PromptReco-v2_nTuple.root",
        "JetHT_Run2016F-PromptReco-v1_nTuple.root",
        "JetHT_Run2016G-PromptReco-v1_nTuple.root",
        ],
    "gjet": [
        "GJets_HT-40To100_nTuple.root",
        "GJets_HT-100To200_nTuple.root",
        "GJets_HT-200To400_nTuple.root",
        "GJets_HT-400To600_nTuple.root",
        "GJets_HT-600ToInf_nTuple.root",
        ],
    "qcd": [
        "QCD_HT200to300_nTuple.root",
        "QCD_HT300to500_nTuple.root",
        "QCD_HT500to700_nTuple.root",
        "QCD_HT700to1000_nTuple.root",
        "QCD_HT1000to1500_nTuple.root",
        "QCD_HT1500to2000_nTuple.root",
        "QCD_HT2000toInf_nTuple.root",
        ],
    "w": [
        "WJetsToLNu_HT-100To200_nTuple.root",
        "WJetsToLNu_HT-1200To2500_nTuple.root",
        "WJetsToLNu_HT-200To400_nTuple.root",
        "WJetsToLNu_HT-2500ToInf_nTuple.root",
        "WJetsToLNu_HT-400To600_nTuple.root",
        "WJetsToLNu_HT-600To800_nTuple.root",
        "WJetsToLNu_HT-800To1200_nTuple.root",
        ],
    "znunu": [
        "ZJetsToNuNu_HT-100To200_nTuple.root",
        "ZJetsToNuNu_HT-200To400_nTuple.root",
        "ZJetsToNuNu_HT-400To600_nTuple.root",
        "ZJetsToNuNu_HT-600To800_nTuple.root",
        "ZJetsToNuNu_HT-800To1200_nTuple.root",
        "ZJetsToNuNu_HT-1200To2500_nTuple.root",
        "ZJetsToNuNu_HT-2500ToInf_nTuple.root",
        ],
    "tt": [
        "TTJets_nTuple.root",
        "TTJets_HT-0to600_nTuple.root",
        "TTJets_HT-600to800_ext_nTuple.root",
        "TTJets_HT-800to1200_ext_nTuple.root",
        "TTJets_HT-1200to2500_ext_nTuple.root",
        "TTJets_HT-2500toInf_nTuple.root",
    ],
    "ttg": ["TTGJets_nTuple.root"],
    "tg": ["TGJets_amcatnlo_madspin_nTuple.root"],
    "zg": [
        "ZNuNuGJets_MonoPhoton_PtG-40to130_nTuple.root",
        "ZNuNuGJets_MonoPhoton_PtG-130_nTuple.root",
        ],
    "wg": [
        "WGJets_MonoPhoton_PtG-130_nTuple.root",
        "WGToLNuG-madgraphMLM_nTuple.root",
        "WGToLNuG-madgraphMLM_PtG-0to130_nTuple.root",
        "WGToLNuG-amcatnloFXFX_nTuple.root",
        ],
}
dir = "/user/kiesel/nTuples/v16/"


#############################################
# Select datasets to process
#############################################

# compile only
run.run()

parser = argparse.ArgumentParser()
parser.add_argument('datasets', nargs='+', default=["all"], help="all "+' '.join(ds.keys()))
parser.add_argument('--condor', action="store_true")
args = parser.parse_args()

if args.datasets == ["all"]:
    toProcess = [x for sublist in ds.values() for x in sublist]
else:
    toProcess = ds[args.datasets[0]]
    for n in args.datasets[1:]:
        toProcess += ds[n]

print toProcess

if args.condor:
    for x in toProcess:
        with open("submit","w") as f:
            f.write("""
Universe   = vanilla
Executable = run.sh
Arguments  = /net/scratch_cms1b1/cms/user/kiesel/{0}
Log        = logs/{0}.log
Output     = logs/{0}.out
Error      = logs/{0}.error
Queue
""".format(x))
        subprocess.call(["condor_submit", "submit"])

else: # local processing


    files = [dir+x for x in toProcess]
    files.sort(key=os.path.getsize, reverse=True)

    p = multiprocessing.Pool()
    p.map(run.run, files)
