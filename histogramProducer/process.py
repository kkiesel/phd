#!/usr/bin/env python2
import os
import multiprocessing
import glob

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# own libray
import run

ds={
    "sp": [
        "SinglePhoton_Run2015C_25ns-16Dec2015-v1_nTuple.root",
        "SinglePhoton_Run2015D-16Dec2015-v1_nTuple.root",
        ],
    "jh": [
        "JetHT_Run2015C_25ns-16Dec2015-v1_nTuple.root",
        "JetHT_Run2015D-16Dec2015-v1_nTuple.root",
        ],
    "gjet": [
        "GJets_HT-40To100_nTuple.root",
        "GJets_HT-100To200_nTuple.root",
        "GJets_HT-200To400_nTuple.root",
        "GJets_HT-400To600_nTuple.root",
        "GJets_HT-600ToInf_nTuple.root",
        ],
    "qcd": [
        "QCD_HT100to200_nTuple.root",
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
        "ZJetsToNuNu_HT-600ToInf_nTuple.root",
        ],
    "tt": ["TTJets_nTuple.root"],
    "ttg": ["TTGJets_nTuple.root"],
    "tg": ["TGJets_amcatnlo_madspin_nTuple.root"],
    "zg": ["ZNuNuGJets_MonoPhoton_PtG-130_nTuple.root",
        "ZGTo2LG_nTuple.root",
        "ZGTo2LGmod_nTuple.root"],
    "wg": ["WGJets_MonoPhoton_PtG-130_nTuple.root",
           "WGToLNuG-madgraphMLM_nTuple.root"],
    "signals": ["T5Wg_1550_1500.root","T5Wg_1550_100.root"],
}
dir = "/user/kiesel/nTuples/v12/"


#############################################
# Select datasets to process
#############################################

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('datasets', nargs='+', default=["all"], help="all "+' '.join(ds.keys()))
args = parser.parse_args()

if args.datasets == ["all"]:
    toProcess = [x for sublist in ds.values() for x in sublist]
else:
    toProcess = ds[args.datasets[0]]
    for n in args.datasets[1:]:
        toProcess += ds[n]

print toProcess
#############################################
#############################################

files = [dir+x for x in toProcess]
files.sort(key=os.path.getsize, reverse=True)

# adding signal scan
#files = [ f for f in glob.glob( dir+"T5Wg_*.root") ]

# compile only
run.run()

p = multiprocessing.Pool()
p.map(run.run, files)
