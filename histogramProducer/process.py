#!/usr/bin/env python2
import os
import multiprocessing

# own libray
import run

ds={
    "sp": [
        "SinglePhoton_Run2015C_25ns-05Oct2015-v1_nTuple.root",
        "SinglePhoton_Run2015D-05Oct2015-v1_nTuple.root",
        "SinglePhoton_Run2015D-PromptReco-v4_nTuple.root",
        ],
    "jh": [
        "JetHT_Run2015C_25ns-05Oct2015-v1_nTuple.root",
        "JetHT_Run2015D-05Oct2015-v1_nTuple.root",
        "JetHT_Run2015D-PromptReco-v4_nTuple.root",
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
    "ewk": [
        "DYJetsToLL_M-50_nTuple.root",
        "TTJets_nTuple.root",
        ],
    "gv": [
        "TTGJets_nTuple.root",
        "ZNuNuGJets_MonoPhoton_PtG-130_nTuple.root",
        "WGJets_MonoPhoton_PtG-130_nTuple.root",
        "WGToLNuG-madgraphMLM_nTuple.root",
        ]
}
dir = "/user/kiesel/nTuples/v09/"


#############################################
# Select datasets to process
#############################################

toProcess = ds["gjet"]+ds["qcd"]
toProcess = ds["znunu"]+["ZNuNuGJets_MonoPhoton_PtG-130_nTuple.root"]
#toProcess = ["TTJets_nTuple.root","TTGJets_nTuple.root"]

#############################################
#############################################

files = [dir+x for x in toProcess]
files.sort( key=os.path.getsize, reverse=True )

# compile only
run.run()

p = multiprocessing.Pool()
p.map( run.run, files )
