#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
if sys.version_info[:2] == (2,6):
    print "Initialize correct python version first!"
    sys.exit()

import ConfigParser
import ROOT
import math
import argparse
import re
from random import randint
from sys import maxint

# private libs
import ratio
import style
import multiplot
from rwthColors import rwth
import limitTools

import auxiliary as aux
from datasets import *

def writeSimplifiedDatacard( bkgInts, sigInts, name="tmpDataCard.txt" ):
    # assuming no observation, 50% bkg uncert, 15% signal uncert
    nBins = len(bkgInts)
    header = """
imax %d
jmax 1
kmax 2
bin %s
observation %s\n"""%(nBins," ".join(["b%d"%i for i in range(nBins)])," ".join([str(int(i)) for i in bkgInts]) )
    infos = [["bin", "process", "process", "rate", "signal lnN", "bgUnc lnN"]]
    for iBin in range(len(bkgInts)):
        bkgInfos = ["b%d"%iBin, "bg", 1, bkgInts[iBin], "-", 1.5 ]
        sigInfos = ["b%d"%iBin, "sig", 0, sigInts[iBin], 1.15, "-" ]
        bkgInfos = [ str(i) for i in bkgInfos if i is not None ]
        sigInfos = [ str(i) for i in sigInfos if i is not None ]
        infos.append(bkgInfos)
        infos.append(sigInfos)
    for l in zip(*infos):
        header += " ".join(l) + "\n"

    with open(name,"w") as f:
        f.write(header)

bkgSet = gjets+ qcd+ ttjets+ ttg+ wjets+wg_mg+zg_130+znunu
#sigSet = signal["T5Wg_1550_100"]
sigSet = signal["T5Wg_1300_100"]

bkgHist = bkgSet.getHist("tr/met")
sigHist = sigSet.getHist("tr/met")

binBoarders = [ bkgHist.GetNbinsX()+10 ]

oldR = 1e20 # maxint wouly be nicer


for bin in range( sigHist.GetNbinsX()+2, 0, -1 ):

    bkgInts = [ bkgHist.Integral(binBoarders[i+1],binBoarders[i]) for i in range(len(binBoarders)-1) ]
    sigInts = [ sigHist.Integral(binBoarders[i+1],binBoarders[i]) for i in range(len(binBoarders)-1) ]

    bkgInts.append( bkgHist.Integral(bin,binBoarders[-1]) )
    sigInts.append( sigHist.Integral(bin,binBoarders[-1]) )

    if min(bkgInts) <1e-6: continue
    if min(sigInts) <1e-6: continue


    writeSimplifiedDatacard( bkgInts, sigInts )
    r = limitTools.infosFromDatacard("tmpDataCard.txt")["exp"]
    print r
    if (r - oldR)/r>0.05: #change must me larger than 5%
        print "append"
        binBoarders.append(bin)
    oldR = r

print "final r =", r
print binBoarders
metBoarders = [ sigHist.GetBinLowEdge(i) for i in binBoarders ]
metBoarders.reverse()
print metBoarders[0:-1]


