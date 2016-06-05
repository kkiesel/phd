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

def simplifiedSig(s, b):
    if b<1e-6: b = 2
    return 0.5 * s / math.sqrt(s+b+(0.5*b)**2)

def optimizeBinnigs2d(bkgSet, name, sigPoint):
    style.style2d()
    bkgHist = bkgSet.getHist("tr/metRaw_vs_emht")
    sigHist = signal[sigPoint].getHist("tr/metRaw_vs_emht")

    xmax = 800
    ymax = 4000
    xmaxBin = bkgHist.GetXaxis().FindBin(xmax)
    ymaxBin = bkgHist.GetYaxis().FindBin(ymax)

    signifHist = sigHist.Clone(aux.randomName())
    for xbin, ybin in aux.loopH(sigHist):
        if xbin > xmaxBin or ybin > ymaxBin: continue
        b = bkgHist.Integral(xbin, -1, ybin, -1)
        s = sigHist.Integral(xbin, -1, ybin, -1)
        signif = simplifiedSig(s, b)
        signifHist.SetBinContent(xbin, ybin, signif)

    signifHist.GetXaxis().SetRange(0,xmaxBin)
    signifHist.GetYaxis().SetRange(0,ymaxBin)
    signifHist.Draw("colz")
    aux.save("optimizeBinning2d_{}_{}".format(name,sigPoint), log=False)



if __name__ == "__main__":
    optimizeBinnigs2d(data, "data", "T5Wg_1550_1500")
    optimizeBinnigs2d(data, "data", "T5Wg_1550_100")
    optimizeBinnigs2d(gjets+qcd+ttjets+wjets+znunu+ttg+wg_mg+zg_130, "mc", "T5Wg_1550_1500")
    optimizeBinnigs2d(gjets+qcd+ttjets+wjets+znunu+ttg+wg_mg+zg_130, "mc", "T5Wg_1550_100")
