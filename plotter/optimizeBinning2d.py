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

def asimovSig(s, b):
    if not b: b = 1e-10
    return math.sqrt(2*((s+b)*math.log(1+s/b)-s))

def optimizeBinnigs2d(bkgSet, name, signalSet):
    style.style2d()
    bkgHist = bkgSet.getHist("tr/met_vs_emht")
    sigHist = signalSet.getHist("tr/met_vs_emht")

    lumiScale=40.00/27.22
    bkgHist.Scale(lumiScale)
    sigHist.Scale(lumiScale)

    xmax = 800
    ymax = 4000
    xmaxBin = bkgHist.GetXaxis().FindBin(xmax)
    ymaxBin = bkgHist.GetYaxis().FindBin(ymax)

    signifHist = sigHist.Clone(aux.randomName())
    for xbin, ybin in aux.loopH(sigHist):
        if xbin > xmaxBin or ybin > ymaxBin: continue
        b = bkgHist.Integral(xbin, -1, ybin, -1)
        s = sigHist.Integral(xbin, -1, ybin, -1)
        signif = asimovSig(s, b)
        signifHist.SetBinContent(xbin, ybin, signif)

    signifHist.GetXaxis().SetRange(0,xmaxBin)
    signifHist.GetYaxis().SetRange(0,ymaxBin)
    signifHist.SetXTitle("minimum "+signifHist.GetXaxis().GetTitle())
    signifHist.SetYTitle("minimum "+signifHist.GetYaxis().GetTitle())
    signifHist.SetZTitle("Asimov significance")
    signifHist.Draw("colz")
    l = aux.Label(info=signalSet.label, drawAll=False)
    l.lum = ROOT.TLatex( .63, .95, "%.2f fb^{-1} (%s TeV)"%(aux.intLumi/1000.*lumiScale, 13) )
    l.draw()

    aux.save("optimizeBinning2d_{}_{}".format(name,signalSet.names[0]), log=False)

def optimizeMetBinning(bkgSet, signalSet, name, minHt, maxHt):
    style.defaultStyle()
    bkgHist = bkgSet.getHist("tr/met_vs_emht")
    sigHist = signalSet.getHist("tr/met_vs_emht")

    lumiScale=40.00/27.22
    bkgHist.Scale(lumiScale)
    sigHist.Scale(lumiScale)

    minHtBin = bkgHist.GetYaxis().FindFixBin(minHt)
    maxHtBin = bkgHist.GetYaxis().FindFixBin(maxHt-1e6)
    bkg1d = bkgHist.ProjectionX(aux.randomName(), minHtBin, maxHt)
    sig1d = sigHist.ProjectionX(aux.randomName(), minHtBin, maxHt)

    signifHist = sig1d.Clone(aux.randomName())
    binning = [sig1d.GetNbinsX()+2]
    while binning[-1] > 1:
        maxSig = 0
        maxSigBin = 0
        for bin in range(binning[-1], -1, -1):
            b = bkg1d.Integral(bin, binning[-1])
            s = sig1d.Integral(bin, binning[-1])
            signif = asimovSig(s, b)
            if signif > maxSig:
                maxSig = signif
                maxSigBin = bin
        binning.append(maxSigBin)
    binning = binning[-1:1:-1]
    print [sig1d.GetBinLowEdge(b) for b in binning]



if __name__ == "__main__":
    #optimizeBinnigs2d(data, "data", t5wg_1600_100)
    optimizeMetBinning(data, t5wg_1600_100, "data_highHt", 2000, 5000)
    optimizeMetBinning(data, t5wg_1600_100, "data_mediumHt", 1200, 2000)
    optimizeMetBinning(data, t5wg_1600_100, "data_lowHt", 0, 1200)
    #optimizeBinnigs2d(gjets+qcd+ttjets+wjets+znunu+ttg+wg_mg+zg_130, "mc", "T5Wg_1550_100")
