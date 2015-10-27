#!/usr/bin/env python2
# -*- coding: utf-8 -*-

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

import auxiliary as aux

def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def drange(start, stop, n):
    out = [start]
    step = 1.*(stop-start)/n
    while out[-1] < stop:
        out.append( out[-1] + step )
    return out



def save( name, folder="plots/", endings=[".pdf"] ):
    for ending in endings:
        ROOT.gPad.GetCanvas().SaveAs( folder + name+ ending )

def drawSame( fullName, fastName, hname, binning=[] ):
    fullh = aux.getFromFile( fullName, hname )
    fasth = aux.getFromFile( fastName, hname )
    if not fullh.Integral() or not fasth.Integral(): return
    if binning:
        fullh = aux.rebin( fullh, binning )
        fasth = aux.rebin( fasth, binning )

    for h in fullh, fasth:
        if "pt" in hname: h.GetXaxis().SetRangeUser(0, 150 )
        if "cIsoWorst" in hname: h.GetXaxis().SetRangeUser(0, 15 )

    fullh.drawOption_ = "hist"

    fasth.SetMarkerColor( ROOT.kRed )
    fasth.SetLineColor( ROOT.kRed )
    fasth.SetMarkerStyle(20)
    fasth.SetMarkerSize(0.4)
    fasth.drawOption_ = "pe"

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add( fullh, "FullSim" )
    m.add( fasth, "FastSim" )

    m.Draw()
    r = ratio.Ratio( "Fast/Full", fasth, fullh )
    r.draw()

    save( "fastSimStudies_DY_"+hname )
    can.SetLogy(1)
    save( "fastSimStudies_DY_"+hname+"_log" )




def main():
    fullName = "../fastSimComparison/DY_FullSim_hists.root"
    fastName = "../fastSimComparison/DY_FastSim_hists.root"


    names = aux.getObjectNames( fullName )

    for name in names: drawSame( fullName, fastName, name )

    drawSame( fullName, fastName, "eta_loose", drange(0, 2.5, 100) )
    drawSame( fullName, fastName, "pt_loose_eb" )
    drawSame( fullName, fastName, "pt_loose_ee" )
    drawSame( fullName, fastName, "cIso_loose", drange(0, 10, 100 ) )
    drawSame( fullName, fastName, "nIso_loose", drange(0, 10, 100 ) )
    drawSame( fullName, fastName, "pIso_loose", drange(0, 10, 100 ) )


if __name__ == "__main__":
    main()

