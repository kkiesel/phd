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

def save( name, folder="plots/", endings=[".pdf"] ):
    for ending in endings:
        ROOT.gPad.GetCanvas().SaveAs( folder + name+ ending )

def drawSame( fullName, fastName, hname ):
    fullh = aux.getFromFile( fullName, hname )
    fasth = aux.getFromFile( fastName, hname )
    if not fullh.Integral() or not fasth.Integral(): return

    for h in fullh, fasth:
        if "pt" in hname: h.GetXaxis().SetRangeUser(0, 150 )
        if "cIsoWorst" in hname: h.GetXaxis().SetRangeUser(0, 15 )

    fullh.SetLineColor( ROOT.kRed )
    fullh.drawOption_ = "hist"

    fasth.SetMarkerStyle(20)
    fasth.SetMarkerSize(0.3)
    fasth.drawOption_ = "e"

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add( fullh, "FullSim" )
    m.add( fasth, "FastSim" )

    m.Draw()
    r = ratio.Ratio( "Fast/Full", fasth, fullh )
    r.draw()

    save( "fastSimStudies_DY_"+hname )



def main():
    fullName = "../../DY_FullSim_hists.root"
    fastName = "../../DY_FastSim_hists.root"


    names = aux.getObjectNames( fullName )

    for name in names:
        drawSame( fullName, fastName, name )


if __name__ == "__main__":
    main()

