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

def infoText( name ):
    infoReplacements = {
        "_noSel":"no photon id",
        "_loose":"",
        "_loose_eb":"EB",
        "_loose_ee":"EE",
        "_loose_eb_genPhoton":"EB,gen #gamma",
        "_loose_ee_genPhoton":"EE,gen #gamma",
        "_loose_eb_genElectron":"EB,gen e",
        "_loose_ee_genElectron":"EE,gen e"
    }
    for n, repl in infoReplacements.iteritems():
        if name.endswith(n):
            return repl
    print "no replacement for", name
    return ""

def datasetAbbr( filename ):
    import re
    m = re.match( ".*/(.*)_F.*.root", filename )
    if m:
        return m.group(1)
    else:
        print "No dataset name found in file", filename
        return "Dataset"


def save( name, folder="plots/", endings=[".pdf"] ):
    for ending in endings:
        ROOT.gPad.GetCanvas().SaveAs( folder + name+ ending )

def drawSame( fullName, fastName, hname, binning=[] ):
    ds = datasetAbbr( fullName )

    fullh = aux.getFromFile( fullName, hname )
    fasth = aux.getFromFile( fastName, hname )
    if not fullh.Integral() or not fasth.Integral(): return


    if hname in ["hasPixelSeed_loose_eb_genElectron"]:
        print "$f_\\text{{\\fullsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fullh.GetBinContent(1)/fullh.Integral()*100. )
        print "$f_\\text{{\\fastsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fasth.GetBinContent(1)/fasth.Integral()*100. )
    if hname in ["passEVeto_loose_eb_genElectron"]:
        print "$f_\\text{{\\fullsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fullh.GetBinContent(2)/fullh.Integral()*100. )
        print "$f_\\text{{\\fastsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fasth.GetBinContent(2)/fasth.Integral()*100. )

    fasth.Scale( fullh.GetEntries()/fasth.GetEntries() )
    if binning:
        fullh = aux.rebin( fullh, binning )
        fasth = aux.rebin( fasth, binning )

    if ds == "QCD":
        if hname in ["pt_loose", "eta_loose", "pt_loose_eb", "pt_loose_ee"]:
            for h in fullh, fasth:
                h.Rebin(4)

    if "hoe" in hname or "Iso" in hname:
        for h in fullh, fasth:
            h.GetXaxis().SetRangeUser(-0.1, h.GetXaxis().GetXmax())


    fullh.drawOption_ = "hist"

    fasth.SetMarkerColor( ROOT.kRed )
    fasth.SetLineColor( ROOT.kRed )
    fasth.SetMarkerStyle(20)
    fasth.SetMarkerSize(0.6)
    fasth.drawOption_ = "pe"

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add( fullh, "#font[132]{F#scale[0.75]{ULL}S#scale[0.75]{IM}}" )
    m.add( fasth, "#color[2]{#font[132]{F#scale[0.75]{AST}S#scale[0.75]{IM}}}" )

    m.Draw()
    r = ratio.Ratio( "#font[132]{#color[2]{F#scale[0.75]{AST}}/F#scale[0.75]{ULL}}", fasth, fullh )
    r.draw()


    info = ROOT.TLatex( 0.4, .95, infoText(hname)+"  "+ds )
    info.SetNDC()
    info.Draw()

    save( "fastSimStudies_"+ds+"_"+hname )
    can.SetLogy(1)
    save( "fastSimStudies_"+ds+"_"+hname+"_log" )




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample', default=["DY"], nargs="+" )
    args = parser.parse_args()

    for ds in args.sample:

        fullName = "../fastSimComparison/%s_FullSim_hists.root"%ds
        fastName = "../fastSimComparison/%s_FastSim_hists.root"%ds


        names = aux.getObjectNames( fullName )
        for name in names: drawSame( fullName, fastName, name )
        drawSame( fullName, fastName, "eta_loose", drange(0, 2.5, 100) )
        drawSame( fullName, fastName, "cIso_loose", drange(0, 3, 10 ) )
        drawSame( fullName, fastName, "nIso_loose", drange(0, 8, 10 ) )
        drawSame( fullName, fastName, "pIso_loose", drange(0, 5, 10 ) )


if __name__ == "__main__":
    main()

