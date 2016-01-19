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

def infoText( name ):
    infoReplacements = {
#        "_noSel":"no photon id",
#        "_loose":"",
#        "_loose_eb":"EB",
#        "_loose_ee":"EE",
#        "_loose_eb_genPhoton":"EB,gen #gamma",
#        "_loose_ee_genPhoton":"EE,gen #gamma",
#        "_loose_eb_genElectron":"EB,gen e",
#        "_loose_ee_genElectron":"EE,gen e",

        "_loose_genPhoton":"gen #gamma",
        "_loose_eb_genPhoton":"EB,gen #gamma",
        "_loose_ee_genPhoton":"EE,gen #gamma",
        "_medium_genPhoton":"mediumID,gen #gamma",
        "_medium_eb_genPhoton":"mediumID,EB,gen #gamma",
        "_medium_ee_genPhoton":"mediumID,EE,gen #gamma",
        "_tight_genPhoton":"tightID,gen #gamma",
        "_tight_eb_genPhoton":"tightID,EB,gen #gamma",
        "_tight_ee_genPhoton":"tightID,EE,gen #gamma",
        "_fake":"fake",
        "_fake_eb":"EB,fake",
        "_fake_ee":"EE,fake",
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


def drawSame( fullName, fastName, hname, binning=[] ):
#    if hname not in ["pt_loose_ee_genPhoton","pt_loose_eb_genPhoton"]: return
    ds = datasetAbbr( fullName )

    fullNgen = aux.getNgen(fullName)
    fastNgen = aux.getNgen(fastName)

    fullh = aux.getFromFile( fullName, hname )
    fasth = aux.getFromFile( fastName, hname )
    if not fullh.Integral() or not fasth.Integral(): return


    if hname in ["hasPixelSeed_loose_eb_genElectron"]:
        print "$f_\\text{{\\fullsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fullh.GetBinContent(1)/fullh.Integral()*100. )
        print "$f_\\text{{\\fastsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fasth.GetBinContent(1)/fasth.Integral()*100. )
    if hname in ["passEVeto_loose_eb_genElectron"]:
        print "$f_\\text{{\\fullsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fullh.GetBinContent(2)/fullh.Integral()*100. )
        print "$f_\\text{{\\fastsim}}(e\\rightarrow\gamma)={:.2f}\\%$".format( fasth.GetBinContent(2)/fasth.Integral()*100. )

    fasth.Scale( 1./fastNgen )
    fullh.Scale( 1./fullNgen )

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
    r = ratio.Ratio( "#font[132]{F#scale[0.75]{ULL}/#color[2]{F#scale[0.75]{AST}}}", fullh, fasth )
    r.ratioStat.SetLineColor(fasth.GetLineColor())
    r.ratio.SetMarkerStyle(7)
    r.draw()

    scale = fullh.Integral(0,-1)/fasth.Integral(0,-1)
    fastInt, fastErr = aux.integralAndError(fasth)
    scaleError = scale * fastErr/fastInt

    info = ROOT.TLatex( 0.25, .95, "{}  {}  #scale[0.75]{{#Deltas=({:.2f}#pm{:.2f})%}}".format(infoText(hname),ds, 100*(scale-1),100*(scaleError)) )
    info.SetNDC()
    info.Draw()

    aux.save( "fastSimStudies_"+ds+"_"+hname )
    can.SetLogy(1)
    aux.save( "fastSimStudies_"+ds+"_"+hname+"_log" )

def scaleFactors( fullName, fastName ):
    ds = datasetAbbr( fullName )

    fullNGen = aux.getNgen( fullName )
    fastNGen = aux.getNgen( fastName )

    fullh = aux.getFromFile( fullName, "isLoose_loose_eb_genPhoton" )
    fasth = aux.getFromFile( fastName, "isLoose_loose_eb_genPhoton" )

    print "full/fast", fullh.GetEntries()/fasth.GetEntries() * fastNGen/fullNGen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample', default=["DY"], nargs="+" )
    args = parser.parse_args()

    for ds in args.sample:

        fullName = "../fastSimComparison/%s_FullSim_hists.root"%ds
        fastName = "../fastSimComparison/%s_FastSim_hists.root"%ds


        names = aux.getObjectNames( fullName )
        for name in names: drawSame( fullName, fastName, name )
        #drawSame( fullName, fastName, "cIso_loose", aux.drange(0, 3, 10 ) )
        #drawSame( fullName, fastName, "nIso_loose", aux.drange(0, 8, 10 ) )
        #drawSame( fullName, fastName, "pIso_loose", aux.drange(0, 5, 10 ) )

        scaleFactors( fullName, fastName )

if __name__ == "__main__":
    main()

