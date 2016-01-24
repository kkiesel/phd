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
    if isinstance( fullh, ROOT.TH3 ): return
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

def getYbinCenterOfShiftedBin( h3, val, shift ):
    return h3.GetYaxis().GetBinCenter(h3.GetYaxis().FindFixBin( val ) + shift )

def h3analyser( fullName, fastName ):
    genName = "h3_pt_eta_nV_gen"
    full3gen = aux.getFromFile( fullName, genName )
    fast3gen = aux.getFromFile( fastName, genName )

    ebEtaMax = getYbinCenterOfShiftedBin( full3gen, 1.4442, -1 )
    eeEtaMin = getYbinCenterOfShiftedBin( full3gen, 1.556, +1 )
    eeEtaMax = getYbinCenterOfShiftedBin( full3gen, 2.5, -1 )

    cuts = [
        ( "", 0, 5000, 0, 5000, 0, 5000 ),
        ( "EB", 0, 5000, 0, ebEtaMax, 0, 5000 ),
        ( "EE", 0, 5000, eeEtaMin, eeEtaMax, 0, 5000 ),
        ( "pt40", 40, 5000, 0, 5000, 0, 5000 ),
        ( "pt500", 500, 5000, 0, 5000, 0, 5000 ),
        ]



    for compareName in ["h3_pt_eta_nV_loose","h3_pt_eta_nV_medium","h3_pt_eta_nV_tight","h3_pt_eta_nV_fake"]:
        idName = compareName.split("_")[-1]
        full3com = aux.getFromFile( fullName, compareName )
        fast3com = aux.getFromFile( fastName, compareName )

        for axis in "XYZ":

            for cutn, xmin,xmax,ymin,ymax,zmin,zmax in cuts:
                can = ROOT.TCanvas()

                for h in full3gen,full3com,fast3gen,fast3com:
                    h.GetXaxis().SetRangeUser(xmin,xmax)
                    h.GetYaxis().SetRangeUser(ymin,ymax)
                    h.GetZaxis().SetRangeUser(zmin,zmax)

                full1gen = full3gen.Project3D(axis+"gen")
                full1com = full3com.Project3D(axis+"full")
                fast1gen = fast3gen.Project3D(axis+"gen")
                fast1com = fast3com.Project3D(axis+"fast")
                if axis == "X":
                    binning = range(0,1000,20)+range(1000,2000,50)+range(2000,3000,100)
                    full1gen = aux.rebin(full1gen, binning, False)
                    full1com = aux.rebin(full1com, binning, False)
                    fast1gen = aux.rebin(fast1gen, binning, False)
                    fast1com = aux.rebin(fast1com, binning, False)

                fullEff = ROOT.TEfficiency( full1com, full1gen )
                fastEff = ROOT.TEfficiency( fast1com, fast1gen )
                fastEff.SetLineColor(ROOT.kRed)

                fullEff.Draw()
                fastEff.Draw("same")

                can.Update()

                fullgr = fullEff.GetPaintedGraph()
                fastgr = fastEff.GetPaintedGraph()
                n = fullgr.GetN()
                x = fullgr.GetX()
                fully = fullgr.GetY()
                fullyEU = fullgr.GetEYlow()
                fullyED = fullgr.GetEYhigh()

                fasty = fastgr.GetY()
                fastyEU = fastgr.GetEYlow()
                fastyED = fastgr.GetEYhigh()

                rGr = fullgr.Clone( aux.randomName() )

                for i in range(n):
                    if not fasty[i] or not fully[i]:
                        rGr.RemovePoint(i)
                        continue
                    r = fully[i]/fasty[i]
                    eu = r * math.sqrt( abs( (fullyEU[i]/fully[i])**2 - (fastyEU[i]/fasty[i])**2 ) )
                    ed = r * math.sqrt( abs( (fullyED[i]/fully[i])**2 - (fastyED[i]/fasty[i])**2 ) )
                    rGr.SetPoint(i,x[i],r)
                    rGr.SetPointError(i,0,0,eu, ed)

                ratio.clearXaxisCurrentPad()
                p = ratio.createBottomPad(0.5)

                rGr.SetMaximum(1.05)
                rGr.SetMinimum(0.95)
                rGr.Draw("ap")
                rGr.GetYaxis().SetTitle("Full/Fast")

                aux.save("fastSimStudies_{}_{}_{}".format(idName,axis,cutn) )
                return



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample', default=["DY"], nargs="+" )
    args = parser.parse_args()

    for ds in args.sample:

        fullName = "../fastSimComparison/%s_FullSim_hists.root"%ds
        fastName = "../fastSimComparison/%s_FastSim_hists.root"%ds


        #names = aux.getObjectNames( fullName )
        #for name in names: drawSame( fullName, fastName, name )
        h3analyser( fullName, fastName )

        #scaleFactors( fullName, fastName )

if __name__ == "__main__":
    main()

