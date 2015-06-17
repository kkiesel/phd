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
from datasets import *

intLumi = 10000 # /pb

def randomName():
    """
    Generate a random string. This function is useful to give ROOT objects
    different names to avoid overwriting.
    """
    return "%x"%(randint(0, maxint))

def readHist( filename, histoname ):
    f = ROOT.TFile( filename )
    h = f.Get( histoname )
    if not h:
        print "Histogram {} not found in file {}".format(histoname, filename)
        return
    h = ROOT.gROOT.CloneObject( h )
    if not h.GetSumw2N():
        h.Sumw2()
    h.drawOption_ = ""
    return h

def getObjectNames( filename, path ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )

    outList = []
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()
        if isinstance( obj, ROOT.TH1 ):
            outList.append( element.GetName() )
        elif isinstance( obj, ROOT.TObjString ) or isinstance( obj, ROOT.TTree ):
            pass
        else:
            print "do not know what to do with", element.GetName()

    return outList

def getNgenDQM( filename ):
    h = readHist( filename, "DQMData/Run 1/Generator/Run summary/GenParticles/nEvt" )
    nGen = h.GetBinContent( 1 )
    return nGen

def getValAndError( val, err, sig=2 ):
    from math import floor, log10
    digit = sig - int(floor(log10(err))) - 1
    return "{} #pm {}".format( round(val,digit), round(err,digit) )

def getOwnStatBox( h, x1, y1 ):

    text = ROOT.TLatex()
    text.SetTextColor( h.GetLineColor() )

    mean = h.GetMean()
    mean_err = h.GetMeanError()

    text.DrawLatexNDC( x1, y1, "#mu = "+getValAndError( mean, mean_err ) )

def absHistWeighted( origHist ):
    origNbins = origHist.GetNbinsX()
    origXmin = origHist.GetBinLowEdge(1)
    origXmax = origHist.GetBinLowEdge(origHist.GetNbinsX()+1)
    if origXmin + origXmax:
        if origXmin:
            print "cant handle assymetric histograms"
        # else: print "already symetric?"
        return origHist

    h = ROOT.TH1F( "", origHist.GetTitle(), int(math.ceil(origNbins/2.)), 0, origXmax )

    for origBin in range( origNbins+2 ):
        newBin = int(abs(origBin - (origNbins+1.)/2)) + 1

        c1 = origHist.GetBinContent( origBin )
        e1 = origHist.GetBinError( origBin )
        c2 = h.GetBinContent( newBin )
        e2 = h.GetBinError( newBin )

        if e1 and e2:
            h.SetBinContent( newBin, ( c1*e1**-2 + c2*e2**-2 )/(e1**-2 + e2**-2) )
            h.SetBinError( newBin, 1./math.sqrt( e1**-2 + e2**-2 ) )

        else:
            h.SetBinContent( newBin, origHist.GetBinContent(origBin) )
            h.SetBinError( newBin, origHist.GetBinError(origBin) )

    return h

def getMean( h, eb=True ):
    (xmin, xmax) = (0, 1.5) if eb else (1.5, 5)

    if not h.Integral( h.FindBin(xmin), h.FindBin(xmax) ): return None
    h.Fit( "pol0", "0fcq", "", xmin, xmax )
    f = h.GetFunction( "pol0" )
    if not f: return 0
    return f.GetParameter(0)

def draw( files, path, name, config ):

    c = ROOT.TCanvas()

    processTex = ""
    match = re.match( ".*__RelVal(.*)_13__.*", files[0] )
    if match:
        regex = match.group(1)
        if regex == "H130GGgluonfusion":
            processTex = "H#rightarrow#gamma#gamma"
            processName = "Hgg"
        elif regex == "ZEE":
            processTex = "Z#rightarrowee"
            processName = "Zee"
        else: print "does not know what to do with", regex
    else:
        processTex = "Single e^{#minus}"
        processName = "closure"
    processTex += "  FullSim #color[2]{FastSim}"
    if len(files) > 2:
        processTex += " #color[4]{Mod}"

    hists = []
    for file in files:
        h = readHist( file, "{}/{}".format( path, name ) )
        if not h:
            return
        if not round( h.GetEntries()) or not round( h.Integral() ):
            print "no entries in {}".format( name )
            return

        if isinstance( h, ROOT.TH2 ):
            h = h.ProfileX( randomName() )
        if "VsEta" in name:
            h = absHistWeighted( h )

        hists.append( h )

    for h in hists:
        if not h.GetXaxis().GetTitle():
            h.SetXTitle( h.GetName() )

    hists[0].SetName("FullSim")
    if len(hists) > 1:
        hists[1].SetName("FastSim")
        hists[1].SetLineColor(2)
    if len(hists) > 2:
        hists[2].SetName("FastSim+Mod")
        hists[2].SetLineColor( ROOT.kBlue )

    for h in hists:
        # constumize histograms
        xmin = h.GetXaxis().GetXmin()
        xmax = h.GetXaxis().GetXmax()
        if config.has_option( name, "xmin" ): xmin = config.getfloat( name, "xmin")
        if config.has_option( name, "xmax" ): xmax = config.getfloat( name, "xmax" )
        h.GetXaxis().SetRangeUser( xmin, xmax )

        if config.has_option( name, "title" ):
            h.SetTitle( config.get( name, "title" )+ "  " )
        elif h.GetTitle():
            h.SetXTitle( h.GetTitle() +"    "+ h.GetXaxis().GetTitle() )
            h.SetTitle( "" )
        else:
            pass

        if config.has_option( name, "rebin" ): h.Rebin( config.getint( name, "rebin" ) )

        h.drawOption_ = "hist e"
        h.SetMarkerSize(0)


    for h in hists:
        if not isinstance( h, ROOT.TProfile ) and not "VsEta" in name:
            h.Scale( 1./h.GetEntries() )

    m = multiplot.Multiplot()
    for h in hists:
        m.add( h )
    m.Draw()



    label = ROOT.TLatex()
    label.DrawLatexNDC( .01, .96, "#font[61]{CMS} #scale[0.8]{#it{Simulation}}  "+processTex )
    label.DrawLatexNDC( .18, .88, "Private Work" )

    if len(hists) == 2:
        r = ratio.Ratio( "Full/Fast", hists[0], hists[1] )
    else:
        r = ratio.Ratio( "Full/Mod", hists[0], hists[2] )
    r.draw()
    if len(hists)>2:
        processName = "modIncl_"+processName

    if "VsEta" in name:
        ebFullMean = getMean( hists[0] )
        eeFullMean = getMean( hists[0], False )
        eb2Mean = getMean( hists[-1] )
        ee2Mean = getMean( hists[-1], False )

        agreementEB, agreementEE = 0, 0
        if ebFullMean and eb2Mean:
            agreementEB =  100*( abs(eb2Mean/ebFullMean) - 1 )
        if eeFullMean and ee2Mean:
            agreementEE =  100*( abs(ee2Mean/eeFullMean) - 1 )


        agreementLeg = ROOT.TLegend(.2, .3, .5, .5)
        agreementLeg.SetFillColor(0)
        agreementLeg.SetTextSize( hists[0].GetXaxis().GetLabelSize() )
        agreementLeg.SetTextFont( hists[0].GetXaxis().GetLabelFont() )
        agreementLeg.SetHeader("Agreement Mod,Full")
        agreementLeg.AddEntry( 0, "EB %.1f%%"%agreementEB, "" )
        agreementLeg.AddEntry( 0, "EE %.1f%%"%agreementEE, "" )
        agreementLeg.Draw()



    if name == "h_ele_PoPtrueVsEta":
        for bin in range( 1, r.ratio.GetNbinsX()+1 ):
            x = r.ratio.GetBinLowEdge( bin+1 )
            y = r.ratio.GetBinContent( bin )
            ey1 = r.ratio.GetBinError( bin )
            ey2 = r.totalUncert.GetBinError(bin)
            ey = ROOT.TMath.Sqrt( ey1**2 + ey2**2 )
            #if abs(y-1) < ey: y = 1.
            print "else if( genEta < %s ) { scale = %s; }"%(x,y)

    ROOT.gPad.GetCanvas().SaveAs("plots/%s_%s.pdf"%(processName, name ))

def getNprocessed( filename ):
    f = readHist( filename, "hCutFlow" )
    return int(f.GetBinContent(1))

def getHistoFromDataset( dataset, name ):
    h0 = None
    for i in range( len(dataset.files) ):
        h = readHist( dataset.files[i], name )
        h.Scale( dataset.xsecs[i]*dataset.ngens[i]/intLumi )
        h.SetLineColor( dataset.color )
        h.SetMarkerColor( dataset.color )

        if h0: h0.Add( h )
        else: h0 = h

    return h0

def save( name, folder="plots/", endings=[".pdf"] ):
    for ending in endings:
        ROOT.gPad.SaveAs( folder + name+ ending )

def compare( datasets, name, saveName ):
    m = multiplot.Multiplot()

    for d in datasets:
        h = getHistoFromDataset( d, name )
        if not h.Integral(): continue
        h.Scale( 1./h.Integral() )
        m.add( h, d.label )

    m.Draw()

    save( "compare%s_%s"%(saveName,name) )

def drawH2( dataset, name ):
    h = getHistoFromDataset( dataset, name )
    h.Draw("colz")
    save( "simpleH2_%s_%s"%(dataset.label,name) )


def compareAll( saveName="test", *datasets ):
    names = getObjectNames( datasets[0].files[0], "" )

    for name in names:
        if name.startswith("h_"):
            compare( datasets, name, saveName )
        #if name.startswith("h2_"):
        #    for d in datasets:
        #        drawH2( d, name )

def drawSameHistogram( saveName, name, data, bkg, additional ):

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    for d in bkg[-1::-1]:
        h = getHistoFromDataset( d, name )
        #if not h.Integral(): continue
        #h.Scale( 1./h.Integral() )
        m.addStack( h, d.label )

    m.Draw()

    save( "sameHistogram%s_%s"%(saveName,name) )
    can.SetLogy()
    save( "sameHistogram%s_%s_log"%(saveName,name) )


def drawSameHistograms( saveName="test", data=None, bkg=[], additional=[] ):
    names = getObjectNames( bkg[0].files[0], "" )

    names = ["h_met_loose"]

    for name in names:
        if name.startswith("h_"):
            drawSameHistogram( saveName, name, data, bkg, additional )


def main():
    #compareAll( "_all", gjets400, gjets600, znn400, znn600 )
    #compareAll( "_GjetsVsZnn", gjets, znn )
    #compareAll( "_allMC", gjets, znn, qcd, wjets )
    drawSameHistograms( "_allMC", bkg=[gjets, qcd, ttjets, wjets, znn ] )


if __name__ == "__main__":
    main()

