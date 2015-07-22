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
from datasets import *
intLumi = 5000 # /pb

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

def getObjectNames( filename, path="", objects=[ROOT.TH1] ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )

    outList = []
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()

        if any( [ isinstance( obj, o ) for o in objects ] ):
            outList.append( element.GetName() )

    return outList

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

def getNprocessed( filename ):
    f = readHist( filename, "hCutFlow" )
    return int(f.GetBinContent(1))

def getHistoFromDataset( dataset, name ):
    h0 = None
    for i in range( len(dataset.files) ):
        h = readHist( dataset.files[i], name )
        h.Scale( intLumi * dataset.xsecs[i] / dataset.ngens[i] )
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
        h.SetYTitle( aux.getYAxisTitle( h ) )
        #if not h.Integral(): continue
        #h.Scale( 1./h.Integral() )
        m.addStack( h, d.label )

    m.Draw()

    l = aux.Label()
    save( "sameHistogram%s_%s"%(saveName,name) )
    can.SetLogy()
    save( "sameHistogram%s_%s_log"%(saveName,name) )


def drawSameHistograms( saveName="test", data=None, bkg=[], additional=[] ):
    names = getObjectNames( bkg[0].files[0], "" )

#    names = ["h_met_loose"] # test plot

    for name in names:
        if name.startswith("h_"):
            drawSameHistogram( saveName, name, data, bkg, additional )

def getProjections( h2, alongX=True ):
    hs = []
    label = h2.GetYaxis().GetTitle()

    for ybin in range( h2.GetNbinsY()+2 ):

        ylow = h2.GetYaxis().GetBinLowEdge(ybin)
        yhigh = h2.GetYaxis().GetBinUpEdge(ybin)
        name = "{} #leq {} < {}".format( ylow, label, yhigh )
        if ybin == 0: name = "{} < {}".format( label, yhigh )
        if ybin == h2.GetNbinsY()+1: name = "{} #leq {}".format( ylow, label )


        h = h2.ProjectionX( name, ybin, ybin )
        h.SetLineColor( ybin+2)
        if h.GetEntries():
            h.Scale( 1./h.GetEntries() )
            hs.append( h )

    return hs


def drawRazor( dataset ):
    h2 = getHistoFromDataset( dataset, "h2_razorPlane" )
    h2.Rebin2D( 1, 20 )
    razorFit = ROOT.TF2("razorFitFunc", "[0]*( [1]*(x[0]-[2])*(x[1]-[3]) - 1 ) * exp( -[1]*(x[0]-[2])*(x[1]-[3]) )", 0, 2000, 0, 0.5 )
    razorFit.SetParameters( h2.GetEntries(), 0.0005, 170, 0.00001 )
    razorFit.FixParameter( 2, 170 )
    fr = h2.Fit( "razorFitFunc" )
    h2.Draw("cont2")
    razorFit.Draw("same")
    save( "razorPlane" )

    pX = getProjections( h2 )
    pX[0].Draw()
    for h in pX[1:]: h.Draw("same")

    leg = ROOT.TLegend( .7, .7, .95, .95 )
    for h in pX: leg.AddEntry( h, h.GetName(), "l" )
    leg.Draw()

    save( "razorAlongX" )


def ewkClosure( dataset ):
    names = getObjectNames( dataset.files[0] )

    gSet = "loose_genElectron"
    eSet = "looseElectron"

    for name in names:
        if gSet not in name: continue
        m = multiplot.Multiplot()

        h = getHistoFromDataset( dataset, name )
        h.SetLineColor(1)
        h.SetMarkerColor(1)
        h.SetMarkerSize(20)
        m.add( h, "#gamma (gen e)" )

        h = getHistoFromDataset( dataset, name.replace( gSet, eSet ) )
        h.Scale( 0.01 )
        m.add( h, "0.01 #times #gamma_{pixel}" )

        m.Draw()

        save( "ewkClosure_"+name )




def main():
    #compareAll( "_all", gjets400, gjets600, znn400, znn600 )
    #compareAll( "_GjetsVsZnn", gjets, znn )
    #compareAll( "_allMC", gjets, znn, qcd, wjets )
    #drawSameHistograms( "_allMC", bkg=[gjets, qcd, ttjets, wjets] )
    #drawRazor( ttjets )

    ewkClosure( ttjets )


if __name__ == "__main__":
    main()

