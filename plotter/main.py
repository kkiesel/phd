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

def getNprocessed( filename ):
    f = aux.getFromFile( filename, "hCutFlow" )
    return int(f.GetBinContent(1))

def getHistoFromDataset( dataset, name ):
    h0 = None
    for i in range( len(dataset.files) ):
        h = aux.getFromFile( dataset.files[i], name )
        if isinstance( h, ROOT.TH1 ):
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
    names = aux.getObjectNames( datasets[0].files[0], "" )

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
    names = aux.getObjectNames( bkg[0].files[0], "" )

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


def ewkClosure( dataset, samplename="" ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TH1F] )

    gSet = "loose_genElectron"
    eSet = "looseElectron"

    for name in names:
        if gSet not in name: continue
        m = multiplot.Multiplot()

        h = getHistoFromDataset( dataset, name )
        h.SetLineColor(1)
        h.SetMarkerColor(1)
        m.add( h, "#gamma (gen e)" )

        h = getHistoFromDataset( dataset, name.replace( gSet, eSet ) )
        h.Scale( 0.01 )
        h.drawOption_ = "hist"
        m.add( h, "0.01 #times #gamma_{pixel}" )

        m.Draw()

        l = aux.Label()
        save( "ewkClosure_"+name+samplename )

def drawROCs():
    a = aux.getROC( getHistoFromDataset( gjets, "h_g_mva_base" ), getHistoFromDataset( qcd, "h_g_mva_base" ) )
    a.Draw()
    a.GetXaxis().SetRangeUser(0.1,1)
    a.GetYaxis().SetRangeUser(0.1,1)
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetLogx()
    save("test")

def efficiencies( dataset ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TEfficiency] )

    for name in names:
        h = getHistoFromDataset( dataset, name )
        if h.UsesWeights(): h.SetStatisticOption( ROOT.TEfficiency.kFNormal )
        h.Draw()

        l = aux.Label()
        save( "efficiency_"+name )





def main():
    #compareAll( "_all", gjets400, gjets600, znn400, znn600 )
    #compareAll( "_GjetsVsZnn", gjets, znn )
    #compareAll( "_allMC", gjets, znn, qcd, wjets )
    #drawSameHistograms( "_allMC", bkg=[gjets, qcd, ttjets, wjets] )
    #drawRazor( ttjets )

    #ewkClosure( ttjets, "_tt" )
    #ewkClosure( wjets, "_w" )
    #ewkClosure( wjets+ttjets, "_ewk" )

    efficiencies( ttjets+qcd+gjets+wjets )



if __name__ == "__main__":
    main()

