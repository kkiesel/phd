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
import rebinner
intLumi = 2000 # /pb
intLumi = 166.37 # /pb

def getNprocessed( filename ):
    f = aux.getFromFile( filename, "hCutFlow" )
    return int(f.GetBinContent(1))

def getHistoFromDataset( dataset, name ):
    h0 = None
    for i in range( len(dataset.files) ):
        h = aux.getFromFile( dataset.files[i], name )
        if isinstance( h, ROOT.TH1 ):
            if dataset.xsecs[i]:
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

def drawSameHistogram( saveName, name, data, bkg, additional, binning=None ):

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    for d in bkg[-1::-1]:
        h = getHistoFromDataset( d, name )
        if binning:
            h = aux.rebin( h, binning )
        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        if not h.Integral(): continue
        #h.Scale( 1./h.Integral() )
        m.addStack( h, d.label )

    for d in additional:
        h = getHistoFromDataset( d, name )
        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption = "ep"
            h.SetMarkerStyle(20)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
        else:
            h.drawOption_ = "hist"
        m.add( h, d.label )

    if m.Draw():

        l = aux.Label()
        save( "sameHistogram%s_%s"%(saveName,name) )
        can.SetLogy()
        save( "sameHistogram%s_%s_log"%(saveName,name) )


def drawSameHistograms( saveName="test", data=None, bkg=[], additional=[] ):
    names = aux.getObjectNames( bkg[0].files[0], "" )

    #names = ["h_met__tr"] # test plot
    for name in names:
        if name.startswith("h_"):
            drawSameHistogram( saveName, name, data, bkg, additional )
            hname = re.match( "(.*)__tr", name ).group(1)
            for binningName, binnings in rebinner.cfg.items( hname ):
                binnings = [ float(x) for x in binnings.split(" ") ]
                drawSameHistogram( saveName+"_bin"+binningName, name, data, bkg, additional, binnings )

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

def qcdClosure( dataset, samplename="" ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TH1F] )

    gSet = "loose"
    eSet = "looseJet"

    for name in names:
        can = ROOT.TCanvas()
        if eSet not in name: continue
        m = multiplot.Multiplot()

        h1 = getHistoFromDataset( dataset, name.replace( eSet, gSet ) )
        h1.SetLineColor(1)
        h1.SetMarkerColor(1)
        m.add( h1, "#gamma" )

        h = getHistoFromDataset( dataset, name )
        if h.Integral(): h.Scale( h1.Integral()/h.Integral() )
        h.drawOption_ = "hist"
        m.add( h, "#gamma-like" )

        m.Draw()

        l = aux.Label()
        save( "qcdClosure_"+name+samplename )
        can.SetLogy()
        save( "qcdClosure_"+name+samplename+"_log" )

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

    for hname in ["h_gCol_mva__base", "h_g_cIso__base", "h_g_nIso__base", "h_g_pIso__base", "h_g_sigmaIetaIeta__base", "h_g_hOverE__base"]:
        c = ROOT.TCanvas()

        effSig = getHistoFromDataset( gjets, hname )
        effBkg = getHistoFromDataset( qcd, hname )

        a = aux.getROC( effSig.GetPassedHistogram(), effBkg.GetTotalHistogram(), "mva" in hname )
        if not a: continue
        a.Draw()
        a.GetXaxis().SetRangeUser(0.01,1)
        a.GetYaxis().SetRangeUser(0.01,1)
        l = aux.Label()
        save("roc_"+hname)

def efficiencies( dataset, savename="" ):
    names = aux.getObjectNames( dataset.files[0], "", [ROOT.TEfficiency] )

    for name in names:
        h = getHistoFromDataset( dataset, name )
        if h.UsesWeights(): h.SetStatisticOption( ROOT.TEfficiency.kFNormal )

        h_pas = h.GetPassedHistogram()
        h_tot = h.GetTotalHistogram()
        if "hlt" in name:
            eff = ROOT.TEfficiency( h_pas, h_tot )
        else:
            eff = h
        eff.Draw()

        if name == "eff_hlt_ht__base":
            cutValue = 600
        elif name == "eff_hlt_pt__base":
            cutValue = 100
        else:
            cutValue = None

        if cutValue:
            bin = h_pas.FindFixBin( cutValue )
            passed = int(h_pas.Integral( bin, -1 ))
            total = int(h_tot.Integral( bin, -1 ))
            conf = 0.682689492137
            e = 1.*passed/total
            e_up = ROOT.TEfficiency.ClopperPearson( total, passed, conf, True )
            e_dn = ROOT.TEfficiency.ClopperPearson( total, passed, conf, False )
            print "ε = {:.1%} + {:.1%} - {:.1%}".format(e, e_up-e,e-e_dn )

            # graphical representation
            l = ROOT.TLine()
            l.SetLineWidth(2)
            l.SetLineColor( ROOT.kRed )
            xmax = eff.CreateGraph().GetHistogram().GetXaxis().GetXmax()
            l.DrawLine( cutValue, e, xmax, e )

            l.SetLineStyle(2)
            ymin = eff.CreateGraph().GetHistogram().GetYaxis().GetXmin()
            ymax = eff.CreateGraph().GetHistogram().GetYaxis().GetXmax()
            l.DrawLine( cutValue, ymin, cutValue, ymax )


        l = aux.Label()
        save( "efficiency_"+savename+name )





def main():
    #compareAll( "_all", gjets400, gjets600, znn400, znn600 )
    #compareAll( "_GjetsVsZnn", gjets, znn )
    #compareAll( "_allMC", gjets, znn, qcd, wjets )
    drawSameHistograms( "_mc_data", bkg=[gjets, qcd, ttjets, wjets, ttg], additional=[data])
    #drawSameHistograms( "_allMC", bkg=[gjets, qcd, ttjets, wjets, ttg], additional=[t2ttgg, t5gg, t5hg])
    #drawSameHistograms( "_QCD", bkg=[ qcd2000, qcd1500, qcd1000, qcd700, qcd500, qcd300 ] )
    #drawRazor( ttjets )

    #ewkClosure( ttjets, "_tt" )
    #ewkClosure( wjets, "_w" )
    #ewkClosure( wjets+ttjets, "_ewk" )
    #qcdClosure( qcd+gjets, "_qcd-gjets" )

    #efficiencies( ttjets+qcd+gjets+wjets, "allMC_" )
    #efficiencies( data, "singlePhoton_" )
    #efficiencies( jetHt, "jetHt_" )

    #drawROCs()



if __name__ == "__main__":
    main()

