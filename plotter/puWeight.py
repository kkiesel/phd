#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import ROOT
import os


# private libs
import ratio
import style
import multiplot

import auxiliary as aux

if __name__ == "__main__":
    path = os.environ['CMSSW_BASE']+'/src/TreeWriter/PUreweighting/'

    h_data = aux.getFromFile( path+"nTrueVertexData.root", "pileup" )
    h_mc = aux.getFromFile( path+"nTrueVertexMC.root", "mix_2015_25ns_Startup_PoissonOOTPU" )
    h_weight = aux.getFromFile( path+"data/puWeights.root", "pileupWeight_mix_2015_25ns_Startup_PoissonOOTPU" )

    for h in h_data, h_mc,h_weight:
        h.SetTitle(";number of pileup events;Events")
        h.GetXaxis().SetRangeUser(0, 30 )


    h_data.drawOption_ = "ep"
    h_data.SetLineColor( ROOT.kBlack )
    h_data.SetMarkerColor( ROOT.kBlack )
    h_data.SetMarkerStyle( 20 )
    h_data.SetMarkerSize( 0.5 )

    h_mc.Scale( h_data.Integral() )
    h_mc.SetLineColor( ROOT.kBlue )


    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.leg.SetY1(.8)
    m.leg.SetX1(.7)
    m.leg.SetX2(.99)
    m.add( h_data, "Data" )
    m.add( h_mc, "Simulation" )

    m.Draw()

    # just create the pad, the ratio is not used
    r = ratio.Ratio( "Data/Sim.", h_data, h_mc )
    r.draw()

    h_weight.SetMaximum(3)
    h_weight.SetTitle( ";number of pileup events;Data/Sim." )
    h_weight.Draw("hist")


    l = aux.Label()

    can.SaveAs("plots/pu_weight.pdf")

