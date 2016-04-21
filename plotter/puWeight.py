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
    h_mc2 = aux.getFromFile( path+"nTrueVertexMC.root", "mix_2015_25ns_FallMC_matchData_PoissonOOTPU" )
    h_weight = aux.getFromFile( path+"data/puWeights.root", "pileupWeight_mix_2015_25ns_Startup_PoissonOOTPU" )
    h_weight2 = aux.getFromFile( path+"data/puWeights.root", "pileupWeight_mix_2015_25ns_FallMC_matchData_PoissonOOTPU" )

    for h in h_data, h_mc, h_mc2, h_weight, h_weight2:
        h.SetTitle(";number of pileup events;Events")
        h.GetXaxis().SetRangeUser(0, 30 )

    h_data.drawOption_ = "ep"
    h_data.SetLineColor( ROOT.kBlack )
    h_data.SetMarkerColor( ROOT.kBlack )
    h_data.SetMarkerStyle( 20 )
    h_data.SetMarkerSize( 0.5 )

    h_mc.Scale( h_data.Integral() )
    h_mc.SetLineColor( ROOT.kBlue )

    h_mc2.Scale( h_data.Integral() )
    h_mc2.SetLineColor( ROOT.kRed )

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.leg.SetY1(.7)
    m.leg.SetX1(.6)
    m.leg.SetX2(.99)
    m.add( h_data, "Data" )
    m.add( h_mc2, "76X (SM)" )
    m.add( h_mc, "74X (Signal)" )

    m.Draw()

    # just create the pad, the ratio is not used
    r = ratio.Ratio( "Data/Sim.", h_data, h_mc )
    r.draw()

    h_weight.SetMaximum(3)
    h_weight.SetTitle( ";number of pileup events;Data/Sim." )
    h_weight.SetLineColor(h_mc.GetLineColor())
    h_weight.Draw("hist")
    h_weight2.SetLineColor(h_mc2.GetLineColor())
    h_weight2.Draw("hist same")


    l = aux.Label()

    can.SaveAs("plots/pu_weight.pdf")

