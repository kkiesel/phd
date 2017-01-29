#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from include import *

if __name__ == "__main__":
    path = os.environ['CMSSW_BASE']+'/src/TreeWriter/PUreweighting/'

    puName = "mix_2016_25ns_Moriond17MC_PoissonOOTPU"

    h_data = aux.getFromFile(path+"nTrueVertexData.root", "pileup")
    h_mc = aux.getFromFile(path+"nTrueVertexMC.root", puName)
    h_mc2 = aux.getFromFile(path+"nTrueVertexMC.root", "mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU")
    h_weight = aux.getFromFile(path+"data/puWeights.root", "pileupWeight_"+puName)
    h_weightUp = aux.getFromFile(path+"data/puWeights.root", "pileupWeightUp_"+puName)
    h_weightDn = aux.getFromFile(path+"data/puWeights.root", "pileupWeightDown_"+puName)

    for h in h_data, h_mc, h_mc2, h_weight:
        h.SetTitle(";;Events")

    aux.drawOpt(h_data, "data")
    h_mc.Scale(h_data.Integral())
    h_mc.SetLineColor(ROOT.kRed)
    h_mc2.Scale(h_data.Integral())
    h_mc2.SetLineColor(ROOT.kBlue)

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.leg.SetY1(.7)
    m.leg.SetX1(.6)
    m.leg.SetX2(.99)
    m.add(h_data, "Data")
    m.add(h_mc, "Simulation")
    m.add(h_mc2, "Signal simulation")
    m.Draw()

    # just create the pad, the ratio is not used
    r = ratio.Ratio( "Data/Sim.", h_data, h_mc )
    r.draw()

    h_weight.SetMaximum(2)
    h_weight.SetTitle(";number of interactions per event;Data/Sim.")
    h_weight.SetLineColor(h_mc.GetLineColor())
    h_weight.Draw("hist")

    hsys = aux.getSystFromDifference(h_weightUp, h_weightDn)
    aux.drawOpt(hsys, "sysUnc")
    hsys.Draw("same e2")


    l = aux.Label()
    can.SaveAs("plots/pu_weight.pdf")

