#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

fname = "../fakeRate/DYJetsToLL_M-50_ext_fake.root"

rfile = ROOT.TFile(fname)
for key in rfile.GetListOfKeys():
    name = key.GetName()
    eff = key.ReadObj()
    if name.startswith("gen_"):
        c = ROOT.TCanvas()
        if eff.GetDimension()>1: continue
        g = eff.CreateGraph()
        g.SetMaximum(0.05)
        g.SetMinimum(0)
        g.Draw("ap")
        aux.save("fakeRate_"+name, log=False)
    else:
        h2Num = eff.GetCopyPassedHisto()
        h2Den = eff.GetCopyTotalHisto()

        hOut = h2Den.ProjectionY()
        hOut.Reset("ICES")
        for bin in range(0, h2Den.GetNbinsY()+2):
            ax = h2Den.GetYaxis()
            yMin, yMax = ax.GetBinLowEdge(bin), ax.GetBinUpEdge(bin)
            if not yMin - int(yMin): yMin = int(yMin)
            if not yMax - int(yMax): yMax = int(yMax)
            infoText = "{} < {} < {}".format(yMin, ax.GetTitle(), yMax)
            hNum = h2Num.ProjectionX(aux.randomName(), bin, bin)
            hDen = h2Den.ProjectionX(aux.randomName(), bin, bin)
            xMin = hNum.GetXaxis().FindFixBin(60)
            xMax = hNum.GetXaxis().FindFixBin(120)
            num = hNum.Integral(xMin, xMax)
            den = hDen.Integral(xMin, xMax)
            if den:
                hOut.SetBinContent(bin, 100*num/den)
                hOut.SetBinError(bin, 100*math.sqrt(num)/den)
        hOut.SetTitle("")
        hOut.SetYTitle("f_{e#rightarrow#gamma} (%)")
        hOut.SetTitleOffset(0.9,"y")
        aux.drawOpt(hOut, "data")
        hOut.SetMaximum(5)
        hOut.SetMinimum(0)
        c = ROOT.TCanvas()
        style.defaultStyle()
        hOut.Draw("e hist")
        l = aux.Label(info="")
        aux.save("fakeRate_gen_vs_{}".format(name), log=False)


