#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

fname = "../fakeRate/DYJetsToLL_M-50_ext_fake.root"

rfile = ROOT.TFile(fname)
for key in rfile.GetListOfKeys():
    name = key.GetName()
    eff = key.ReadObj()
    c = ROOT.TCanvas()
    info = ""
    if eff.GetDimension() == 1:
        c = ROOT.TCanvas()
        g = eff.CreateGraph()
        g.SetMaximum(0.05)
        g.SetMinimum(0)
        g.Draw("ap")
        if name.startswith("gen_"): info = "Single #gamma gen match"
    elif eff.GetDimension() == 2:
        h2Num = eff.GetCopyPassedHisto()
        h2Den = eff.GetCopyTotalHisto()
        hOut = h2Den.ProjectionY()
        hOut.Reset("ICES")
        for bin in range(0, h2Den.GetNbinsY()+2):
            hNum = h2Num.ProjectionX(aux.randomName(), bin, bin)
            hDen = h2Den.ProjectionX(aux.randomName(), bin, bin)
            xMin = hNum.GetXaxis().FindFixBin(60)
            xMax = hNum.GetXaxis().FindFixBin(120)
            xMin = hNum.GetXaxis().FindFixBin(-10) # temporary try out taking all
            xMax = hNum.GetXaxis().FindFixBin(1e8)
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
        style.defaultStyle()
        hOut.Draw("e hist")
        if "_gen" in name: info="Probe gen matched"
    elif eff.GetDimension() == 3:
        h3Num = eff.GetCopyPassedHisto()
        h3Den = eff.GetCopyTotalHisto()

        h2Out = h3Num.Project3D("yz")
        h2Out.Reset("ICES")
        for ybin in range(h3Num.GetNbinsY()+2):
            for zbin in range(h3Num.GetNbinsY()+2):
                hNum = h3Num.ProjectionX(aux.randomName(), ybin, ybin, zbin, zbin)
                hDen = h3Den.ProjectionX(aux.randomName(), ybin, ybin, zbin, zbin)
                xMin = hNum.GetXaxis().FindFixBin(60)
                xMax = hNum.GetXaxis().FindFixBin(120)
                num = hNum.Integral(xMin, xMax)
                den = hDen.Integral(xMin, xMax)
                if den:
                    h2Out.SetBinContent(ybin, zbin, 100*num/den)
                    h2Out.SetBinError(ybin, zbin, 100*math.sqrt(num)/den)
    else:
        print "Do not know what do do with {} dimensions".format(eff.GetDimension())
    l = aux.Label(info=info, sim=True)
    aux.save("fakeRate_DY_"+name, log=False)


