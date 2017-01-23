#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def ratioH2(h1, h2):
    h2.Scale(h1.Integral()/h2.Integral())
    c = ROOT.TCanvas()
    h1 = aux.rebin2d(h1, range(0,1000,50), [0.01*i for i in range(100)])
    h2 = aux.rebin2d(h2, range(0,1000,50), [0.01*i for i in range(100)])
    #h1.Divide(h2)
    for xbin in aux.iterate(h1):
        for ybin in aux.iterate(h1, "y"):
            c1 = h1.GetBinContent(xbin, ybin)
            c2 = h2.GetBinContent(xbin, ybin)
            e1 = max(1.14, h1.GetBinError(xbin,ybin))
            e2 = h2.GetBinError(xbin,ybin)
            h1.SetBinContent(xbin, ybin, (c1-c2)/math.sqrt(e1**2+e2**2))

    ROOT.gStyle.SetPalette(104)
    h1.Draw("colz")
    h1.SetMaximum(5)
    h1.SetMinimum(-5)
    aux.save("razor_ratio")

def getH2(selection="default"):
    h2 = None
    for f in dataHt.files:
        h = aux.getFromFile(f.replace("hists", "razor"), selection)
        if h2: h2.Add(h)
        else: h2 = h
    return h2

def razorFit(h):
    """ parameters:
    0 : normalization
    1 : b
    2 : mr0
    3 : r0
    4 : n
    """
    f = ROOT.TF2("razor", "[0]* ( [1]*pow((x[0]-[2])*(x[1]-[3]),1./[4]) -1 ) * exp( - [1]*pow((x[0]-[2])*(x[1]-[3]),1./[4]))", 800, 1000, 0.1, 2)
    f.SetParameters(h.Integral(), 3., 500, 0.1, 8)
    h.Fit("razor")


if __name__ == "__main__":
    h2data = getH2()
    razorFit(h2data)
    lnide
    h2jets = getH2("jControl")
    """

    style.style2d()
    c = ROOT.TCanvas()
    h2data.Draw("colz")
    aux.save("razor_sr")

    c = ROOT.TCanvas()
    h2jets.Draw("colz")
    aux.save("razor_cr")

    ratioH2(h2data, h2jets)
    """

    h2xBin = aux.rebin2d(h2data, range(0,1000,50), [0,.1,.2,.3,.4,.5])

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    for h in aux.getProjections(h2xBin):
        h.drawOption_ = "hist e"
        m.add(h, h.GetName())
    m.Draw()
    aux.save("projx", log=False)

    h2yBin = aux.rebin2d(h2data, range(0,1000,200), [0.01*i for i in range(100)] )

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    for h in aux.getProjections(h2yBin, "y"):
        h.drawOption_ = "hist e"
        m.add(h, h.GetName())
    m.Draw()
    aux.save("projy", log=False)

