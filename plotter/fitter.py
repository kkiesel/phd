#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("cFunctions/RooCMSShape_cc.so")
ROOT.gSystem.Load("cFunctions/ExpGaussExp_cc.so")

import ROOT
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

def getFromFile(fname, hname):
    f = ROOT.TFile(fname)
    x = f.Get(hname)
    x.SetDirectory(0)
    return x

def getHist(eff, num=True):
    return eff.GetCopyPassedHisto() if num else eff.GetCopyTotalHisto()


def getSigHist(hname="pt", num=True):
    f = ROOT.TFile("SingleElectron_Run2016-PromptReco_fake.root")
    eff = f.Get(hname)
    y = x.ProjectionX()
    y.SetDirectory(0)
    return y

def fitHist(name, h, infoText=""):
    var = "m_{e#gamma}" if "num" in name else "m_{ee+e#gamma}"
    x = ROOT.RooRealVar("x", var, 60, 120, "GeV")

    # signal - breit wigner
    zBosonMass = ROOT.RooRealVar("mZ","", 91.1876)
    zBosonWidth = ROOT.RooRealVar("BZ","", 2.4952)
    bw = ROOT.RooLandau("breitwigner","Breit Wigner", x, zBosonMass, zBosonWidth)

    # signal - smearing function
    p0 = ROOT.RooRealVar("p0", "p0", 90, 80, 100)
    p1 = ROOT.RooRealVar("p1", "p1", 8, 0, 20)
    p2 = ROOT.RooRealVar("p2", "p2", 3.5, 0, 10)
    p3 = ROOT.RooRealVar("p3", "p3", 3.5, 0, 10)
    signal_smear = ROOT.ExpGaussExp("signal_norm", "ExpGaussExp", x, p0, p1, p2, p3)
    x.setBins(10000, "cache")
    signal_norm = ROOT.RooFFTConvPdf("bwcb","Convolution", x, bw, signal_smear)

    nSig = ROOT.RooRealVar("nSig","number of events", 1000, 0, 20000)
    signal = ROOT.RooExtendPdf("signal", "", signal_norm, nSig)

    # background
    bkg_alpha = ROOT.RooRealVar("alpha", "alpha", 50, 0, 200)
    bkg_beta = ROOT.RooRealVar("beta", "beta", .02, 0, 20)
    bkg_gamma = ROOT.RooRealVar("gamma", "gamma", .05, 0, 10)
    bkg_peak = ROOT.RooRealVar("peak", "peak", 90, 0, 200)
    bkg_norm = ROOT.RooCMSShape("background_norm", "CMSShape", x, bkg_alpha, bkg_beta, bkg_gamma, bkg_peak)

    nBkg = ROOT.RooRealVar("nBkg", "number of events", 1000, 0, 100000)
    bkg = ROOT.RooExtendPdf("bkg", "", bkg_norm, nBkg)

    total = ROOT.RooAddPdf("total", "sig+bkg", ROOT.RooArgList(signal, bkg))

    # ranges
    x.setRange("belowZ", 60, 75)
    x.setRange("aboveZ", 105, 120)
    x.setRange("onZ", 80, 100)

    # import histogram
    dh = ROOT.RooDataHist("db", "db", ROOT.RooArgList(x), ROOT.RooFit.Import(h))

    # fit
    #signal.fitTo(dh, ROOT.RooFit.Range("onZ"))
    #bkg.fitTo(dh, ROOT.RooFit.Range("belowZ,aboveZ"))
    #total.fitTo(dh, ROOT.RooFit.Range("onZ"))

    # draw
    frame = x.frame(ROOT.RooFit.Title(" "))
    dh.plotOn(frame)
    #total.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
    #total.plotOn(frame, ROOT.RooFit.Components(ROOT.RooArgSet(bkg)), ROOT.RooFit.LineColor(ROOT.kGray), ROOT.RooFit.LineStyle(ROOT.kDashed))
    #total.plotOn(frame, ROOT.RooFit.Components(ROOT.RooArgSet(signal)), ROOT.RooFit.LineColor(ROOT.kGreen))

    c = ROOT.TCanvas()
    frame.Draw()
    l = aux.Label(info=infoText, sim=False)
    aux.save(name)
    return nSig.getVal()

fname = "../fakeRate/SingleElectron_Run2016-PromptReco_fake.root"

def inclusive():
    eff = getFromFile(fname, "vtx")
    num2d = getHist(eff)
    num1d = num2d.ProjectionX()
    nNum = fitHist("inclusive_num", num1d)
    den2d = getHist(eff, False)
    den1d = den2d.ProjectionX()
    nDen = fitHist("inclusive_den", den1d)
    print nNum, nDen

def binned(hname):
    eff = getFromFile(fname, hname)
    num2d = getHist(eff)
    den2d = getHist(eff, False)
    for bin in range(1, num2d.GetNbinsY()+1):
        yMin, yMax = num2d.GetYaxis().GetBinLowEdge(bin), num2d.GetYaxis().GetBinUpEdge(bin)
        infoText = "{} < {} < {}".format(yMin, num2d.GetYaxis().GetTitle(), yMax)
        num1d = num2d.ProjectionX("{}_{}_num".format(hname,bin), bin, bin)
        den1d = den2d.ProjectionX("{}_{}_den".format(hname,bin), bin, bin)
        nNum = fitHist("binned_{}_{}_num".format(bin,hname), num1d, infoText)
        nDen = fitHist("binned_{}_{}_den".format(bin,hname), den1d, infoText)
        print bin, nNum, nDen


inclusive()
binned("pt")

