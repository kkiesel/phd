#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("cFunctions/RooCMSShape_cc.so")
ROOT.gSystem.Load("cFunctions/ExpGaussExp_cc.so")

import ROOT
#ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


fnameData = "../fakeRate/SingleElectron_Run2016-PromptReco_fake.root"
fnameMC = "../fakeRate/DYJetsToLL_M-50_ext_fake.root"


def getFromFile(fname, hname):
    f = ROOT.TFile(fname)
    x = f.Get(hname)
    x.SetDirectory(0)
    return x

def getHist(eff, num=True):
    return eff.GetCopyPassedHisto() if num else eff.GetCopyTotalHisto()

def fitHist(name, hData, hSig, hBkg, infoText=""):
    totalInt = hData.Integral(0,-1)
    var = "m_{e#gamma}" if "num" in name else "m_{ee+e#gamma}"
    x = ROOT.RooRealVar("x", var, 60, 120, "GeV")

    # signal - breit wigner
    zBosonMass = ROOT.RooRealVar("mZ","", 91.1876)
    zBosonWidth = ROOT.RooRealVar("BZ","", 2.4952)
    bw = ROOT.RooLandau("breitwigner","Breit Wigner", x, zBosonMass, zBosonWidth)

    dhMC = ROOT.RooDataHist("mc", "mc", ROOT.RooArgList(x), ROOT.RooFit.Import(hMC))
    bw = ROOT.RooHistPdf("histpdf1", "histpdf1", ROOT.RooArgSet(x), dhMC, 0)

    # signal - smearing function
    p0 = ROOT.RooRealVar("p0", "p0", 0, -5, 5)
    p1 = ROOT.RooRealVar("p1", "p1", 8, 0, 20)
    p2 = ROOT.RooRealVar("p2", "p2", 3.5, 0, 10)
    p3 = ROOT.RooRealVar("p3", "p3", 3.5, 0, 10)
    signal_smear = ROOT.ExpGaussExp("signal_norm", "ExpGaussExp", x, p0, p1, p2, p3)

    signal_smear = ROOT.RooGaussian("gaus", "gaus", x, p0, p1)


    x.setBins(10000, "cache")
    signal_norm = ROOT.RooFFTConvPdf("bwcb","Convolution", x, bw, signal_smear)

    nSig = ROOT.RooRealVar("nSig","number of events", totalInt, 0, 2*totalInt)
    signal = ROOT.RooExtendPdf("signal", "", signal_norm, nSig)

    # background
    bkg_alpha = ROOT.RooRealVar("alpha", "alpha", 0, -200, 200)
    bkg_beta = ROOT.RooRealVar("beta", "beta", .02, 0, 20)
    bkg_gamma = ROOT.RooRealVar("gamma", "gamma", .05, 0, 10)
    bkg_peak = ROOT.RooRealVar("peak", "peak", 90, 0, 200)
    bkg_norm = ROOT.RooCMSShape("background_norm", "CMSShape", x, bkg_alpha, bkg_beta, bkg_gamma, bkg_peak)

    bkg_norm = ROOT.RooExponential("exp","", x, bkg_alpha)

    nBkg = ROOT.RooRealVar("nBkg", "number of events", totalInt/10, 0, 2*totalInt)
    bkg = ROOT.RooExtendPdf("bkg", "", bkg_norm, nBkg)

    total = ROOT.RooAddPdf("total", "sig+bkg", ROOT.RooArgList(signal, bkg))

    # ranges
    x.setRange("belowZ", 60, 75)
    x.setRange("aboveZ", 105, 120)
    x.setRange("onZ", 80, 100)

    # import histogram
    dh = ROOT.RooDataHist("db", "db", ROOT.RooArgList(x), ROOT.RooFit.Import(hData))

    # fit
    #signal.fitTo(dh, ROOT.RooFit.Range("onZ"))
    #bkg.fitTo(dh, ROOT.RooFit.Range("belowZ,aboveZ"))
    #total.fitTo(dh, ROOT.RooFit.Range("onZ"))
    total.fitTo(dh)

    # draw
    frame = x.frame(ROOT.RooFit.Title(" "))
    dh.plotOn(frame)
    signal.plotOn(frame)
    #total.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
    #total.plotOn(frame, ROOT.RooFit.Components(ROOT.RooArgSet(bkg)), ROOT.RooFit.LineColor(ROOT.kGray), ROOT.RooFit.LineStyle(ROOT.kDashed))
    #total.plotOn(frame, ROOT.RooFit.Components(ROOT.RooArgSet(signal)), ROOT.RooFit.LineColor(ROOT.kGreen))

    residuals = frame.residHist()

    frameRes = x.frame(ROOT.RooFit.Title(" "))
    frameRes.addPlotable(residuals, "P")

    c = ROOT.TCanvas()
    frame.Draw()
    ratio.createBottomPad()
    frameRes.Draw()
    l = aux.Label(info=infoText, sim=False)
    aux.save("fakeRate_peaks_"+name, log=False)
    return nSig.getVal()

def getHist(name, den=False, data=True, yBin1=0, yBin2=-1):
    fname = fnameData if data else fnameMC
    f = ROOT.TFile(fname)
    eff = f.Get(name)
    h2 = eff.GetCopyTotalHisto() if den else eff.GetCopyPassedHisto()
    h1 = h2.ProjectionX(aux.randomName(), yBin1, yBin2)
    h1.SetDirectory(0)
    return h1


hData = getHist("pt")
hBkg = getHist("pt_bkg", True, data=False)
hSig = getHist("pt", data=False)

fitHist("inclusive", hData, hSig, hBkg)

"""
def inclusive():
    for i in range(2,12):
        eff = getFromFile(fnameData, "pt")
        num2d = getHist(eff, False)
        num1d = num2d.ProjectionX("pt1", i, i)

        effMC = getFromFile(fnameMC, "pt_bkg")
        num2dMC = getHist(effMC, False)
        num1dMC = num2dMC.ProjectionX("ldien")
        num1dMC.Draw()
        aux.save("test")

        effMC2 = getFromFile(fnameMC, "pt")
        num2dMC2 = getHist(effMC2, False)
        num1dMC2 = num2dMC2.ProjectionX("ilnde", i, i)

        num1d.Draw("e")
        if num1dMC.GetMaximum():
            num1dMC.Scale(num1d.GetMaximum()/num1dMC.GetMaximum())
        num1dMC.Draw("hist same")

        if num1dMC2.GetBinContent(1):
            num1dMC2.Scale(num1d.GetBinContent(1)/num1dMC2.GetBinContent(1))
        num1dMC2.Draw("same")


        aux.save("test_%s"%i)

        #nNum = fitHist("inclusive_num", num1d, num1dMC)



    #den2d = getHist(eff, False)
    #den1d = den2d.ProjectionX()
    #nDen = fitHist("inclusive_den", den1d)
    #print nNum, nDen

def binned(hname):
    eff = getFromFile(fname, hname)
    num2d = getHist(eff)
    den2d = getHist(eff, False)
    for bin in range(1, num2d.GetNbinsY()+1):
        yMin, yMax = num2d.GetYaxis().GetBinLowEdge(bin), num2d.GetYaxis().GetBinUpEdge(bin)
        if not yMin - int(yMin): yMin = int(yMin)
        if not yMax - int(yMax): yMax = int(yMax)
        infoText = "{} < {} < {}".format(yMin, num2d.GetYaxis().GetTitle(), yMax)
        num1d = num2d.ProjectionX("{}_{}_num".format(hname,bin), bin, bin)
        den1d = den2d.ProjectionX("{}_{}_den".format(hname,bin), bin, bin)
        if not num1d.Integral() or not den1d.Integral(): continue
        nNum = fitHist("binned_{}_{}_num".format(bin,hname), num1d, infoText)
        nDen = fitHist("binned_{}_{}_den".format(bin,hname), den1d, infoText)
        print bin, nNum, nDen


#inclusive()
#binned("pt")
#binned("vtx")
#binned("jets")
#binned("met")
#binned("emht")

"""
