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

def getHistFromEff(eff, num=True):
    return eff.GetCopyPassedHisto() if num else eff.GetCopyTotalHisto()

def getHist(name, den=False, data=True, yBin1=0, yBin2=-1):
    fname = fnameData if data else fnameMC
    f = ROOT.TFile(fname)
    eff = f.Get(name)
    h2 = eff.GetCopyTotalHisto() if den else eff.GetCopyPassedHisto()
    h1 = h2.ProjectionX(aux.randomName(), yBin1, yBin2)
    h1.SetDirectory(0)
    return h1



def fitHist(name, hData, hSig, hBkg, infoText=""):
    totalInt = hData.Integral(0,-1)
    var = "m_{e#gamma}" if "num" in name else "m_{ee+e#gamma}"
    x = ROOT.RooRealVar("x", var, 60, 120, "GeV")

    dhSig = ROOT.RooDataHist("dhSig", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hSig))
    dhBkg = ROOT.RooDataHist("dhBkg", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hBkg))
    dh = ROOT.RooDataHist("dhData", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hData))

    pdfSig = ROOT.RooHistPdf("histpdf1", "", ROOT.RooArgSet(x), dhSig, 0)
    pdfBkg = ROOT.RooHistPdf("histpdf2", "", ROOT.RooArgSet(x), dhBkg, 0)

    meanSig = ROOT.RooRealVar("meanSig", "meanSig", 0, -5, 5)
    widthSig = ROOT.RooRealVar("widthSig", "widthSig", 2, 0, 10)
    a1Sig = ROOT.RooRealVar("a1Sig", "", 1, 1, 25)
    a2Sig = ROOT.RooRealVar("a2Sig", "", 1, 1, 25)
    smearSig = ROOT.RooGaussian("gausSig", "", x, meanSig, widthSig)
    smearSig = ROOT.ExpGaussExp("gausSig", "", x, meanSig, widthSig, a1Sig, a2Sig)

    meanBkg = ROOT.RooRealVar("meanBkg", "meanBkg", 0, -5, 5)
    widthBkg = ROOT.RooRealVar("widthBkg", "widthBkg", 2, 0, 10)
    smearBkg = ROOT.RooGaussian("gausBkg", "", x, meanBkg, widthBkg)

    x.setBins(10000, "cache")
    #x.setMin("cache", 40)
    #x.setMax("cache", 130)

    smearedSig = ROOT.RooFFTConvPdf("smearedSig","", x, pdfSig, smearSig)
    smearedBkg = ROOT.RooFFTConvPdf("smearedBkg","", x, pdfBkg, smearBkg)

    nSig = ROOT.RooRealVar("nSig","number of events", totalInt, 0, 2*totalInt)
    signal = ROOT.RooExtendPdf("signal", "", smearedSig, nSig)

    nBkg = ROOT.RooRealVar("nBkg", "number of events", totalInt/10, 0, 2*totalInt)
    bkg = ROOT.RooExtendPdf("bkg", "", smearedBkg, nBkg)

    total = ROOT.RooAddPdf("total", "sig+bkg", ROOT.RooArgList(signal, bkg))

    # fit
    #signal.fitTo(dh, ROOT.RooFit.Range("onZ"))
    #bkg.fitTo(dh, ROOT.RooFit.Range("belowZ,aboveZ"))
    #total.fitTo(dh, ROOT.RooFit.Range("onZ"))
    #total.fitTo(dh)

    # draw
    frame = x.frame(ROOT.RooFit.Title(" "))
    dh.plotOn(frame)
    #total.plotOn(frame)
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


hData = getHist("pt")
hBkg = getHist("pt_bkg")
hSig = getHist("pt", data=False)
#fitHist("inclusive", hData, hSig, hBkg)

hData.Draw("ep")
hSig.Scale(hData.GetMaximum()/hSig.GetMaximum())
hSig.Draw("hist same")
hBkg.Scale(hData.Integral()/hBkg.Integral()/10)
hBkg.Draw("hist same")
aux.save("test", log=False)

def binned(hname):
    eff = getFromFile(fnameMC, hname)
    h2 = getHistFromEff(eff)
    hOut = h2.ProjectionY()
    for bin in range(1, h2.GetNbinsY()+1):
        ax = h2.GetYaxis()
        yMin, yMax = ax.GetBinLowEdge(bin), ax.GetBinUpEdge(bin)
        if not yMin - int(yMin): yMin = int(yMin)
        if not yMax - int(yMax): yMax = int(yMax)
        infoText = "{} < {} < {}".format(yMin, ax.GetTitle(), yMax)
        hDataNum = getHist(hname, data=False, yBin1=bin, yBin2=bin)
        hDataDen = getHist(hname, True, data=False, yBin1=bin, yBin2=bin)
        r1, r2 = hDataNum.FindFixBin(80), hDataNum.FindFixBin(100)
        num = hDataNum.Integral(r1, r2)
        den = hDataDen.Integral(r1, r2)
        if den:
            hOut.SetBinContent(bin, 100*num/den)
            hOut.SetBinError(bin, 100*math.sqrt(num)/den)
    hOut.SetTitle("")
    hOut.SetYTitle("f_{e#rightarrow#gamma} (%)")
    hOut.SetMaximum(5)
    hOut.SetMinimum(0)
    hOut.Draw("hist e")
    aux.save("fakeRate_vs_{}".format(hname), log=False)

#inclusive()
binned("pt")
binned("vtx")
binned("jets")
binned("met")
binned("emht")
binned("eta")

