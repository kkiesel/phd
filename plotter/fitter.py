#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("cFunctions/RooCMSShape_cc.so")
ROOT.gSystem.Load("cFunctions/ExpGaussExp_cc.so")

import ROOT
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
#ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.NumericIntegration)

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

def drawUnfitted(name, hData, hSig, hBkg, infoText=""):
    var = "m_{e#gamma}" if "num" in name else "m_{ee+e#gamma}"
    aux.drawOpt(hData, "data")
    hData.SetTitle(";{} (GeV);Pairs [abitrary normalization]".format(var))
    hData.Draw("ep")
    hSig.SetLineColor(ROOT.kRed)
    hSig.Scale(hData.GetMaximum()/hSig.GetMaximum())
    hSig.Draw("hist same")
    hBkg.SetLineColor(ROOT.kGreen)
    hBkg.Scale(hData.Integral()/hBkg.Integral())
    hBkg.Draw("hist same")
    l = aux.Label(info=infoText, sim=False)
    leg = ROOT.TLegend(.6, .6, .93, .92)
    leg.SetFillStyle(0)
    leg.AddEntry(hData, "Data")
    leg.AddEntry(hSig, "Signal shape")
    leg.AddEntry(hBkg, "Background shape")
    leg.Draw()
    aux.save(name+"_raw", log=False)

def fitHist(name, hData, hSig, hBkg, infoText=""):
    if not hData.GetEntries() or not hSig.GetEntries() or not hBkg.GetEntries(): return 0
    drawUnfitted(name, hData, hSig, hBkg, infoText)
    totalInt = hData.Integral(0,-1)
    var = "m_{e#gamma}" if "num" in name else "m_{ee+e#gamma}"
    x = ROOT.RooRealVar("x", var, 60, 120, "GeV")
    x.setRange("peak", 80, 100)

    dhSig = ROOT.RooDataHist("dhSig", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hSig))
    dhBkg = ROOT.RooDataHist("dhBkg", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hBkg))
    dh = ROOT.RooDataHist("dhData", "", ROOT.RooArgList(x), ROOT.RooFit.Import(hData))

    pdfSig = ROOT.RooHistPdf("histpdf1", "", ROOT.RooArgSet(x), dhSig, 2)
    pdfBkg = ROOT.RooHistPdf("histpdf2", "", ROOT.RooArgSet(x), dhBkg, 2)

    meanSig = ROOT.RooRealVar("meanSig", "meanSig", 0, -5, 5)
    widthSig = ROOT.RooRealVar("widthSig", "widthSig", 2, 0, 10)
    a1Sig = ROOT.RooRealVar("a1Sig", "", 1, 0, 25)
    a2Sig = ROOT.RooRealVar("a2Sig", "", 1, 0, 25)
    smearSig = ROOT.RooGaussian("gausSig", "", x, meanSig, widthSig)
    #smearSig = ROOT.ExpGaussExp("gausSig", "", x, meanSig, widthSig, a1Sig, a2Sig)

    meanBkg = ROOT.RooRealVar("meanBkg", "meanBkg", 0, -5, 5)
    widthBkg = ROOT.RooRealVar("widthBkg", "widthBkg", 2, 0, 10)
    smearBkg = ROOT.RooGaussian("gausBkg", "", x, meanBkg, widthBkg)

    x.setBins(10000, "cache")
    #x.setMin("cache", 10.3)
    #x.setMax("cache", 200.42)

    smearedSig = ROOT.RooFFTConvPdf("smearedSig","", x, pdfSig, smearSig)
    smearedBkg = ROOT.RooFFTConvPdf("smearedBkg","", x, pdfBkg, smearBkg)
    smearedBkg = ROOT.RooExponential("expo", "exponential PDF", x, meanBkg)

    nSig = ROOT.RooRealVar("nSig","number of events", totalInt, 0, 2*totalInt)
    signal = ROOT.RooExtendPdf("signal", "", smearedSig, nSig)

    nBkg = ROOT.RooRealVar("nBkg", "number of events", totalInt/10, 0, 2*totalInt)
    bkg = ROOT.RooExtendPdf("bkg", "", smearedBkg, nBkg)

    total = ROOT.RooAddPdf("total", "sig+bkg", ROOT.RooArgList(signal, bkg))

    # fit
    #signal.fitTo(dh)
    #bkg.fitTo(dh, ROOT.RooFit.Range("belowZ,aboveZ"))
    #total.fitTo(dh, ROOT.RooFit.Range("onZ"))
    total.fitTo(dh, ROOT.RooFit.Range("peak"))

    # draw
    frame = x.frame(ROOT.RooFit.Title(" "))
    dh.plotOn(frame)
    frame.SetYTitle(frame.GetYaxis().GetTitle().replace("Events", "Pairs"))
    total.plotOn(frame)
    #total.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
    total.plotOn(frame, ROOT.RooFit.Components("expo"), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed))
    total.plotOn(frame, ROOT.RooFit.Components("smearedSig"), ROOT.RooFit.LineColor(ROOT.kRed))

    #residuals = frame.residHist()

    #frameRes = x.frame(ROOT.RooFit.Title(" "))
    #frameRes.addPlotable(residuals, "P")

    c = ROOT.TCanvas()
    frame.Draw()
    #ratio.createBottomPad()
    #frameRes.Draw()
    l = aux.Label(info=infoText, sim=False)
    aux.save(name+"_fit", log=False)
    return nSig.getVal()


def binned(hname, isData=True):
    eff = getFromFile(fnameMC, hname)
    h2 = getHistFromEff(eff)
    hOut = h2.ProjectionY()
    for bin in range(1, h2.GetNbinsY()+1):
        ax = h2.GetYaxis()
        yMin, yMax = ax.GetBinLowEdge(bin), ax.GetBinUpEdge(bin)
        if not yMin - int(yMin): yMin = int(yMin)
        if not yMax - int(yMax): yMax = int(yMax)
        infoText = "{} < {} < {}".format(yMin, ax.GetTitle(), yMax)
        hNum = getHist(hname, data=isData, yBin1=bin, yBin2=bin)
        hNumMC = getHist(hname, data=False, yBin1=bin, yBin2=bin)
        hDen = getHist(hname, True, data=isData, yBin1=bin, yBin2=bin)
        hDenMC = getHist(hname, True, data=False, yBin1=bin, yBin2=bin)
        hBkg = getHist(hname+"_bkg", data=isData, yBin1=bin, yBin2=bin)
        num = fitHist("fakeRate_vs_{}_bin{}_num".format(hname,bin), hNum, hNumMC, hBkg, infoText)
        den = fitHist("fakeRate_vs_{}_bin{}_den".format(hname,bin), hDen, hDenMC, hBkg, infoText)
        if den:
            hOut.SetBinContent(bin, 100*num/den)
            hOut.SetBinError(bin, 100*math.sqrt(num)/den)
    hOut.SetTitle("")
    hOut.SetYTitle("f_{e#rightarrow#gamma} (%)")
    hOut.SetMaximum(5)
    hOut.SetMinimum(0)
    c = ROOT.TCanvas()
    hOut.Draw("hist e")
    aux.save("fakeRate_vs_{}".format(hname), log=False)

def inclusive():
    hNum = getHist("pt")
    hDen = getHist("pt", True)
    hBkg = getHist("pt_bkg")
    hNumMC = getHist("pt", data=False)
    hDenMC = getHist("pt", True, data=False)
    num = fitHist("fakeRate_inclusive_num", hNum, hNumMC, hBkg)
    den = fitHist("fakeRate_inclusive_den", hDen, hDenMC, hBkg)
    print "Inclusive fake-rate: {:.2f}%".format(100.*num/den)

inclusive()
binned("pt")
#binned("vtx")
#binned("jets")
#binned("met")
#binned("emht")
#binned("eta")
