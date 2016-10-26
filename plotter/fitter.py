#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("cFunctions/RooCMSShape_cc.so")
ROOT.gSystem.Load("cFunctions/ExpGaussExp_cc.so")

import ROOT
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setSilentMode(True)

fnameData = "../fakeRate/SingleElectron_Run2016-PromptReco_fake.root"
fnameMC = "../fakeRate/DYJetsToLL_M-50_ext_fake.root"

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
    #smearSig = ROOT.ExpGaussExp("gausSig", "", x, meanSig, widthSig, a1Sig, a2Sig) # double sided crystal ball

    meanBkg = ROOT.RooRealVar("meanBkg", "meanBkg", 0, -5, 5)
    widthBkg = ROOT.RooRealVar("widthBkg", "widthBkg", 2, 0, 10)
    smearBkg = ROOT.RooGaussian("gausBkg", "", x, meanBkg, widthBkg)

    smearedSig = ROOT.RooFFTConvPdf("smearedSig","", x, pdfSig, smearSig)
    smearedBkg = ROOT.RooFFTConvPdf("smearedBkg","", x, pdfBkg, smearBkg)

    nSig = ROOT.RooRealVar("nSig","number of events", totalInt, 0, 2*totalInt)
    signal = ROOT.RooExtendPdf("signal", "", smearedSig, nSig)

    nBkg = ROOT.RooRealVar("nBkg", "number of events", totalInt/10, 0, 2*totalInt)
    bkg = ROOT.RooExtendPdf("bkg", "", smearedBkg, nBkg)

    total = ROOT.RooAddPdf("total", "sig+bkg", ROOT.RooArgList(signal, bkg))

    # fit
    x.setBins(10000, "cache") # for the convolution
    x.setRange("peak", 60, 120) # range can be decreased to avoid boundary effects
    total.fitTo(dh, ROOT.RooFit.Range("peak"), ROOT.RooFit.PrintEvalErrors(-1))

    # draw
    frame = x.frame(ROOT.RooFit.Title(" "))
    dh.plotOn(frame)
    frame.SetYTitle(frame.GetYaxis().GetTitle().replace("Events", "Pairs"))
    total.plotOn(frame)
    total.plotOn(frame, ROOT.RooFit.Components("smearedBkg"), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed))
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


def binned(hname, selectionInfo="", intOnly=False):
    binning = None
    if hname.startswith("met"): binning = range(0,50,5)+range(50,120,10)
    if hname.startswith("pt"): binning = range(30, 90, 5) + [90,100,120,150,200]
    if hname.startswith("cIsoWorst"): binning = range(0, 8) + [8, 10, 15, 20]
    if hname.startswith("pIso"): binning = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.5, 2]
    if hname.startswith("hoe"): binning = list(aux.frange(0, 0.05, 0.005))
    effMC = aux.getFromFile(fnameMC, hname)
    h2NumMC = effMC.GetCopyPassedHisto()
    h2NumMC = aux.rebin2d(h2NumMC, None, binning)
    h2DenMC = effMC.GetCopyTotalHisto()
    h2DenMC = aux.rebin2d(h2DenMC, None, binning)
    eff = aux.getFromFile(fnameData, hname)
    h2Num = eff.GetCopyPassedHisto()
    h2Num = aux.rebin2d(h2Num, None, binning)
    h2Den = eff.GetCopyTotalHisto()
    h2Den = aux.rebin2d(h2Den, None, binning)
    effBkg = aux.getFromFile(fnameData, hname+"_bkg")
    h2Bkg = effBkg.GetCopyTotalHisto()
    h2Bkg = aux.rebin2d(h2Bkg, None, binning)

    hOut = h2Den.ProjectionY()
    hOut.Reset("ICES")
    for bin in range(0, h2Den.GetNbinsY()+2):
        ax = h2Den.GetYaxis()
        yMin, yMax = ax.GetBinLowEdge(bin), ax.GetBinUpEdge(bin)
        if not yMin - int(yMin): yMin = int(yMin)
        if not yMax - int(yMax): yMax = int(yMax)
        infoText = "{} < {} < {}".format(yMin, ax.GetTitle(), yMax)
        hNum = h2Num.ProjectionX(aux.randomName(), bin, bin)
        hNumMC = h2NumMC.ProjectionX(aux.randomName(), bin, bin)
        hDen = h2Den.ProjectionX(aux.randomName(), bin, bin)
        hDenMC = h2DenMC.ProjectionX(aux.randomName(), bin, bin)
        hBkg = h2Bkg.ProjectionX(aux.randomName(), bin, bin)
        if intOnly:
            xMin = hNum.GetXaxis().FindFixBin(80)
            xMax = hNum.GetXaxis().FindFixBin(100)
            num = hNum.Integral(xMin, xMax)
            den = hDen.Integral(xMin, xMax)
        else:
            num = fitHist("fakeRate_vs_{}_bin{}_num".format(hname,bin), hNum, hNumMC, hBkg, infoText)
            den = fitHist("fakeRate_vs_{}_bin{}_den".format(hname,bin), hDen, hDenMC, hBkg, infoText)
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
    l = aux.Label(info=selectionInfo)
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

#inclusive()
infos = {"":"", "_40pt":"p_{T}>40GeV", "_EB":"Barrel", "_EB_40pt":"Barrel, p_{T}>40GeV", "_EB_01eta_40pt":"Barrel, p_{T}>40GeV, |#eta|>0.1", "_EB_100pt":"Barrel, p_{T}>100GeV"}

variables = "pt", "jets", "met", "emht", "vtx", "eta", "sie", "sip", "hoe", "r9", "cIso", "nIso", "pIso", "cIsoWorst"
selections = "", "_40pt", "_EB", "_EB_40pt", "_EB_01eta_40pt", "_EB_100pt"
selections = ["_EB_40pt"]
variables = ["vtx"]

for var in variables:
    for sel in selections:
        binned("{}{}".format(var,sel), selectionInfo=infos[sel], intOnly=False)


