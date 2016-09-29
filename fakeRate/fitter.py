#!/usr/bin/env python2
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gSystem.Load("RooCMSShape_cc.so")
ROOT.gSystem.Load("ExpGaussExp_cc.so")

def getSigHist():
    f = ROOT.TFile("SingleElectron_Run2016D-PromptReco-v2_fake.root")
    h2 = f.Get("eg_sig_pt")
    h = h2.ProjectionX("px", 2,2)
    h = h2.ProjectionX()
    h.SetDirectory(0)
    return h

def main():
    x = ROOT.RooRealVar("x", "m_{ll} (GeV)", 60, 120)

    # signal
    # breit wigner
    zBosonMass = ROOT.RooRealVar("mZ","", 91.1876)
    zBosonWidth = ROOT.RooRealVar("BZ","", 2.4952)
    bw = ROOT.RooLandau("breitwigner","Breit Wigner", x, zBosonMass, zBosonWidth)

    p0 = ROOT.RooRealVar("p0", "p0", 90, 80, 100)
    p1 = ROOT.RooRealVar("p1", "p1", 8, 0, 20)
    p2 = ROOT.RooRealVar("p2", "p2", 3.5, 0, 10)
    p3 = ROOT.RooRealVar("p3", "p3", 3.5, 0, 10)

    signal_smear = ROOT.ExpGaussExp("signal_norm", "ExpGaussExp", x, p0, p1, p2, p3)
    x.setBins(10000,"cache")

    signal_norm = ROOT.RooFFTConvPdf("bwcb","Convolution", x, bw, signal_smear)
    nSig = ROOT.RooRealVar("nSig","number of events", 1000, 0, 20000)
    signal = ROOT.RooExtendPdf("signal", "", signal_smear, nSig)

    # background
    bkg_alpha = ROOT.RooRealVar("alpha", "alpha", 50, 0, 200)
    bkg_beta = ROOT.RooRealVar("beta", "beta", .02, 0, 20)
    bkg_gamma = ROOT.RooRealVar("gamma", "gamma", .05, 0, 10)
    bkg_peak = ROOT.RooRealVar("peak", "peak", 90, 0, 200)

    nBkg = ROOT.RooRealVar("nBkg", "number of events", 1000, 0, 100000)
    bkg_norm = ROOT.RooCMSShape("background_norm", "CMSShape", x, bkg_alpha, bkg_beta, bkg_gamma, bkg_peak)
    bkg = ROOT.RooExtendPdf("bkg", "", bkg_norm, nBkg)

    total = ROOT.RooAddPdf("total","sig+bkg", ROOT.RooArgList(signal, bkg))

    # ranges
    x.setRange("belowZ", 60, 75)
    x.setRange("aboveZ", 105, 120)
    x.setRange("onZ", 80, 100)

    # import histogram
    h = getSigHist()
    dh = ROOT.RooDataHist("db", "db", ROOT.RooArgList(x), ROOT.RooFit.Import(h))

    # fit
    #signal.fitTo(dh, ROOT.RooFit.Range("onZ"))
    #bkg.fitTo(dh, ROOT.RooFit.Range("belowZ,aboveZ"))
    total.fitTo(dh, ROOT.RooFit.Range("onZ"))

    # draw
    frame = x.frame()
    dh.plotOn(frame)
    total.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
    total.plotOn(frame, ROOT.RooFit.Components(ROOT.RooArgSet(bkg)), ROOT.RooFit.LineColor(ROOT.kGray), ROOT.RooFit.LineStyle(ROOT.kDashed))
    total.plotOn(frame, ROOT.RooFit.Components(ROOT.RooArgSet(signal)), ROOT.RooFit.LineColor(ROOT.kGreen))


    c = ROOT.TCanvas()
    frame.Draw()
    ROOT.gPad.SaveAs("test.pdf")
    print nSig.getVal(), nBkg.getVal()

main()


