#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
if sys.version_info[:2] == (2,6):
    print("Initialize correct python version first!")
    sys.exit()

import ROOT
import math
import argparse
import re
from random import randint

# private libs
import ratio
import style
import multiplot

import auxiliary as aux

intLumi = 2.26e3 # https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/2577.html


def combinatoricalBkg():
    inFileName = "../histogramProducer/DYJetsToLL_M-50_fakeRate.root"
    effTot = aux.getFromFile( inFileName, "eta" )
    effTotBkg = aux.getFromFile( inFileName, "eta_bkg" )
    h2P = effTot.GetPassedHistogram()
    h2T = effTot.GetTotalHistogram()
    h2PBkg = effTotBkg.GetPassedHistogram()
    h2TBkg = effTotBkg.GetTotalHistogram()

    hP = h2P.ProjectionX()
    hT = h2T.ProjectionX()
    hPBkg = h2PBkg.ProjectionX("bkgP")
    hTBkg = h2TBkg.ProjectionX("bkgT")


    for h in hP, hPBkg: h.SetLineColor(ROOT.kRed)

    crDn, crUp = hP.FindFixBin(60), hP.FindFixBin(70)
    hPBkg.Scale( hP.Integral(crDn,crUp) / hPBkg.Integral(crDn,crUp) )
    hTBkg.Scale( hT.Integral(crDn,crUp) / hTBkg.Integral(crDn,crUp) )

    #for h in hP, hPBkg, hT, hTBkg: h.Scale(1./h.Integral())

    hT.Draw()
    hP.Draw("same")
    hPBkg.Draw("same hist")
    hTBkg.Draw("same hist")

    ROOT.gPad.SetLogy(1)
    aux.save("fakeRate_combinatoricalBkg")


"""To obtain the uncertainty on the fake rate coming from the fits, we fit the signal with a
Gaussian while modeling the background with an error function times an exponential. Also
we fit the signal with a Crystal Ball convolved with a Breit-Wigner and then fit the background
with an exponential function only.
"""
def makeFit( hist, name=None ):

    #return aux.integralAndError(hist, 80, 100, False)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    # variable
    x = ROOT.RooRealVar( "x", "x", 80., 100. )
    x.setBins(10000,"cache")
    x.setMin("cache",70.5)
    x.setMax("cache",110.5)

    # breit wiegner
    m0 = ROOT.RooRealVar( "m0", "m0", 91.1876 )
    width = ROOT.RooRealVar( "width", "width", 2.4952 )
    bw = ROOT.RooBreitWigner( "bw", "bw", x, m0, width )

    # crystal ball
    mean = ROOT.RooRealVar( "mean", "mean", 0., -1, 1 )
    sigma = ROOT.RooRealVar( "sigma", "sigma", 1.6, 0, 2 )
    alpha = ROOT.RooRealVar( "alpha", "alpha", 1.2, 0, 10 )
    n = ROOT.RooRealVar( "n", "n", 0.81, 0, 5 )
    cb = ROOT.RooCBShape( "cb", "cb", x, mean, sigma, alpha, n )


    # convolute breit wiegner with crystal ball
    sig = ROOT.RooFFTConvPdf( "sig", "sig", x, bw, cb )
    #sig = ROOT.RooVoigtian( "voigtian","vog", x, mean, width, sigma )

    # define background function
    a0 = ROOT.RooRealVar("a0","a0",10, 0, 1000)
    a1 = ROOT.RooRealVar("a1","a1",0,-2,0)
    a2 = ROOT.RooRealVar("a2","a2",1)
    bkg = ROOT.RooPolynomial("p2","p2",x,ROOT.RooArgList(a0,a1),0)


    # add signal + background
    nSig = ROOT.RooRealVar("nSig","nSig",0,1e5)
    nBkg = ROOT.RooRealVar("nBkg","nBkg",0,1e5)

    sig_ext = ROOT.RooExtendPdf("esig","extended signal p.d.f",sig,nSig,"signalRange")
    bkg_ext = ROOT.RooExtendPdf("ebkg","extended background p.d.f",bkg,nBkg,"signalRange")

    model = ROOT.RooAddPdf("model","sig+bkg",ROOT.RooArgList(sig_ext, bkg_ext))

    # import data
    dh = ROOT.RooDataHist("dh","dh",ROOT.RooArgList(x),ROOT.RooFit.Import(hist))

    # do the fit
    r = model.fitTo(dh, ROOT.RooFit.Save())
    #r.Print()

    plot = x.frame()
    plot.SetTitle(";m (GeV)")
    dh.plotOn( plot )
    model.plotOn( plot )
    plot.Draw()
    ROOT.gPad.SetLogy(0)
    if name: aux.save( "fakeRate_rooFit_{}".format(name) )
    return nSig.getVal(), nSig.getPropagatedError(r)



def fakeRate1d( inFileName, hname ):
    effTot = aux.getFromFile( inFileName, hname )
    h2P = effTot.GetPassedHistogram()
    h2T = effTot.GetTotalHistogram()

    hfakeP = h2P.ProjectionY("py1")
    hfakeT = h2T.ProjectionY("py2")
    hfakeP.Reset()
    hfakeT.Reset()

    for yBin in range(h2T.GetNbinsY()+2):
        hP = h2P.ProjectionX("px1", yBin, yBin)
        if hP.Integral()<1: continue
        cP, eP = makeFit( hP, "p_bin{}".format(yBin) )
        hfakeP.SetBinContent(yBin,cP)
        hfakeP.SetBinError(yBin,eP)

        hT = h2T.ProjectionX("px1", yBin, yBin)
        if hT.Integral()<1: continue
        cT, eT = makeFit( hT, "t_bin{}".format(yBin) )
        hfakeT.SetBinContent(yBin,cT)
        hfakeT.SetBinError(yBin,eT)

    c = ROOT.TCanvas()
    eff = ROOT.TGraphAsymmErrors()
    eff.BayesDivide( hfakeP, hfakeT )
    eff.Draw()
    aux.save("fakeRate_{}".format(hname))

def plotBoth( inFileName, hname ):
    effTot = aux.getFromFile( inFileName, hname )
    h2P = effTot.GetPassedHistogram()
    h2T = effTot.GetTotalHistogram()

    hP = h2P.ProjectionX("px1")
    hT = h2T.ProjectionX("px2")

    hP.SetLineColor(ROOT.kRed)

    hT.Draw()
    hT.SetMinimum(10)
    hP.Draw("same")

    ROOT.gPad.SetLogy(1)
    aux.save("fakeRate_both")

    c1, e1 = makeFit( hP, "test" )




def main():
    inFileName = "../histogramProducer/DYJetsToLL_M-50_fakeRate.root"
    #fakeRate1d( inFileName, "pt" )
    plotBoth( inFileName, "eta" )


if __name__ == "__main__":
    main()



