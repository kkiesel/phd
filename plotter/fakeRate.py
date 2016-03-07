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
def makeFit( hist ):
    ROOT.gPad.SetLogy(0)

    # variable
    x = ROOT.RooRealVar( "x", "x", 00., 220. )
    x.setBins(10000,"cache")
    x.setMin("cache",70.5)
    x.setMax("cache",150.5)

    # breit wiegner
    m0 = ROOT.RooRealVar( "m0", "m0", 90.5 )
    width = ROOT.RooRealVar( "width", "width", 2.49 )
    bw = ROOT.RooBreitWigner( "bw", "bw", x, m0, width )

    # crystal ball
    mean = ROOT.RooRealVar( "mean", "mean", 0. )
    sigma = ROOT.RooRealVar( "sigma", "sigma", 2.6 )
    alpha = ROOT.RooRealVar( "alpha", "alpha", 1.2 )
    n = ROOT.RooRealVar( "n", "n", 0.81 )
    cb = ROOT.RooCBShape( "cb", "cb", x, mean, sigma, alpha, n )

    # convolute breit wiegner with crystal ball
    pdf = ROOT.RooFFTConvPdf( "pdf", "pdf", x, bw, cb )

    # define exponentional
    gamma = ROOT.RooRealVar( "gamma", "gamma", -.9 )
    exp = ROOT.RooExponential( "exp", "exp", x, gamma )

    # crystal ball
    mean_bk = ROOT.RooRealVar( "mean_bk", "mean_bk", 80 )
    sigma_bk = ROOT.RooRealVar( "sigma_bk", "sigma_bk", 200 )
    alpha_bk = ROOT.RooRealVar( "alpha_bk", "alpha_bk", 1.2 )
    n_bk = ROOT.RooRealVar( "n_bk", "n_bk", 0.81 )
    cb_bk = ROOT.RooCBShape( "cb_bk", "cb_bk", x, mean_bk, sigma_bk, alpha_bk, n_bk )



    # add signal + background
    sig1frac = ROOT.RooRealVar("sig1frac","fraction of component 1 in signal",0.8)
    #model = ROOT.RooAddPdf("model","sig+bkg",ROOT.RooArgList(exp,pdf), sig1frac)
    model = ROOT.RooAddPdf("model","sig+bkg",pdf, cb_bk, sig1frac)


    # import data
    dh = ROOT.RooDataHist("dh","dh",ROOT.RooArgList(x),ROOT.RooFit.Import(hist))

    # do the fit
    r = pdf.fitTo(dh, ROOT.RooFit.Save())
    r.Print()

    plot = x.frame()
    plot.SetTitle(";m (GeV)")
    dh.plotOn( plot )
    model.plotOn( plot )
    plot.Draw()
    aux.save( "fakeRate_test" )




def main():
    inFileName = "../histogramProducer/DYJetsToLL_M-50_fakeRate.root"
    effTot = aux.getFromFile( inFileName, "eta" )
    h2P = effTot.GetPassedHistogram()
    h2T = effTot.GetTotalHistogram()

    hP = h2P.ProjectionX()
    hT = h2T.ProjectionX()

    hP.SetLineColor(ROOT.kRed)

    hT.Draw()
    hT.SetMinimum(10)
    hP.Draw("same")

    ROOT.gPad.SetLogy(1)
    aux.save("fakeRate_both")

    makeFit( hP )



if __name__ == "__main__":
    main()



