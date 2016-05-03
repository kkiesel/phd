#!/usr/bin/env python2

import ROOT
import glob
import re

ROOT.gROOT.SetBatch()



def h2s():
    title = ";m_{#tilde{g}} (GeV);m_{#tilde{#chi}^{0}_{1}} (GeV)"
    h = {}
    h["met"] = ROOT.TH2F("",title, 6, 900, 2100, 9, 100, 1900 )
    h["pt"] = ROOT.TH2F("",title, 6, 900, 2100, 9, 100, 1900 )
    h["eta"] = ROOT.TH2F("",title, 6, 900, 2100, 9, 100, 1900 )
    h["ht"] = ROOT.TH2F("",title, 6, 900, 2100, 9, 100, 1900 )
    h["nJet"] = ROOT.TH2F("",title, 6, 900, 2100, 9, 100, 1900 )


    for file in glob.glob("SMS-T5gg_mGluino-*_mNeutralino-*.root"):
        m = re.match("SMS-T5gg_mGluino-(\d+)_mNeutralino-(\d+).root", file)
        mg, mn = m.groups()
        mg, mn = float(mg), float(mn)

        ch = ROOT.TChain("signalTree")
        ch.AddFile( file )

        for name, hist in h.iteritems():
            if name == "eta": name = "abs(eta)"
            ch.Draw(name)
            hist.Fill( mg, mn, ROOT.htemp.GetMean() )
            del ROOT.htemp


    f = ROOT.TFile( "h2.root", "RECREATE");
    for name, hist in h.iteritems():
        hist.Write( name, ROOT.TObject.kWriteDelete )

def h1s(points):
    for file in glob.glob("SMS-T5gg_mGluino-*_mNeutralino-*.root"):



def main():
    h2s()
    h1s([(1600,200),(1600,1400),(2000,200),(2000,1800)])

main()
