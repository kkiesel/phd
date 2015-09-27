#!/usr/bin/env python2

import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filenames', nargs="+" )

args = parser.parse_args()

for fname in args.filenames:
    hname = fname.split('/')[-1].split(".")[0]
    ch = ROOT.TChain("TreeWriter/eventTree")
    ch.AddFile( fname )
    ch.Draw("mc_weight>>%s(1,0,1.01)"%hname, "1", "goff")
    h = ROOT.gROOT.Get( hname )
    print hname,
    print h.GetBinContent(0) / h.GetEntries()
