#!/usr/bin/env python2

import ROOT
import style
ROOT.gROOT.SetBatch()

diag = ROOT.TF1("diagonal", "x", 0, 900 )
diag_mt = ROOT.TF1("diag2", "x-173", 0, 900 )

l = ROOT.TLatex()

ms = "m_{#tilde{t}}"
mc = "m_{#tilde{#chi}_{1}^{0}}"

for f in diag, diag_mt:
    f.SetLineColor( ROOT.kBlack )
    f.SetTitle(";m_{#tilde{t}} (GeV);m_{#tilde{#chi}_{1}^{0}} (GeV)")
diag_mt.SetLineStyle(2)

diag.Draw()
diag_mt.Draw("same")

l.DrawLatexNDC( .3, .7, "forbidden: %s < %s"%(ms,mc) )

l.DrawLatexNDC( .58, .5, "T6bbxg" )
l.DrawLatexNDC( .7, .3, "T2ttxg" )

ROOT.gPad.SaveAs("stop_plane.pdf")
