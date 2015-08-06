#!/usr/bin/env python

# For more information, see
# Patrick Meade, Matthew Reece and David Shih:
# "Prompt Decays of General Neutralino NLSPs at the Tevatron"


import math
import ROOT

import style
ROOT.gStyle.SetPadTopMargin(0.08)
ROOT.gStyle.SetPadRightMargin(0.05)

ROOT.gROOT.SetBatch()


# natural constants
sw2 = 0.23120 # = sin^2(theta_weinberg)
cw2 = 1 - sw2
mz = 91.1876 # GeV
mh = 125. # GeV


def func_bToG( x ):
    rz = ( 1 - min([1,mz**2/x[0]**2]) )**4
    return cw2 / ( cw2 + sw2 * rz )

def func_bToZ( x ):
    rz = ( 1 - min([1,mz**2/x[0]**2]) )**4
    return sw2 * rz / ( cw2 + sw2 * rz )

def func_wToG( x ):
    rz = ( 1 - min([1,mz**2/x[0]**2]) )**4
    return sw2 / ( sw2 + cw2 * rz )

def func_wToZ( x ):
    rz = ( 1 - min([1,mz**2/x[0]**2]) )**4
    return cw2 * rz / ( sw2 + cw2 * rz )

xMin = 0
xMax = 1000

bToG = ROOT.TF1("bToG", func_bToG, xMin, xMax )
bToZ = ROOT.TF1("bToZ", func_bToZ, xMin, xMax )
wToG = ROOT.TF1("wToG", func_wToG, xMin, xMax )
wToZ = ROOT.TF1("wToZ", func_wToZ, xMin, xMax )

fsw2 = ROOT.TF1("sw2", str(sw2), xMin, xMax )
fcw2 = ROOT.TF1("cw2", str( 1-sw2 ), xMin, xMax )

for f in bToG, bToZ, wToG, wToZ:
    f.SetMinimum(0)
    f.SetMaximum(1.01)
    f.SetTitle(";m_{#tilde{#chi}^{0}_{1}} (GeV);#tilde{#chi}^{0}_{1} branching fraction       ")

for f in bToG, wToG, fcw2:
    f.SetLineColor( ROOT.kRed )

for f in bToZ, wToZ, fsw2:
    f.SetLineColor( ROOT.kBlue )

for f in fcw2, fsw2:
    f.SetLineColor(1)
    f.SetLineStyle(2)


bToG.Draw()
for f in bToZ, fcw2, fsw2:
    f.Draw("same")

weinbergScale=0.8

label = ROOT.TLatex()
label.DrawLatexNDC( .2, .67, "#scale[%s]{cos^{2}(#theta_{W})}"%weinbergScale )
label.DrawLatexNDC( .2, .33, "#scale[%s]{sin^{2}(#theta_{W})}"%weinbergScale )

label.DrawLatexNDC( .5, .8, "#color[2]{#tilde{#chi}^{0}_{1}#rightarrow#tilde{G} + #gamma}" )
label.DrawLatexNDC( .5, .2, "#color[4]{#tilde{#chi}^{0}_{1}#rightarrow#tilde{G} + Z}" )

label.DrawLatexNDC( .45, .94, "Bino-like #tilde{#chi}^{0}_{1}" )

ROOT.gPad.SaveAs("binoBranching.pdf")


wToG.Draw()
for f in wToZ, fcw2, fsw2:
    f.Draw("same")

label = ROOT.TLatex()
label.DrawLatexNDC( .3, .75, "#scale[%s]{cos^{2}(#theta_{W})}"%weinbergScale )
label.DrawLatexNDC( .3, .25, "#scale[%s]{sin^{2}(#theta_{W})}"%weinbergScale )

label.DrawLatexNDC( .65, .35, "#color[2]{#tilde{#chi}^{0}_{1}#rightarrow#tilde{G} + #gamma}" )
label.DrawLatexNDC( .65, .65, "#color[4]{#tilde{#chi}^{0}_{1}#rightarrow#tilde{G} + Z}" )

label.DrawLatexNDC( .45, .94, "Wino-like #tilde{#chi}^{0}_{1}" )

ROOT.gPad.SaveAs("winoBranching.pdf")




"""
def higgsinos( x, p ):
    M1, M2, tanb, eta = p
    sb = math.sin( math.atan( tanb ) )
    cb = math.cos( math.atan( tanb ) )


    bg = 0.5 * (sb + eta * cb )**2 * ( math.sqrt(sw2*cw2)*(M1 - M2) * mz / (M1*M2) )**2

    rz = ( 1 - min([1,mz**2/x[0]**2]) )**4
    bz = 0.5 * (sb + eta * cb )**2 * rz

    rh = ( 1 - min([1,mh**2/x[0]**2]) )**4
    bh = 0.5 * (sb - eta * cb )**2 * rh

    return bg, bz, bh


def func_hToG( x, p ):
    bg, bz, bh = higgsinos( x, p )
    return bg / ( bg + bz + bh )

def func_hToZ( x, p ):
    bg, bz, bh = higgsinos( x, p )
    return bz / ( bg + bz + bh )

def func_hToH( x, p ):
    bg, bz, bh = higgsinos( x, p )
    return bh / ( bg + bz + bh )

hToG = ROOT.TF1("hToG", func_hToG, xMin, xMax, 4 )
hToZ = ROOT.TF1("hToZ", func_hToZ, xMin, xMax, 4 )
hToH = ROOT.TF1("hToH", func_hToH, xMin, xMax, 4 )

for f in hToG, hToZ, hToH:
    f.SetMinimum(0)
    f.SetMaximum(1.01)
    f.SetParameters( 500, 1000, 1.5, -1 )
hToG.SetLineColor( ROOT.kBlue )
hToH.SetLineColor( ROOT.kMagenta )

hToG.Draw()
hToZ.Draw("same")
hToH.Draw("same")
ROOT.gPad.SaveAs("test.pdf")
"""

