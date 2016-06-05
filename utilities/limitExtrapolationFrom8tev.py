#!/usr/bin/env python2

import ROOT
import re
from math import sqrt

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadRightMargin(0.01)
ROOT.gStyle.SetPadTopMargin(0.02)
ROOT.gStyle.SetLegendBorderSize(0)
ROOT.gROOT.SetBatch()

def readXsecSMS( name ):
    with open( name ) as f:
        lines = f.readlines()

    out = {}
    for line in lines:
        if line[0] == "#": continue
        m, x, uncrt = line.split(" ")
        m = int(m)
        x = float(x)
        uncrt = float(uncrt)
        out[m] = x
    return out

def getObj( fname, oname ):
    f = ROOT.TFile( fname )
    obj = f.Get( oname )
    obj = ROOT.gROOT.CloneObject( obj )
    return obj

def getLastObj( fname ):
    f = ROOT.TFile( fname )
    for i in f.GetListOfKeys():
        obj = i.ReadObj()
    obj = ROOT.gROOT.CloneObject( obj )
    return obj

def getGr( h_acc_multi, lumi, xsecs, b_for_20, b_uncert_rel, col=ROOT.kBlack, sigma=2 ):
    # b is normed to 20/fb

    h_acc = h_acc_multi[-1]

    h_limit = h_acc.Clone("limit")

    for x in range( h_limit.GetNbinsX()+2 ):
        for y in range( h_limit.GetNbinsY()+2 ):
            acc = h_acc.GetBinContent(x,y)
            mg = int( h_limit.GetXaxis().GetBinCenter(x) )
            mchi = int( h_limit.GetYaxis().GetBinCenter(y) )
            if mg not in xsecs:
                h_limit.SetBinContent(x,y,0)
                continue
            xsec = xsecs[mg]

            #if mg in [ 1100, 1500]:
            #    print mg, xsec, b, s, r*2

            b = b_for_20 *lumi/19712
            s = lumi*xsec*acc
            r = 1./sigma *s/sqrt(s+b+(b*b_uncert_rel)**2)

            #r = s/( 2+2*sqrt(1+b+(b*b_uncert_rel)**2) )


            h_limit.SetBinContent( x, y, r )

    h_limit.SetMaximum(2)
    h_limit.SetContour(2)
    c = ROOT.TCanvas()
    h_limit.Draw("cont z list")
    ROOT.gPad.Update()
    gr = ROOT.gROOT.GetListOfSpecials().FindObject("contours")[0].First()
    gr = ROOT.gROOT.CloneObject( gr )
    gr.SetLineColor( col )

    return gr

def interpolateAcceptance( hin ):

    newXmax = 2025
    xmin = hin.GetXaxis().GetXmin()
    xmax = hin.GetXaxis().GetXmax()
    xbins = hin.GetNbinsX()
    ymin = hin.GetYaxis().GetXmin()
    ymax = hin.GetYaxis().GetXmax()
    ybins = hin.GetNbinsY()

    newXbins = int((newXmax-xmin)/(xmax-xmin)*xbins)

    newYmax = newXmax - 25
    newYbins = int((newYmax-ymin)/(ymax-ymin)*ybins)


    hout = ROOT.TH2F("hacc2", "", newXbins, xmin, newXmax, newYbins, ymin, newYmax )

    for xbin in range(xbins+2):
        for ybin in range(ybins+2):
            x = hin.GetXaxis().GetBinCenter( xbin )
            y = hin.GetYaxis().GetBinCenter( ybin )

            bin = hout.FindBin( x, y )
            hout.SetBinContent( bin, hin.GetBinContent(xbin,ybin) )

    for xbin in range(18, newXbins+2):
        for ybin in range( xbin+16):
            hout.SetBinContent( xbin, ybin, 0.6)


    return hout


def interpolate( h ):


    for x in range(h.GetNbinsX()+2):
        for y in range(h.GetNbinsY()+2):
            if not h.GetBinContent(x,y):
                neighbourEntries = [
                    h.GetBinContent(x,y+1),
                    h.GetBinContent(x,y-1),
                    h.GetBinContent(x+1,y),
                    h.GetBinContent(x-1,y)
                ]
                filledNeighbourEntries = filter(lambda x: x != 0, neighbourEntries)
                if len(filledNeighbourEntries)>2:
                    mean = sum(filledNeighbourEntries) / len(filledNeighbourEntries)
                    h.SetBinContent( x,y, mean )


            if x > 17 and y < x+16:
                blub = h.GetBinContent(x-1,y-1)
                if y < x+8:
                    blub = h.GetBinContent(x-1,y)
                h.SetBinContent( x,y, blub )



    return h



def getAcceptances( name = "eventYieldT5gg-2015-01-09.txt" ):
    with open( name ) as f:
        lines = f.readlines()

    out = [ ROOT.TH2F("accbin6", ";M_{#tilde{g}} (GeV); M_{#tilde{#chi}^{0}_{1}} (GeV); Acceptance", 21, 775, 1825, 34, 0.0, 1700 ) ]*6
    #out = [ ROOT.TH2F("accbin6", ";M_{#tilde{g}} (GeV); M_{#tilde{#chi}^{0}_{1}} (GeV); Acceptance", 25, 775, 2025, 40, 0.0, 2000 ) ]*6
    for iline in range(0,len(lines), 12):
        mg = int(lines[iline].split("=")[-1])
        mc = int(lines[iline+1].split("=")[-1])
        ngen = int(lines[iline+2].split("=")[-1])
        yields = lines[iline+5].split("=")[-1]
        normYields = [ float(i)/ngen for i in yields.split(" ")[1:-1] ]
        out_bin = out[0].FindBin(mg, mc)
        for i, acc in enumerate(normYields):
            out[i].SetBinContent( out_bin, normYields[-1] )


    out = [ interpolate( h ) for h in out ]
    return out


original_gr = getObj("../../master/singlePhoton/plotTree/SMS_T5gg.root", "graph_Exp" )
original_gr.SetLineColor(2)
original_gr.SetLineWidth(2)
original_gr.SetLineStyle(2)
original_gr.SetMarkerColor(2)

xsecs13tev_gluglu = readXsecSMS( "../../master/infos/simplifiedModel_13TeV_GluGlu.xsec" )
xsecs8tev_gluglu = readXsecSMS( "../../master/infos/simplifiedModel.xsec" )

h_acc = getLastObj( "../../master/singlePhoton/PlotsSMS/config/SUS14004/2015-01-09-limits/SMS_T5gg/ROOT/SMS_T5gg_gluino_chi1_Acceptance.root" )
h_acc.UseCurrentStyle()
h_acc.SetTitle(";M_{#tilde{g}} (GeV);M_{#tilde{#chi^{0}_{1}}} (GeV);Acceptance" )
#h_acc = interpolateAcceptance( h_acc )
h_acc = getAcceptances()
h_acc[-1].Draw("colz")
ROOT.gPad.SaveAs("acceptance.pdf")

leg = ROOT.TLegend(.13, .13, .4, .5 )
leg.SetFillColor(0)
leg.SetHeader( "#scale[2]{T5gg}" )

bkg_20 = 29
rel_bkg_uncert = 9./29
bkg_20_13tev = 2.6 * bkg_20

c = ROOT.TCanvas()
h_acc[0].Draw("axis")
original_gr.Draw("l")
leg.AddEntry( original_gr, " 8 TeV 19.8fb^{-1} expected excl", "l" )

gr8_20 = getGr( h_acc, 19712, xsecs8tev_gluglu, bkg_20, rel_bkg_uncert, 2 )
c.cd()
#gr8_20.Draw("l")
#leg.AddEntry( gr8_20, "8TeV  20fb^{-1}", "l" )
grList = []

for lumi, col, sigma in [
        (2000.,ROOT.kBlue,2),
        (30000.,ROOT.kMagenta,2),
        (2e5,ROOT.kGreen+3,2),
        (2000.,ROOT.kBlue,3),
        (30000.,ROOT.kMagenta,3),
    ]:

    gr = getGr( h_acc, lumi, xsecs13tev_gluglu, bkg_20_13tev, rel_bkg_uncert, col, sigma )
    grList.append(gr)
    if sigma == 3:
        gr.SetLineStyle(3)
        gr.SetLineWidth(3)

        leg.AddEntry( gr, "13 TeV {:.0f}fb^{{-1}} 3#sigma sensitivity".format(lumi/1000), "l" )
    else:
        leg.AddEntry( gr, "13 TeV {:.0f}fb^{{-1}}".format(lumi/1000), "l" )
    c.cd()
    gr.Draw("l")

leg.Draw()

ROOT.gPad.SaveAs("limitInterpolation.pdf")

