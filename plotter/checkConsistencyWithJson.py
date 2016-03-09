#!/usr/bin/env python2
import sys
if sys.version_info[:2] == (2,6):
    print "Initialize correct python version first!"
    sys.exit()

import ROOT
import argparse
import re

# private libs
import style
import multiplot

import auxiliary as aux



from math import sqrt


def clearXaxisCurrentPad():
    # Delete label and title of all histograms in the current pad
    for ding in ROOT.gPad.GetListOfPrimitives():
        if isinstance( ding, ROOT.TH1 ) or isinstance( ding, ROOT.THStack ):
            xaxis = ding.GetXaxis()
            xaxis.SetLabelSize(0)
            xaxis.SetLabelColor(0)
            xaxis.SetLabelOffset(1000)
            xaxis.SetTitle("")
            xaxis.SetTitleColor(0)
            xaxis.SetTitleSize(0)

def createBottomPad( r=.5 ):
    # r: the ratio in which the pad is splitted
    ROOT.gPad.SetBottomMargin( r + (1-r)*ROOT.gPad.GetBottomMargin() - r*ROOT.gPad.GetTopMargin() )
    rPad = ROOT.TPad( "rPad", "ratio", 0, 0, 1, 1 )
    rPad.SetTopMargin( (1-r) - (1-r)*rPad.GetBottomMargin() + r*rPad.GetTopMargin() )
    rPad.SetFillStyle(3955)
    rPad.Draw()
    rPad.cd()
    rPad.SetLogy(0)
    return rPad

def convertToTH1( profile ):
    h = profile.ProjectionX( aux.randomName() )
    h.SetLineColor( profile.GetLineColor() )
    return h

class Ratio:
    def __init__( self, title, numerator, denominator, sysHisto=None ):

        # convcert TProfiles to histograms
        if isinstance( numerator, ROOT.TProfile ): numerator = convertToTH1( numerator )
        if isinstance( denominator, ROOT.TProfile ): denominator = convertToTH1( denominator )

        self.title = title
        self.numerator = numerator
        self.denominator = denominator
        self.sysHisto = sysHisto
        self.ratio = numerator.Clone( aux.randomName() )
        self.ratioStat = denominator.Clone( aux.randomName() )
        self.ratioSys = denominator.Clone( aux.randomName() )
        self.totalUncert = denominator.Clone( aux.randomName() )
        self.allowUnsymmetricYaxis = False

    def calculateRatio( self ):
        self.ratio.Divide( self.denominator )
        self.ratioStat.Divide( self.denominator )
        if self.sysHisto:
            self.ratioSys.Divide( self.denominator )
            self.totalUncert.Add( self.ratioSys, self.ratioStat, .5, .5 )


    def getYrange( self ):
        # If no minimum or maximum is specified, choose a minimum from 0 to .5
        # and a maximum from 1.5 to 50
        yMin = 0
        minimum = self.ratio.GetBinContent(self.ratio.GetMaximumBin())
        for bin in range( self.ratio.GetNbinsX()+2 ):
            minInBin = self.ratio.GetBinContent(bin)
            if minInBin < minimum and minInBin > 0:
                minimum = minInBin
        yMin = minimum*.95

        from math import ceil
        yMax = min( max(1.5, ceil(self.ratio.GetBinContent(self.ratio.GetMaximumBin()))), 50 )
        yMax = self.ratio.GetBinContent(self.ratio.GetMaximumBin())*1.05


        yValues = [ self.ratio.GetBinContent(bin) for bin in range( self.ratio.GetNbinsX()+2 ) ]
        yValues = filter( lambda a: a != 0, yValues )
        if self.allowUnsymmetricYaxis:
            return yMin, yMax
        else:
            yValuesAbsDiff = [ abs(x-1) for x in yValues ]
            yValuesAbsDiff.sort()
            maxYDiff = yValuesAbsDiff[-1]
            y = 0.5
            if maxYDiff < 0.05:
                y = 0.05
            if maxYDiff < 0.01:
                y = 0.01

            return 1-y, 1+y


    def draw( self, yMin=.9, yMax=1.1 ):
        self.calculateRatio()

        #yMin, yMax = self.getYrange()

        # Set ratio properties
        for hist in [ self.ratio, self.ratioSys, self.ratioStat, self.totalUncert ]:
            hist.GetYaxis().SetNdivisions(2, 5, 2)
            hist.SetTitleOffset(1.2, "Y")
            hist.SetYTitle( self.title )
            hist.SetMinimum( yMin )
            hist.SetMaximum( yMax )

        self.ratioSys.SetFillStyle(3254)
        self.ratioSys.SetMarkerSize(0)
        self.ratioSys.SetFillColor(46)

        self.totalUncert.SetFillStyle(3002)
        self.totalUncert.SetMarkerSize(0)
        self.totalUncert.SetMarkerStyle(0)
        self.totalUncert.SetFillColor( self.denominator.GetLineColor() )
        self.totalUncert.SetLineColor(0)

        self.ratioStat.SetLineWidth(5)
        self.ratioStat.SetMarkerStyle(0)
        self.ratioStat.SetLineColor(ROOT.kGray)

        clearXaxisCurrentPad()
        p = createBottomPad()

        self.totalUncert.Draw("e2")
        self.ratioStat.Draw("same e")
        if self.sysHisto:
            self.ratioSys.Draw("e2 same")
        self.ratio.Draw("e0 same")

        if yMin < 1 and yMax > 1:
            oneLine = ROOT.TLine()
            oneLine.SetLineStyle(2)
            axis = self.ratio.GetXaxis()
            #oneLine.DrawLine( axis.GetXmin(), 1.0, axis.GetXmax(), 1.0 )



def infoFromBrilcalcOutput( text ):
    lines = text.split("\n")
    out = {}
    for l in lines:
        p = l.split("|")
        if len(p)<2: continue
        m = re.match("\s*(\d+):(\d+)\s*", p[1] )
        if m: out[int(m.group(1))] = float(p[-2])
    return out

def getEventsPerRun( files ):
    chain = ROOT.TChain("TreeWriter/eventTree")
    for file in files: chain.AddFile( file )

    runs = {}
    for event in chain:
        if event.HLT_Photon90_CaloIdL_PFHT500_v:
            if event.runNo not in runs: runs[event.runNo] = 0
            runs[event.runNo] += 1

    return runs

def getLumisPerRun( json ):
    lumiText = ""#subprocess.check_output( "brilcalc lumi -u /pb --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/OfflineNormtagV2.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt", shell=True )
    if not lumiText:
        with open("bril.log") as f:
            lumiText = f.read()

    return infoFromBrilcalcOutput( lumiText )


def checkConsistencyWithJson( files, json ):

    #events = getEventsPerRun( files )
    events = {260627: 198576, 259626: 15159, 259637: 18285, 257461: 3926, 257613: 94231, 257614: 1095, 258157: 4965, 258129: 6279, 258403: 19722, 258136: 4727, 254906: 1846, 258655: 631, 258656: 33492, 259681: 2480, 259683: 9548, 259685: 67719, 259686: 32763, 257645: 78153, 258158: 134313, 258159: 32551, 256630: 1157, 259861: 7037, 258177: 132653, 258694: 19561, 259721: 14821, 258702: 37433, 258703: 39724, 254232: 114, 258706: 65973, 258712: 44177, 258713: 12887, 258714: 5290, 256673: 33, 256674: 24, 256675: 9293, 256676: 11655, 256677: 20032, 258214: 18968, 258215: 476, 258741: 5766, 258742: 75101, 258745: 26110, 257722: 1003, 257723: 7528, 258749: 55849, 258750: 17975, 257735: 728, 257751: 33919, 259809: 17791, 259810: 11864, 258427: 9786, 259813: 869, 259817: 364, 259818: 15875, 259820: 14957, 259821: 19796, 259822: 39777, 258287: 16640, 258213: 14292, 257804: 388, 257805: 21009, 260373: 13111, 259862: 54324, 254231: 39, 257816: 30021, 257819: 18782, 256801: 10992, 259884: 8008, 259890: 11267, 259891: 11335, 254790: 13319, 260424: 77719, 260425: 27532, 260426: 50669, 256843: 47047, 260431: 40944, 259811: 9229, 256866: 43, 256867: 5635, 256868: 27867, 256869: 1916, 258705: 9604, 257682: 16403, 258425: 12705, 258426: 881, 257599: 6250, 258428: 14458, 258432: 310, 258434: 37897, 254852: 1089, 258440: 55484, 258444: 2551, 258445: 20551, 258446: 9152, 258448: 44652, 256926: 1887, 254879: 2126, 256941: 10857, 257968: 19655, 257969: 48203, 260532: 78153, 260533: 1308, 260534: 36323, 260536: 16381, 260538: 25630, 254907: 1306, 256842: 23, 254914: 1186, 260427: 18582, 258211: 8093, 260575: 1982, 260576: 18273, 260577: 9008, 260593: 39803, 260541: 2087, 257531: 10565}
    lumis = getLumisPerRun( json )

    allRuns = sorted(lumis.keys() + list(set(events.keys()) - set(lumis.keys())))
    allInfo = {}

    h_lumi = ROOT.TH1F("","",len(lumis),0,len(lumis))
    h_data = ROOT.TH1F("","",len(lumis),0,len(lumis))

    for irun, run in enumerate(allRuns):
        e = events[run] if run in events else 0
        l = lumis[run] if run in lumis else 0
        allInfo[run] = (e,l)
        for h in h_lumi,h_data: h.GetXaxis().SetBinLabel(irun+1, str(run))
        h_lumi.SetBinContent(irun+1,l)
        h_lumi.SetBinError(irun+1,0.05*l)
        h_data.SetBinContent(irun+1,e)

    # draw
    c = ROOT.TCanvas("","",1200,400)

    h_lumi.SetLineColor(ROOT.kRed)

    h_lumi.Scale( h_data.Integral()/h_lumi.Integral() )

    h_lumi.Draw("hist")
    h_data.Draw("same ep")
    r = Ratio( "Data/Lumi", h_data, h_lumi )
    r.draw(0,2)

    aux.save("test")




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+', default=[] )
    parser.add_argument('--json', default="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt" )
    args = parser.parse_args()


    checkConsistencyWithJson( args.files, args.json )


if __name__ == "__main__":
    main()
