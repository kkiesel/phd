import ROOT
import auxiliary as aux
from math import sqrt
import auxiliary as aux


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

def createBottomPad( r=.2 ):
    # r: the ratio in which the pad is splitted
    ROOT.gPad.SetBottomMargin( r + (1-r)*ROOT.gPad.GetBottomMargin() - r*ROOT.gPad.GetTopMargin() )
    rPad = ROOT.TPad( "rPad", "ratio", 0, 0, 1, 1 )
    rPad.SetTopMargin( (1-r) - (1-r)*rPad.GetBottomMargin() + r*rPad.GetTopMargin() )
    rPad.SetFillStyle(3955)
    rPad.Draw()
    rPad.cd()
    rPad.SetLogy(0)

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
        createBottomPad()

        self.totalUncert.Draw("e2")
        self.ratioStat.Draw("same e")
        if self.sysHisto:
            self.ratioSys.Draw("e2 same")
        self.ratio.Draw("e0 same")

        if yMin < 1 and yMax > 1:
            oneLine = ROOT.TLine()
            oneLine.SetLineStyle(2)
            axis = self.ratio.GetXaxis()
            oneLine.DrawLine( axis.GetXmin(), 1.0, axis.GetXmax(), 1.0 )


