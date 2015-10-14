import ROOT
import re
import math
from math import *

def getXsecInfoSMS( mother_mass, pklfilename ):
    import pickle

    info = 0
    with open( pklfilename, 'rb') as f:
        data = pickle.load( f )
        if mother_mass in data:
            info = data[mother_mass]
        else:
            print "could not find %s in "%mother_mass
    return info

def getXsecSMSglu( mother_mass ):
    return getXsecInfoSMS( mother_mass, "data/xSec_SMS_Gluino_13TeV.pkl" )[0]
def getXsecSMSsq( mother_mass ):
    return getXsecInfoSMS( mother_mass, "data/xSec_SMS_Squark_13TeV.pkl" )[0]
def getXsecSMSstop( mother_mass ):
    return getXsecInfoSMS( mother_mass, "data/xSec_SMS_StopSbottom_13TeV.pkl" )[0]

def getXsecFromName( name ):
    m = re.match( "[^_]*_(\d+)_[^\d]*.*", name )
    firstNumber = int( m.groups()[0] ) if m else 0

    if "T5" in name:
        return getXsecSMSglu( firstNumber )

    if "T2ttgg" in name:
        return getXsecSMSstop( firstNumber )

    print "cout not guess xsec from file name"
    return 0


def getFromFile( filename, histoname ):
    f = ROOT.TFile( filename )
    h = f.Get( histoname )
    if not h:
        print "Object {} not found in file {}".format(histoname, filename)
        return
    h = ROOT.gROOT.CloneObject( h )
    if isinstance( h, ROOT.TH1 ) and not h.GetSumw2N():
        h.Sumw2()
    h.drawOption_ = ""
    return h

def getNgen( filename ):
    h = getFromFile( filename, "hCutFlow" )
    return h.GetBinContent(1)

def getObjectNames( filename, path="", objects=[ROOT.TH1] ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )

    outList = []
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()

        if any( [ isinstance( obj, o ) for o in objects ] ):
            outList.append( element.GetName() )

    return outList

def rebin( h, binEdges ):
    import array
    binEdgesArr = array.array( 'd', binEdges )
    hnew = h.Rebin( len(binEdges)-1, "new", binEdgesArr )
    try: hnew.drawOption_ = h.drawOption_
    except: pass
    hnew.Scale( 1., "width" )
    return hnew

def absHistWeighted( origHist ):
    origNbins = origHist.GetNbinsX()
    origXmin = origHist.GetBinLowEdge(1)
    origXmax = origHist.GetBinLowEdge(origHist.GetNbinsX()+1)
    if origXmin + origXmax > 0.1:
        if origXmin:
            print "cant handle assymetric histograms"
        # else: print "already symetric?"
        return origHist

    h = origHist.Clone()
    newN = int(math.ceil(origNbins/2.))
    # TODO: doesn not work, only first bin in filled
    h = rebin( h, [ origXmax *i/newN for i in range(newN+1)] )

    for origBin in range( origNbins+2 ):
        newBin = int(abs(origBin - (origNbins+1.)/2)) + 1

        c1 = origHist.GetBinContent( origBin )
        e1 = origHist.GetBinError( origBin )
        c2 = h.GetBinContent( newBin )
        e2 = h.GetBinError( newBin )

        if e1 and e2:
            h.SetBinContent( newBin, ( c1*e1**-2 + c2*e2**-2 )/(e1**-2 + e2**-2) )
            h.SetBinError( newBin, 1./math.sqrt( e1**-2 + e2**-2 ) )

        else:
            h.SetBinContent( newBin, origHist.GetBinContent(origBin) )
            h.SetBinError( newBin, origHist.GetBinError(origBin) )

    return h

def randomName():
    """
    Generate a random string. This function is useful to give ROOT objects
    different names to avoid overwriting.
    """
    from random import randint
    from sys import maxint
    return "%x"%(randint(0, maxint))


def TH1F_binning( name, title, binEdges ):
    # Wrapper for the TH1F constructor for variable binning
    import array
    binEdgesArray = array.array( "d", binEdges )
    return ROOT.TH1F( name, ";%s;Events"%xVar.title, len(binEdgesArray)-1, binEdgesArray )

def createHistoFromTree( tree, variable, weight, nBins=20, xmin=0, xmax=0 ):
    name = randomName()
    if isinstance( nBins, list ):
        result = TH1F_binning( name, variable, nBins )
    else:
        result = ROOT.TH1F( name, variable, nBins, xmin, xmax )
    tree.Draw("%s>>%s"%(variable, name), weight, "goff")
    return result

def sumSq( *items ):
    return sqrt( sum( [i**2 for i in items ] ) )

def mergeBins( h, dest, source ):
    h.SetBinContent( dest, h.GetBinContent( dest ) + h.GetBinContent( source ) )
    h.SetBinError( dest, sumSq( h.GetBinError( dest ), h.GetBinError( source ) ) )

def appendFlowBin( h, under=True, over=True ):
    if under:
        mergeBins( h, 1, 0 )
    if over:
        mergeBins( h, h.GetNbinsX(), h.GetNbinsX()+1 )

def getValAndError( val, err, sig=2 ):
    from math import floor, log10
    digit = sig - int(floor(log10(err))) - 1
    return ( round(val,digit), round(err,digit) )

def getValAndErrorStr( val, err, sig=2 ):
    return "{} #pm {}".format( getValAndError( val, err, sig ) )


def getYAxisTitle( histo ):
    # returns e.g.: "Events / 10 GeV"
    yTitle = "Events"

    binW1 = histo.GetXaxis().GetBinWidth(1)
    binW2 = histo.GetXaxis().GetBinWidth(histo.GetNbinsX()+1)
    unit = "GeV" if "GeV" in histo.GetXaxis().GetTitle() else None

    if binW1 == binW2: #assume constant bin size
        if binW1 == 1:
            return yTitle

        # get two significant digits
        binW = getValAndError( 0, binW1 )[1]
        if binW.is_integer():
            binW = int(binW)
        if unit:
            return yTitle + " / " + str(binW) + " " + unit
        else:
            return yTitle + " / " + str(binW)
    else: # assume variable bin size
        if unit:
            return yTitle + " / " + unit
        else:
            return yTitle

def getROC( hSig, hBkg, highX=True ):
    # highX: signal is at high values of the variable

    nRocBins = hSig.GetNbinsX()

    sigEff = []
    bkgEff = []

    sigDen = hSig.Integral()
    bkgDen = hBkg.Integral()
    if not sigDen or not bkgDen:
        print "Warning, signal or background histogram has no integral"
        return

    for i in range(1, nRocBins+1 ):
        if  highX:
            sigNum = hSig.Integral(i, nRocBins+1)
            bkgNum = hBkg.Integral(i, nRocBins+1)
        else:
            sigNum = hSig.Integral(1, i)
            bkgNum = hBkg.Integral(1, i)

        sigEff.append( sigNum / sigDen )
        bkgEff.append( bkgNum / bkgDen )

    import numpy
    rocGraph = ROOT.TGraph( nRocBins, numpy.array(bkgEff), numpy.array(sigEff) )
    rocGraph.SetTitle(";#varepsilon_{bkg};#varepsilon_{sig}")
    return rocGraph


def automaticRebinner( hlist, maxEvents=0 ):
    # Ereates an array of bin edges on a list of histograms, such that each histogram
    # has at least 'maxEvents' events in each bin.
    out = []
    tmp = [0]*len(hlist)
    nBins = hlist[0].GetNbinsX()

    # check overflow bin
    overflowList = [ h.GetBinContent( nBins+1 ) for h in hlist ]
    if max(overflowList) > 1e-10:
        out.append( hlist[0].GetBinLowEdge(nBins+1) )

    for bin in range(nBins, -1, -1 ):

        for ih, h in enumerate(hlist):
            x = int(h.GetBinContent(bin))
            if x: print bin, h.GetBinLowEdge(bin), x
            else: print bin, h.GetBinLowEdge(bin)
            tmp[ih] += h.GetBinContent(bin)

        if max( tmp ) > maxEvents:
            if out:
                out.append( hlist[0].GetBinLowEdge(bin) )
            else:
                out.append( hlist[0].GetBinLowEdge(bin+1) )
            tmp = [0]*len(hlist)


    print out[::-1]


class Label:
    # Create labels
    # Usage:
    # * With Labels(), all default labels will be printed
    # * With Labels(False), the method is only initiated and labels can be modified before calling the 'draw' method

    cmsEnergy = 13 #TeV
    from main import intLumi

    def draw( self ):
        varDict = vars( self )
        for varName, obj in varDict.iteritems():
            if isinstance( obj, ROOT.TLatex ):
                obj.SetNDC()
                obj.Draw()

    def __init__( self, drawAll=True, sim=False, status="Private Work" ):
        # todo: include margins, etc
        if sim:
            self.cms = ROOT.TLatex( 0.2, .895, "#font[61]{CMS} #scale[0.76]{#font[52]{Simulation}}" )
        else:
            self.cms = ROOT.TLatex( 0.2, .895, "#font[61]{CMS}" )
        self.pub = ROOT.TLatex( 0.2, .865, "#scale[0.76]{#font[52]{%s}}"%status )
        self.lum = ROOT.TLatex( .63, .95, "%.2f fb^{-1} (%s TeV)"%(self.intLumi/1000., self.cmsEnergy) )

        if drawAll:
            self.draw()
