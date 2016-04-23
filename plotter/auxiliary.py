# -*- coding: utf-8 -*-
import ROOT
import re
import math
from math import *
import ConfigParser

binCfg = ConfigParser.SafeConfigParser()
binCfg.readfp(open('rebin.cfg'))

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
    if not m: m = re.match( ".*mGluino-(\d+)_.*", name )
    firstNumber = int( m.groups()[0] ) if m else 0

    if "T5" in name:
        return getXsecSMSglu( firstNumber )

    if "T2ttgg" in name:
        return getXsecSMSstop( firstNumber )

    print "cout not guess xsec from file name"
    return 0

def getSignalLabel( dataset ):
    m = re.match( "([^_]+)_(\d+)_(\d+)", dataset )
    if m:
        return "{} m(#tilde{{g}})={} m(#tilde{{#chi}}^{{0}}_{{1}})={}".format(*m.groups())
    return dataset

def getDatasetFromKey(key):
    m = re.match( "T([^_]+)_(\d+)_(\d+)", key )
    if m: return "SMS-T{}_mGluino-{}_mNeutralino-{}".format(*m.groups())
    return key


def getProjections( h2, alongX=True ):
    hs = []
    label = h2.GetYaxis().GetTitle()

    for ybin in range( h2.GetNbinsY()+2 ):

        ylow = h2.GetYaxis().GetBinLowEdge(ybin)
        yhigh = h2.GetYaxis().GetBinUpEdge(ybin)
        name = "{} #leq {} < {}".format( ylow, label, yhigh )
        if ybin == 0: name = "{} < {}".format( label, yhigh )
        if ybin == h2.GetNbinsY()+1: name = "{} #leq {}".format( ylow, label )


        h = h2.ProjectionX( name, ybin, ybin )
        h.SetLineColor( ybin+2)
        if h.GetEntries():
            h.Scale( 1./h.GetEntries() )
            hs.append( h )

    return hs


def write2File( obj2Write, name, fname ):
    obj = obj2Write.Clone()
    if isinstance( obj, ROOT.TH1 ):
        for a in obj.GetXaxis(),obj.GetYaxis():
            pass
            #a.UnZoom() #TODO: find out why this affects original histogram as well

    f = ROOT.TFile( fname,"update")
    obj.Write( name, ROOT.TObject.kWriteDelete )
    f.Close()

def writeWeight( obj, name, sampleName ):
    print "Info in <auxiliary::writeWeight>: {} written to file {}".format(name, "weights.root" )
    write2File( obj, name.replace("h_","weight_{}_".format(sampleName)), "weights.root" )

def getFromFile( filename, histoname ):
    f = ROOT.TFile( filename )
    h = f.Get( histoname )
    if not h:
        if ROOT.gErrorIgnoreLevel < ROOT.kBreak:
            print "Object {} not found in file {}".format(histoname, filename)
        return
    h = ROOT.gROOT.CloneObject( h )
    if isinstance( h, ROOT.TH1 ) and not h.GetSumw2N():
        h.Sumw2()
    h.drawOption_ = ""
    return h

def getNgen( filename ):
    h = getFromFile( filename, "hCutFlow" )
    return h.GetBinContent(2)

def getObjectNames( filename, path="", objects=[ROOT.TH1] ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )

    outList = []
    for element in tmpDir.GetListOfKeys():
        if element.GetName() == "rawEff_vs_run": continue
        obj = element.ReadObj()

        if any( [ isinstance( obj, o ) for o in objects ] ):
            outList.append( element.GetName() )

    return outList

def getBinning( axis ):
    binning = []
    for i in range(axis.GetNbins()+1):
        binning.append( axis.GetBinUpEdge(i) )
    return binning

def checkRebinningConsistence( axis, newBinning ):
    oldBinning = getBinning( axis )
    # get rid of unprecise floats:
    oldBinning = [ round(i,5) for i in oldBinning ]
    newBinning = [ round(i,5) for i in newBinning ]
    # ignore new bin edges out of range of old binning
    newBinning = [ i for i in newBinning if i>= oldBinning[0] ]
    for i in newBinning:
        if i not in oldBinning: print "New bin edge is not compatible with old binning", i, "old binning:", oldBinning

def rebin2d( h, binEdgesX, binEdgesY ):
    # Check consistency with old binning
    checkRebinningConsistence( h.GetXaxis(), binEdgesX )
    checkRebinningConsistence( h.GetYaxis(), binEdgesY )

    # Create
    import array
    binEdgesXArr = array.array( 'd', binEdgesX )
    binEdgesYArr = array.array( 'd', binEdgesY )
    hnew = ROOT.TH2F(h.GetName(),h.GetTitle(), len(binEdgesX)-1, binEdgesXArr, len(binEdgesY)-1, binEdgesYArr )

    # GetProperties
    hnew.drawOption_ = h.drawOption_ if hasattr( h, "drawOption_" ) else ""
    hnew.SetTitle("{};{};{}".format(h.GetTitle(),h.GetXaxis().GetTitle(),h.GetYaxis().GetTitle()))

    # Fill
    for xbin in range(h.GetNbinsX()+2):
        x = h.GetXaxis().GetBinCenter(xbin)
        for ybin in range(h.GetNbinsY()+2):
            y = h.GetYaxis().GetBinCenter(ybin)
            newBin = hnew.FindFixBin(x,y)
            hnew.SetBinContent(newBin, hnew.GetBinContent(newBin)+h.GetBinContent(xbin,ybin))
            hnew.SetBinError( newBin, sqrt(hnew.GetBinError(newBin)**2 + h.GetBinError(xbin,ybin)**2 ) )
    return hnew


def rebin( h, binEdges, scale=True ):
    checkRebinningConsistence( h.GetXaxis(), binEdges )
    import array
    binEdgesArr = array.array( 'd', binEdges )
    hnew = h.Rebin( len(binEdges)-1, "new", binEdgesArr )
    hnew.drawOption_ = h.drawOption_ if hasattr( h, "drawOption_" ) else ""
    if scale: hnew.Scale( 1., "width" )
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
    h.SetBinContent( source, 0 )
    h.SetBinError( source, 0 )

def appendFlowBin( h, under=True, over=True ):
    if under:
        mergeBins( h, 1, 0 )
    if over:
        mergeBins( h, h.GetNbinsX(), h.GetNbinsX()+1 )

def integralAndError( h, binx1=0, binx2=-1, bins=True ):
    if not bins:
        binx1 = h.FindFixBin(binx1)
        binx2 = h.FindFixBin(binx2)
    e = ROOT.Double()
    c = h.IntegralAndError(binx1,binx2,e)
    return c,e


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


def automaticRebinner( hlist, minEvents=3 ):
    # Ereates an array of bin edges on a list of histograms, such that each histogram
    # has at least 'maxEvents' events in each bin.
    if minEvents == 0: minEvents = 1e-10
    out = []
    tmp = [0]*len(hlist)
    nBins = hlist[0].GetNbinsX()

    # check overflow bin
    overflowList = [ h.GetBinContent( nBins+1 ) for h in hlist ]
    if min(overflowList) > minEvents:
        out.append( hlist[0].GetBinLowEdge(nBins+1) )

    for bin in range(nBins, -1, -1 ):

        contents = [ h.GetBinContent( bin ) for h in hlist ]
        tmp = [sum(x) for x in zip(tmp, contents)]

        # check upper bondary of last entry
        if not out:
            if contents != [0]*len(hlist):
                out.append( hlist[0].GetBinLowEdge(bin+1) )

        else:
            if min(contents) > minEvents:
                out.append( hlist[0].GetBinLowEdge(bin) )
                tmp = [0]*len(hlist)

    print out[::-1]

def getMinimum( hists ):
    # do not use TH1.GetMinimum(), since it returns the minimum set by SetMinimum()
    return min( [ h.GetBinContent(h.GetMinimumBin()) for h in hists ] )

def setMinMaxForLog():
    allStuff = [ i for i in ROOT.gPad.GetCanvas().GetListOfPrimitives()]
    allH = []
    for h in allStuff:
        if isinstance(h, ROOT.THStack):
            for sh in h.GetHists():
                allH.append(sh)
        elif isinstance(h, ROOT.TH1):
            allH.append(h)
    minC = getMinimum( allH )
    minExp = 1/maxBinWidth(allH[0])
    maxC = max( [ h.GetMaximum() for h in allH ] )
    for i in allH:
        i.SetMaximum( 2.5*maxC )
        i.SetMinimum( .5*max([minC,minExp]) )
    ROOT.gPad.Update()


def save( name, folder="plots/", endings=[".pdf"], normal=True, log=True ):
    if normal:
        for ending in endings:
            ROOT.gPad.GetCanvas().SaveAs( folder+name+ending )
    if log:
        setMinMaxForLog()
        ROOT.gPad.GetCanvas().SetLogy()
        for ending in endings:
            ROOT.gPad.GetCanvas().SaveAs( folder + name + "_log" + ending )


def getBinnigsFromName( name ):
    out = { "": None }
    # get histogram name
    if binCfg.has_section(name):
        for binningName, binning in binCfg.items( name ):
            binning = [ float(x) for x in binning.split(" ") ]
            out[binningName] = binning
    return out

def drange(start, stop, n):
    out = [start]
    step = 1.*(stop-start)/n
    while out[-1] < stop:
        out.append( out[-1] + step )
    return out

def interpolate2D( h ):
    for xbin in range(1,h.GetNbinsX()+1):
        x = h.GetXaxis().GetBinCenter(xbin)
        for ybin in range(1,h.GetNbinsY()+1):
            y = h.GetYaxis().GetBinCenter(ybin)
            c = h.GetBinContent(xbin,ybin)
            if c: continue
            ctop = h.GetBinContent(xbin,ybin+1)
            cbot = h.GetBinContent(xbin,ybin-1)
            crig = h.GetBinContent(xbin+1,ybin)
            clef = h.GetBinContent(xbin-1,ybin)
            intPoints = []
            if cbot and ctop: intPoints.extend( [cbot, ctop])
            if crig and clef: intPoints.extend( [crig, clef])
            newC = sum(intPoints)/len(intPoints) if len(intPoints) else 0
            h.SetBinContent(xbin,ybin, newC )
    return h

def diagonalFlip( original ):
    # original, flipped are both TH2
    flipped = original.Clone(original.GetName()+"flipped")
    flipped.SetTitle("{};{};{}".format(
            original.GetTitle(),
            original.GetYaxis().GetTitle(),
            original.GetXaxis().GetTitle()
        ))

    for xbin in range(original.GetNbinsX()+2):
        for ybin in range(original.GetNbinsY()+2):
            flipped.SetBinContent( ybin, xbin, original.GetBinContent(xbin,ybin) )
            flipped.SetBinError( ybin, xbin, original.GetBinError(xbin,ybin) )
    return flipped

def maxBinWidth( h ):
    return max([ h.GetBinWidth(bin) for bin in range(h.GetNbinsX()+2) ])

def metricPrefix( n ):
    for unit in ["","K","M","G","T"]:
        if abs(n) < 1000:
            return "{:3.1f}{}".format(n, unit)
        n /= 1000
    return "{:3.1f}".format(n)

def loopH2( h2 ):
    return [(xbin,ybin) for xbin in range(h2.GetNbinsX()+2) for ybin in range(h2.GetNbinsY()+2)]

def printH2(h2,flow=True):
    for xbin, ybin in loopH2(h2):
        if not flow and ( not xbin or not xbin or xbin==h2.GetNbinsX()+1 or ybin==h2.GetNbinsY()+1 ): continue
        print xbin, ybin, h2.GetBinContent(xbin,ybin), "±", h2.GetBinError(xbin,ybin)

def stdHist(dataset, name, binning=None, xCut=True, cut1=0, cut2=1e8):
    h = dataset.getHist(name)
    if isinstance(h, ROOT.TH2):
        if xCut: h = h.ProjectionY(randomName(), h.GetXaxis().FindFixBin(cut1), h.GetYaxis().FindFixBin(cut2))
        else:    h = h.ProjectionX(randomName(), h.GetYaxis().FindFixBin(cut1), h.GetXaxis().FindFixBin(cut2))
    if binning: h = rebin(h, binning)
    appendFlowBin(h)
    h.SetYTitle(getYAxisTitle(h))
    return h

def drawOpt(h, style):
    if style == "data":
        h.SetLineColor(1)
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.8)
        h.drawOption_="pe"
    if style == "sys":
        c = h.GetLineColor()
        h.SetFillColor(c)
        h.SetMarkerColor(c)
        h.SetFillStyle(3333)
        h.drawOption_ = "e2"


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

    def __init__( self, drawAll=True, sim=False, status="Private Work", info="" ):
        # todo: include margins, etc
        if sim:
            self.cms = ROOT.TLatex( 0.2, .895, "#font[61]{CMS} #scale[0.76]{#font[52]{Simulation}}" )
        else:
            self.cms = ROOT.TLatex( 0.2, .895, "#font[61]{CMS}" )
        self.pub = ROOT.TLatex( 0.2, .865, "#scale[0.76]{#font[52]{%s}}"%status )
        self.lum = ROOT.TLatex( .63, .95, "%.2f fb^{-1} (%s TeV)"%(self.intLumi/1000., self.cmsEnergy) )
        if info: self.info = ROOT.TLatex( .2, .95, info )

        if drawAll:
            self.draw()
