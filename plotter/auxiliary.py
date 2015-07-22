import ROOT

def getFromFile( filename, objectname ):
    # todo: check if file exists
    f = ROOT.TFile( filename )

    obj = f.Get( objectname )
    if not obj:
        print "ERROR, {} not in {}".format( objectname, filename )
        return
    obj = ROOT.gROOT.CloneObject( obj )

    return obj

def randomName():
    # Returns a random alphanumeric string
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

    def __init__( self, drawAll=True, sim=True, status="Private Work" ):
        # todo: include margins, etc
        if sim:
            self.cms = ROOT.TLatex( 0.2, .895, "#font[61]{CMS} #scale[0.76]{#font[52]{Simulation}}")
        else:
            self.cms = ROOT.TLatex( 0.2, .895, "#font[61]{CMS}")
        self.pub = ROOT.TLatex( 0.2, .865, "#scale[0.76]{#font[52]{%s}}"%status)
        self.lum = ROOT.TLatex( .63, .95, "%s fb^{-1} (%s TeV)"%(self.intLumi/1000., self.cmsEnergy) )

        if drawAll:
            self.draw()
