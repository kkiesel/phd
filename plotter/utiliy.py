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

class Label():
    # Create labels
    # Usage:
    # * With Labels(), all default labels will be printed
    # * With Labels(False), the method is only initiated and labels can be modified before calling the 'draw' method

    cmsEnergy = 13 #TeV
    lumi = 19.7 # fb^{-1}

    cms = ROOT.TLatex( ROOT.gPad.GetLeftMargin()+0.02, .895, "#font[61]{CMS}")
    sim = ROOT.TLatex( ROOT.gPad.GetLeftMargin()+0.02, .89, "#scale[0.76]{#font[52]{Simulation}")
    pub = ROOT.TLatex( ROOT.gPad.GetLeftMargin()+0.02, .885, "#scale[0.76]{#font[52]{Private Work}")
    lum = ROOT.TLatex( .68, .895, "%s fb^{-1} (%s TeV)"%(lumi, cmsEnergy) )

    def draw( self ):
        varDict = vars( self )

        for varName, obj in varDict.iteritems():
            if isInstance( obj, ROOT.TLatex ):
                obj.DrawNDC()

    def __init__( self, drawAll=True ):
        if drawAll:
            self.draw()
