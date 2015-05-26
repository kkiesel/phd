import ROOT
import copy

path = "/home/knut/"

class Dataset:
    names = []
    files = []
    xsecs = []
    ngens = []

    color = None
    label = ""

    def __add__( self, dset ):

        out = copy.deepcopy(dset)

        out.names.extend( dset.names )
        out.files.extend( dset.files )
        out.xsecs.extend( dset.xsecs )
        out.ngens.extend( dset.ngens )

        if not self.color:
            out.color = dset.color
        return out

    def __init__( self, n, xsec, ngen, col ):
        self.names = [ n ]
        self.files = [ path + n + "_hists.root" ]
        self.xsecs = [ xsec ]
        self.ngens = [ ngen ]
        self.color = col
        self.label = n


gjets400 = Dataset( "GJets_HT-400to600", 62.05, 4802346, ROOT.kCyan+2 ) #eff_match = 0.075
gjets600 = Dataset( "GJets_HT-600toInf", 20.87, 4419582, ROOT.kCyan ) # eff_match = 0.063
gjets = gjets600 + gjets400
gjets.label = "#gamma + Jet"

znn400 = Dataset( "ZJetsToNuNu_HT-400to600", 11.99, 4659330, ROOT.kMagenta+2 ) #eff_match = 0.069
znn600 = Dataset( "ZJetsToNuNu_HT-600toInf", 4.113, 4547145, ROOT.kMagenta ) #eff_match = 0.062

znn = znn600 + znn400
znn.label = "Z#rightarrow#nu#nu"



