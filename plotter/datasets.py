import ROOT
import copy

path = "/home/knut/"
path = "../histogramProducer/"

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
        out.label += " + "+dset.label

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

# multijet

gjets400 = Dataset( "GJets_HT-400to600", 62.05, 4802346, ROOT.kCyan+2 ) #eff_match = 0.075
gjets600 = Dataset( "GJets_HT-600toInf", 20.87, 4419582, ROOT.kCyan ) # eff_match = 0.063
gjets = gjets600 + gjets400
gjets.label = "#gamma+Jet"

qcd1000 = Dataset( "QCD_HT_1000ToInf", 769.7, 1130720+333733, ROOT.kBlue ) #eff_match = 0.2358
qcd500 = Dataset( "QCD_HT-500To1000", 26740, 3214312+849033, ROOT.kBlue+2 ) #eff_match = 0.2103
qcd = qcd1000 + qcd500
qcd.label = "Multijet"

# e->g

wjets400 = Dataset( "WJetsToLNu_HT-400to600", 55.61, 4640594, ROOT.kRed+2 ) # 0.075
wjets600 = Dataset( "WJetsToLNu_HT-600toInf", 18.81, 4581841, ROOT.kRed ) # 0.063
wjets = wjets600 + wjets400
wjets.label = "W#rightarrowl#nu"

ttjets = Dataset( "TTJets", 424.5, 25446993, ROOT.kRed-4 ) # 0.208
ttjets.label = "t#bar{t}"

# isr

znn400 = Dataset( "ZJetsToNuNu_HT-400to600", 11.99, 4659330, ROOT.kMagenta+2 ) #eff_match = 0.069
znn600 = Dataset( "ZJetsToNuNu_HT-600toInf", 4.113, 4547145, ROOT.kMagenta ) #eff_match = 0.062
znn = znn600 + znn400
znn.label = "Z#rightarrow#nu#nu"



