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

        out = copy.deepcopy(self)

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

    def __str__( self ):
        return "Dataset: " + self.label + "\ncolor: " + str(self.color) + \
            "\nnames: "+", ".join( self.names ) + \
            "\nfiles: "+", ".join( str(i) for i in self.xsecs ) + \
            "\nngens: "+", ".join( str(i) for i in self.ngens )


""" Phys14 samples
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
"""

# multijet

# /GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
gjets400 = Dataset( "GJets_HT-400To600", 273, 2476770, ROOT.kCyan+2 )

# /GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
gjets600 = Dataset( "GJets_HT-600ToInf", 94.5, 2550765, ROOT.kCyan )

gjets = gjets600 + gjets400
gjets.label = "#gamma+Jet"

# /QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
qcd300 = Dataset( "QCD_HT300to500", 366800, 20086103, ROOT.kBlue+5 )

# /QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
qcd500 = Dataset( "QCD_HT500to700", 29370, 19542847, ROOT.kBlue+4 )

# /QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
qcd700 = Dataset( "QCD_HT700to1000", 6524, 15011016, ROOT.kBlue+3 )

# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
qcd1000 = Dataset( "QCD_HT1000to1500", 1064, 4695055, ROOT.kBlue+2 )

# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8&shown=268307
qcd1500 = Dataset( "QCD_HT1500to2000", 121.5, 3848411, ROOT.kBlue+1 )

# /QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
qcd2000 = Dataset( "QCD_HT2000toInf", 25.42,  1961774, ROOT.kBlue )

qcd = qcd2000 + qcd1500 + qcd1000 + qcd700 + qcd500 + qcd300
qcd = qcd2000 + qcd700 + qcd500 + qcd300
qcd.label = "incomplete Multijet"
# todo: check nEvents

# electroweak

# /TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
ttjets = Dataset( "TTJets", 670.3, 43024978, ROOT.kRed+2 )
ttjets.label = "t#bar{t}"

# /WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
wjets400 = Dataset( "WJetsToLNu_HT-400To600", 47.9, 1901705, ROOT.kRed )

# /WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
wjets600 = Dataset( "WJetsToLNu_HT-600ToInf", 19.9, 1036108, ROOT.kRed-4 )

wjets = wjets600 + wjets400
wjets.label = "W#rightarrowl#nu"


