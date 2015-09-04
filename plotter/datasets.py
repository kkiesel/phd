import ROOT
import copy
import auxiliary as aux

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

    def __init__( self, n, xsec=-1, ngen=-1, col=ROOT.kBlack ):
        if xsec == -1: xsec = aux.getXsecFromName( n )
        if ngen == -1: ngen = aux.getNgen( path + n + "_hists.root" )
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

# /QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
qcd1000 = Dataset( "QCD_HT1000to1500", 1064, 4695055, ROOT.kBlue+2 )

# /QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8&shown=268307
qcd1500 = Dataset( "QCD_HT1500to2000", 121.5, 3848411, ROOT.kBlue+1 )

# /QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
qcd2000 = Dataset( "QCD_HT2000toInf", 25.42,  1961774, ROOT.kBlue )

qcd = qcd2000 + qcd1500 + qcd1000 + qcd700 + qcd500 + qcd300
qcd.label = "Multijet"
# todo: check nEvents

# electroweak

# /TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8
ttjets = Dataset( "TTJets", 670.3, 43024978, ROOT.kRed+2 )
ttjets.label = "t#bar{t}"

# /WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
wjets400 = Dataset( "WJetsToLNu_HT-400To600", 48.98, 1901705, ROOT.kRed )

# /WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8
wjets600 = Dataset( "WJetsToLNu_HT-600ToInf", 18.77, 1036108, ROOT.kRed-4 )

wjets = wjets600 + wjets400
wjets.label = "W#rightarrowl#nu"

# isr

# /TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM
# https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name=TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8&shown=262168
ttg = Dataset( "TTGJets", 3.697, 4784923, ROOT.kOrange )
ttg.label = "#gammat#bar{t}"

# signal samples

t5gg = Dataset( "T5gg_1500_1000", col=ROOT.kMagenta )
t5gg.label = "T5gg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=1000"
t5hg = Dataset( "T5hg_1500_1000", col=ROOT.kMagenta+4 )
t5hg.label = "T5hg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=1000"
t2ttgg = Dataset( "T2ttgg_850_650_fast", col=ROOT.kMagenta+2 )
t2ttgg.label = "T2ttgg m(#tilde{t})=850 m(#tilde{#chi}^{0}_{1})=650"


