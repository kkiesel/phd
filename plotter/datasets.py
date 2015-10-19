import ROOT
import copy
import auxiliary as aux
import os.path

path = "../histogramProducer/"

class Dataset:
    names = []
    files = []
    xsecs = []
    ngens = []
    lname = []

    color = None
    label = ""

    def __add__( self, dset ):

        out = copy.deepcopy(self)

        out.names.extend( dset.names )
        out.files.extend( dset.files )
        out.xsecs.extend( dset.xsecs )
        out.ngens.extend( dset.ngens )
        out.label += " + "+dset.label
        out.lname.extend( dset.lname )

        if not self.color:
            out.color = dset.color
        return out

    def __init__( self, n, xsec=-1, col=ROOT.kBlack, fullname="", ngen=-1 ):
        fname = path+n+"_hists.root"
        if xsec == -1: xsec = aux.getXsecFromName( n )
        if ngen == -1 and os.path.isfile(fname): ngen = aux.getNgen( fname )
        self.names = [ n ]
        self.files = [ fname ]
        self.xsecs = [ xsec ]
        self.ngens = [ ngen ]
        self.color = col
        self.label = n
        self.lname = [ fullname ]

    def mcm( self ):
        for fullname in self.lname:
            print "https://cms-pdmv.cern.ch/mcm/requests?page=0&dataset_name="+fullname


    def __str__( self ):
        return "Dataset: " + self.label + "\ncolor: " + str(self.color) + \
            "\nnames: "+", ".join( self.names ) + \
            "\nfiles: "+", ".join( str(i) for i in self.xsecs ) + \
            "\nngens: "+", ".join( str(i) for i in self.ngens )

    def getNgenFromFile( self ):
        for ngen, file in zip(self.ngens, self.files):
            print file, ngen, aux.getNgen( file )




# data
data = Dataset( "SinglePhoton", 0, ROOT.kBlack )
data.label = "Data(Run D)"


jetHtC = Dataset( "JetHT_Run2015C-PromptReco-v1", 0, ROOT.kBlack )
jetHtC.label = "JetHtC"
jetHtD = Dataset( "JetHT_Run2015D-PromptReco-v3", 0, ROOT.kBlack )
jetHtD.label = "JetHtD"
dataHt = jetHtC + jetHtD
dataHt.label = "Data"

# k-factors from twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns

# multijet

gjets_pt15 = Dataset( "GJet_Pt-15ToInf", 364375, ROOT.kCyan, "GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8" )

gjets40 = Dataset( "GJets_HT-40To100", 20730, ROOT.kCyan-1, "GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets100 = Dataset( "GJets_HT-100To200", 9226, ROOT.kCyan+4, "GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets200 = Dataset( "GJets_HT-200To400", 2281, ROOT.kCyan+3, "GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets400 = Dataset( "GJets_HT-400To600", 273, ROOT.kCyan+2, "GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets600 = Dataset( "GJets_HT-600ToInf", 94.5, ROOT.kCyan, "GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

gjets = gjets600 + gjets400 + gjets200 + gjets100 + gjets40
gjets.label = "#gamma+Jet"

qcd100 = Dataset( "QCD_HT100to200", 27540000, ROOT.kBlue+7, "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd200 = Dataset( "QCD_HT200to300", 1735000, ROOT.kBlue+6, "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd300 = Dataset( "QCD_HT300to500", 366800, ROOT.kBlue+5, "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd500 = Dataset( "QCD_HT500to700", 29370, ROOT.kBlue+4, "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd700 = Dataset( "QCD_HT700to1000", 6524, ROOT.kBlue+3, "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1000 = Dataset( "QCD_HT1000to1500", 1064, ROOT.kBlue+2, "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1500 = Dataset( "QCD_HT1500to2000", 121.5, ROOT.kBlue+1, "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd2000 = Dataset( "QCD_HT2000toInf", 25.42,  ROOT.kBlue, "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

qcd = qcd2000 + qcd1500 + qcd1000 + qcd700 + qcd500 + qcd300 + qcd200 + qcd100
qcd.label = "Multijet"

# electroweak

ttjets = Dataset( "TTJets", 670.3,  ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
ttjets.label = "t#bar{t}"

wjets400 = Dataset( "WJetsToLNu_HT-400To600", 48.98*1.21, ROOT.kRed, "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets600 = Dataset( "WJetsToLNu_HT-600ToInf", 18.77*1.21, ROOT.kRed-4, "WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )

wjets = wjets600 + wjets400
wjets.label = "W#rightarrowl#nu"

# isr

ttg = Dataset( "TTGJets", 3.697, ROOT.kOrange, "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia88" )
ttg.label = "#gammat#bar{t}"

wg = Dataset( "WGToLNuG", 405.271, ROOT.kRed-3 )
wg.label = "#gammaW#rightarrow#gammal#nu"

dy = Dataset( "DYJetsToLL_M-50", 6104, ROOT.kRed+3 )
dy.label = "Z#rightarrowll"

# signal samples

#t5gg = Dataset( "T5gg_1500_1000", col=ROOT.kMagenta )
#t5gg.label = "T5gg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=1000"
#t5hg = Dataset( "T5hg_1500_1000", col=ROOT.kMagenta+4 )
#t5hg.label = "T5hg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=1000"
#t2ttgg = Dataset( "T2ttgg_850_650_fast", col=ROOT.kGreen+4 )
#t2ttgg.label = "T2ttgg #scale[0.7]{#tilde{t}:850 #tilde{#chi}^{0}_{1}:650}"
#t2ttgg.label = "T2ttgg (x10)"
#t2ttgg.xsecs = [ 10.*t2ttgg.xsecs[0] ]

