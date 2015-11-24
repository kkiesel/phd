import ROOT
import copy
import auxiliary as aux
import os.path

from main import intLumi

path = "../histogramProducer/"

class Dataset:
    names = []
    files = []
    xsecs = []
    ngens = []
    lname = []

    color = None
    label = ""

    def __radd__( self, dset ):
        if not dset: return self

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
            "\nfiles: "+", ".join( self.files ) + \
            "\nxsecs: "+", ".join( str(i) for i in self.xsecs ) + \
            "\nngens: "+", ".join( str(i) for i in self.ngens )

    def getNgenFromFile( self ):
        for ngen, file in zip(self.ngens, self.files):
            print file, ngen, aux.getNgen( file )


    def getHist( self, name ):
        h0 = None
        for i in range( len(self.files) ):
            h = aux.getFromFile( self.files[i], name )
            if isinstance( h, ROOT.TH1 ):
                if self.xsecs[i]:
                    h.Scale( intLumi * self.xsecs[i] / self.ngens[i] )
                h.SetLineColor( self.color )
                h.SetMarkerColor( self.color )
            if h0: h0.Add( h )
            else: h0 = h
        return h0


# data
data = Dataset( "SinglePhoton", 0, ROOT.kBlack )

data_2015C = Dataset( "SinglePhoton_Run2015C_25ns-05Oct2015-v1", 0, ROOT.kBlack )
data_2015D = Dataset( "SinglePhoton_Run2015D-05Oct2015-v1", 0, ROOT.kBlack )
data_2015D.label = "05Oct2015 reco"
data_prompt = Dataset( "SinglePhoton_Run2015D-PromptReco-v4", 0, ROOT.kBlack )
data_prompt.label = "PromptReco-v4"

data = data_prompt+data_2015D # +data_2015C
data.label = "Data"

dataHt_2015C = Dataset( "JetHT_Run2015C_25ns-05Oct2015-v1", 0, ROOT.kBlack )
dataHt_2015D = Dataset( "JetHT_Run2015D-05Oct2015-v1", 0, ROOT.kBlack )
dataHt_2015D.label = "05Oct2015 reco"
dataHt_prompt = Dataset( "JetHT_Run2015D-PromptReco-v4", 0, ROOT.kBlack )
dataHt_prompt.label = "Prompt"
dataHt = dataHt_prompt+dataHt_2015D #+dataHt_2015C
dataHt.label = "Data (JetHt)"

# k-factors from twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns

# multijet

gjets_pt15 = Dataset( "GJet_Pt-15ToInf", 364375, ROOT.kCyan-2, "GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8" )

gjets40 = Dataset( "GJets_HT-40To100", 20730, ROOT.kCyan-1, "GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets100 = Dataset( "GJets_HT-100To200", 9226, ROOT.kCyan+4, "GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets200 = Dataset( "GJets_HT-200To400", 2300, ROOT.kCyan+3, "GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets400 = Dataset( "GJets_HT-400To600", 277.4, ROOT.kCyan+2, "GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets600 = Dataset( "GJets_HT-600ToInf", 93.38, ROOT.kCyan, "GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

gjets = gjets40 + gjets100 + gjets200 + gjets400 + gjets600
gjets.label = "#gamma+Jet"

qcd100 = Dataset( "QCD_HT100to200", 27850000, ROOT.kBlue+1, "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd200 = Dataset( "QCD_HT200to300", 1717000, ROOT.kBlue+2, "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd300 = Dataset( "QCD_HT300to500", 361300, ROOT.kBlue+3, "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd500 = Dataset( "QCD_HT500to700", 31630, ROOT.kBlue+4, "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd700 = Dataset( "QCD_HT700to1000", 6802, ROOT.kBlue+3, "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1000 = Dataset( "QCD_HT1000to1500", 1206, ROOT.kBlue+2, "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1500 = Dataset( "QCD_HT1500to2000", 120.4, ROOT.kBlue+1, "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd2000 = Dataset( "QCD_HT2000toInf", 25.24,  ROOT.kBlue, "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

qcd = qcd100 + qcd200 + qcd300 + qcd500 + qcd700 + qcd1000 + qcd1500 + qcd2000
qcd.label = "Multijet"

emqcd15 = Dataset( "QCD_Pt-15to20_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )
emqcd20 = Dataset( "QCD_Pt-20to30_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )
emqcd30 = Dataset( "QCD_Pt-30o50_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )
emqcd50 = Dataset( "QCD_Pt-50to80_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )
emqcd80 = Dataset( "QCD_Pt-80to120_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-80o120_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )
emqcd120 = Dataset( "QCD_Pt-120to170_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-120o170_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )
emqcd170 = Dataset( "QCD_Pt-170to300_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-170o300_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )
emqcd300 = Dataset( "QCD_Pt-300toInf_EMEnriched", 0, ROOT.kBlue+0, "QCD_Pt-300oInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8" )

emqcd = emqcd15 + emqcd20 + emqcd30 + emqcd50 + emqcd80 + emqcd120 + emqcd170 + emqcd300
emqcd.label = "QCD EM enriched"

# electroweak

ttjets = Dataset( "TTJets", 670.3,  ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
ttjets.label = "t#bar{t}"

wjets100 = Dataset( "WJetsToLNu_HT-100To200", 1345., ROOT.kRed-6, "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets200 = Dataset( "WJetsToLNu_HT-200To400", 359.7, ROOT.kRed-5, "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets400 = Dataset( "WJetsToLNu_HT-400To600", 48.91, ROOT.kRed-4, "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets600 = Dataset( "WJetsToLNu_HT-600To800", 12.05, ROOT.kRed-3, "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets800 = Dataset( "WJetsToLNu_HT-800To1200", 5.501, ROOT.kRed-2, "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets1200 = Dataset( "WJetsToLNu_HT-1200To2500", 1.329, ROOT.kRed-1, "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets2500 = Dataset( "WJetsToLNu_HT-2500ToInf", 0.03216, ROOT.kRed, "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets600_inf = Dataset( "WJetsToLNu_HT-600ToInf", 18.77, ROOT.kRed, "WJetsToLNu_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )

# k-factor
wjetsSamples = [ wjets100, wjets200, wjets400, wjets600, wjets800, wjets1200, wjets2500 ]
for ds in wjetsSamples+[wjets600_inf]:
    ds.xsecs = [ ds.xsecs[0] * 1.21 ]

wjets = wjets600_inf + wjets400 + wjets200 + wjets100 # no high HT samples
wjets = sum( wjetsSamples )
wjets.label = "W#rightarrowl#nu"

# isr
ttg = Dataset( "TTGJets", 3.697, ROOT.kOrange, "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia88" )
ttg.label = "#gammat#bar{t}"

wg_mc = Dataset( "WGToLNuG-amcatnloFXFX", 489., ROOT.kRed-2, "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
wg_mg = Dataset( "WGToLNuG-madgraphMLM", 405.271, ROOT.kRed-3, "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wg_pt500 = Dataset( "WGToLNuG_PtG-500", 0.0117887, ROOT.kRed-1, "WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )

# znunu
znunu100 = Dataset( "ZJetsToNuNu_HT-100To200", 280.47, ROOT.kMagenta+3 , "ZJetsToNuNu_HT-100To200_13TeV-madgraph" )
znunu200 = Dataset( "ZJetsToNuNu_HT-200To400", 78.36, ROOT.kMagenta+2 , "ZJetsToNuNu_HT-200To400_13TeV-madgraph" )
znunu400 = Dataset( "ZJetsToNuNu_HT-400To600", 10.94, ROOT.kMagenta+1 , "ZJetsToNuNu_HT-400To600_13TeV-madgraph" )
znunu600 = Dataset( "ZJetsToNuNu_HT-600ToInf", 4.20, ROOT.kMagenta+0 , "ZJetsToNuNu_HT-600ToInf_13TeV-madgraph" )

zg_130 = Dataset( "ZNuNuGJets_MonoPhoton_PtG-130", 0.223, ROOT.kMagenta+1, "ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph" )

znunuSamples = znunu100, znunu200, znunu400, znunu600
for ds in znunuSamples:
    # apply k-factor
    ds.xsecs = [ ds.xsecs[0] * 1.23 ]

znunu = sum( znunuSamples )
znunu.label = "Z#rightarrow#nu#nu"


# cross section from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
dy = Dataset( "DYJetsToLL_M-50", 6025.2, ROOT.kRed+3, "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
dy.label = "Z#rightarrowll"


# signal samples

import collections
class SampleCollection(collections.MutableMapping):
    """Dictionary used to store signal datasets.
    The first time a dataset is requested, it is created"""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        if key not in self.store:
            self.store[key] = Dataset( aux.getDatasetFromKey(key), col=ROOT.kGreen+len(self) )
            self.store[key].label = aux.getSignalLabel( key )
        return self.store[key]

    def __setitem__(self, key, value):
        self.store[key] = value

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

signal = SampleCollection()

