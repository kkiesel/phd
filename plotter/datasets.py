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
data.label = "Data(Run D)"


dataHt = Dataset( "JetHT", 0, ROOT.kBlack )
dataHt.label = "Data"

# k-factors from twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns

# multijet

gjets_pt15 = Dataset( "GJet_Pt-15ToInf", 364375, ROOT.kCyan-2, "GJet_Pt-15ToInf_TuneCUETP8M1_13TeV-pythia8" )

gjets40 = Dataset( "GJets_HT-40To100", 20730, ROOT.kCyan-1, "GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets100 = Dataset( "GJets_HT-100To200", 9226, ROOT.kCyan+4, "GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets200 = Dataset( "GJets_HT-200To400", 2300, ROOT.kCyan+3, "GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets400 = Dataset( "GJets_HT-400To600", 277.4, ROOT.kCyan+2, "GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets600 = Dataset( "GJets_HT-600ToInf", 93.38, ROOT.kCyan, "GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

gjets = gjets600 + gjets400 + gjets200 + gjets100 + gjets40
gjets.label = "#gamma+Jet"

qcd100 = Dataset( "QCD_HT100to200", 27850000, ROOT.kBlue+1, "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd200 = Dataset( "QCD_HT200to300", 1717000, ROOT.kBlue+2, "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd300 = Dataset( "QCD_HT300to500", 361300, ROOT.kBlue+3, "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd500 = Dataset( "QCD_HT500to700", 31630, ROOT.kBlue+4, "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd700 = Dataset( "QCD_HT700to1000", 6802, ROOT.kBlue+3, "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1000 = Dataset( "QCD_HT1000to1500", 1206, ROOT.kBlue+2, "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1500 = Dataset( "QCD_HT1500to2000", 120.4, ROOT.kBlue+1, "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd2000 = Dataset( "QCD_HT2000toInf", 25.24,  ROOT.kBlue, "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

qcd = qcd2000 + qcd1500 + qcd1000 + qcd700 + qcd500 + qcd300 + qcd200 + qcd100
qcd.label = "Multijet"

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
wjetsSamples = [ wjets2500, wjets1200, wjets800, wjets600, wjets400, wjets200, wjets100 ]
for ds in wjetsSamples+[wjets600_inf]:
    ds.xsecs = [ ds.xsecs[0] * 1.21 ]

wjets = wjets600_inf + wjets400 + wjets200 + wjets100 # no high HT samples
wjets = wjets2500 + wjets1200 + wjets800 + wjets600 + wjets400 + wjets200 + wjets100
wjets.label = "W#rightarrowl#nu"

# isr
ttg = Dataset( "TTGJets", 3.697, ROOT.kOrange, "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia88" )
ttg.label = "#gammat#bar{t}"

wg = Dataset( "WGToLNuG", 405.271, ROOT.kRed-3, "" )
wg.label = "#gammaW#rightarrow#gammal#nu"

# cross section from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
dy = Dataset( "DYJetsToLL_M-50", 6025.2, ROOT.kRed+3, "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
dy.label = "Z#rightarrowll"

# signal samples

t5wg_1500_125 = Dataset( "T5wg_1500_125", col=ROOT.kGreen+4 )
t5wg_1500_125.label = "T5wg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=125"

t5wg_1500_1475 = Dataset( "T5wg_1500_1475", col=ROOT.kGreen+3 )
t5wg_1500_1475.label = "T5wg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=1475"


#t5gg = Dataset( "T5gg_1500_1000", col=ROOT.kMagenta )
#t5gg.label = "T5gg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=1000"
#t5hg = Dataset( "T5hg_1500_1000", col=ROOT.kMagenta+4 )
#t5hg.label = "T5hg m(#tilde{g})=1500 m(#tilde{#chi}^{0}_{1})=1000"
#t2ttgg = Dataset( "T2ttgg_850_650_fast", col=ROOT.kGreen+4 )
#t2ttgg.label = "T2ttgg #scale[0.7]{#tilde{t}:850 #tilde{#chi}^{0}_{1}:650}"
#t2ttgg.label = "T2ttgg (x10)"
#t2ttgg.xsecs = [ 10.*t2ttgg.xsecs[0] ]

