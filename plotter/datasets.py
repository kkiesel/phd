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
                    h.Scale( aux.intLumi * self.xsecs[i] / self.ngens[i] )
                h.SetLineColor( self.color )
                h.SetMarkerColor( self.color )
            if h0: h0.Add( h )
            else: h0 = h
        return h0

    def getLatexTableHeader( self ):
        return "\\begin{tabular}{c|c|c}\nPrimary Dataset & cross section (pb) & effective Luminosity (/fb) \\\\\\hline"

    def getLatexTableLine( self ):
        # full samplename & xsec [pb] & effective Luminosity [/fb]
        out = ""
        for lname, xsec, ngen in zip(self.lname,self.xsecs,self.ngens):
            out += "{} & {} & {:.3g} \\\\\n".format( lname.replace("_","\\_"),xsec,0.001*ngen/xsec )
        return out


# data
data = Dataset("SinglePhoton_Run2016B-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("SinglePhoton_Run2016C-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("SinglePhoton_Run2016D-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("SinglePhoton_Run2016E-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("SinglePhoton_Run2016F-PromptReco-v1", 0, ROOT.kBlack ) \
    + Dataset("SinglePhoton_Run2016G-PromptReco-v1", 0, ROOT.kBlack ) \
    + Dataset("SinglePhoton_Run2016H-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("SinglePhoton_Run2016H-PromptReco-v3", 0, ROOT.kBlack )
data.label = "Data"

dataHt = Dataset("JetHT_Run2016B-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016C-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016D-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016E-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016F-PromptReco-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016G-PromptReco-v1", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016H-PromptReco-v2", 0, ROOT.kBlack ) \
    + Dataset("JetHT_Run2016H-PromptReco-v3", 0, ROOT.kBlack )
dataHt.label = "Data (JetHt)"

# k-factors from twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns

# multijet

gjets40 = Dataset( "GJets_HT-40To100", 20790, ROOT.kCyan-1, "GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets100 = Dataset( "GJets_HT-100To200", 9238, ROOT.kCyan+4, "GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
gjets200 = Dataset( "GJets_HT-200To400", 2305, ROOT.kCyan+3, "GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets400 = Dataset( "GJets_HT-400To600", 274.4, ROOT.kCyan+2, "GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
gjets600 = Dataset( "GJets_HT-600ToInf", 93.46, ROOT.kCyan, "GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

gjets = gjets40 + gjets100 + gjets200 + gjets400 + gjets600
gjets = gjets200 + gjets400 + gjets600
gjets.label = "#gamma+Jet"

qcd100 = Dataset( "QCD_HT100to200", 27990000, ROOT.kBlue+1, "QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
qcd200 = Dataset( "QCD_HT200to300", 1712000, ROOT.kBlue+2, "QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd300 = Dataset( "QCD_HT300to500", 347700, ROOT.kBlue+3, "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd500 = Dataset( "QCD_HT500to700", 32100, ROOT.kBlue+4, "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd700 = Dataset( "QCD_HT700to1000", 6831, ROOT.kBlue+3, "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1000 = Dataset( "QCD_HT1000to1500", 1207, ROOT.kBlue+2, "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd1500 = Dataset( "QCD_HT1500to2000", 119.9, ROOT.kBlue+1, "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )
qcd2000 = Dataset( "QCD_HT2000toInf", 25.24,  ROOT.kBlue, "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"  )

qcd = qcd100 + qcd200 + qcd300 + qcd500 + qcd700 + qcd1000 + qcd1500 + qcd2000
qcd = qcd500 + qcd700 + qcd1000 + qcd1500 + qcd2000
qcd.label = "Multijet"

# electroweak
# https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns#TTbar NNLO +19.77 -29.20 +35.06 -35.06
ttjets_nlo = Dataset( "TTJets-amcatnloFXFX", 6.675e+02,  ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
ttjets_nlo.xsecs[0] = 831.76 #+19.77 -29.20 +35.06 -35.06
ttjets = Dataset( "TTJets-madgraphMLM", 5.098e+02,  ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
ttjets.xsecs[0] = 831.76 #+19.77 -29.20 +35.06 -35.06
ttjets.label = "t#bar{t}"

ttjets0 = Dataset("TTJets_HT-0to600", 5.098e+02, ROOT.kRed+2, "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets600 = Dataset("TTJets_HT-600to800_ext", 1.61, ROOT.kRed+2, "TTJets_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets800 = Dataset("TTJets_HT-800to1200_ext", 0.663, ROOT.kRed+2, "TTJets_HT-800to1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets1200 = Dataset("TTJets_HT-1200to2500_ext", 0.12, ROOT.kRed+2, "TTJets_HT-1200to2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets2500 = Dataset("TTJets_HT-2500toInf", 0.00143, ROOT.kRed+2, "TTJets_HT-2500toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
ttjets_ht_samples = ttjets0, ttjets600, ttjets800, ttjets1200, ttjets2500
for ds in ttjets_ht_samples: ds.xsecs = [ ds.xsecs[0] * 831.76/5.098e+02 ] # k-factor
ttjets_ht = sum(ttjets_ht_samples)
ttjets_ht.label = "t#bar{t}"

wjets100 = Dataset( "WJetsToLNu_HT-100To200", 1345., ROOT.kRed-6, "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets200 = Dataset( "WJetsToLNu_HT-200To400", 359.7, ROOT.kRed-5, "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets400 = Dataset( "WJetsToLNu_HT-400To600", 48.91, ROOT.kRed-4, "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets600 = Dataset( "WJetsToLNu_HT-600To800", 12.05, ROOT.kRed-3, "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets800 = Dataset( "WJetsToLNu_HT-800To1200", 5.501, ROOT.kRed-2, "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets1200 = Dataset( "WJetsToLNu_HT-1200To2500", 1.329, ROOT.kRed-1, "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjets2500 = Dataset( "WJetsToLNu_HT-2500ToInf", 0.03216, ROOT.kRed, "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wjetsSamples = [ wjets200, wjets400, wjets600, wjets800, wjets1200, wjets2500 ]
for ds in wjetsSamples: ds.xsecs = [ ds.xsecs[0] * 1.21 ] # k-factor
wjets = sum( wjetsSamples )
wjets.label = "W#rightarrowl#nu"

ttg = Dataset( "TTGJets", 3.697, ROOT.kOrange, "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8" )
ttg.label = "#gammat#bar{t}"

tg = Dataset( "TGJets_amcatnlo_madspin", 2.967, ROOT.kOrange+2, "TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8" )
tg.label = "#gammat"

#wg_pt500 = Dataset( "WGToLNuG_PtG-500", 0.0117887, ROOT.kRed-1, "WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wg_mc = Dataset( "WGToLNuG-amcatnloFXFX_ext", 512.1, ROOT.kRed-2, "WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
wg_mc.label = "#gammaW#rightarrow#gammal#nu (mc@NLO)"
wg_mg = Dataset( "WGToLNuG-madgraphMLM", 377.2*1.34, ROOT.kRed-3, "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wg_mg.label = "#gammaW#rightarrow#gammal#nu (MG)"
wg_mg_0to130 = Dataset( "WGToLNuG-madgraphMLM_PtG-0to130", 377.2*1.34, ROOT.kRed-3, "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wg_mg_0to130.label = "#gammaW p_{T}<130GeV"
wg_130 = Dataset( "WGJets_MonoPhoton_PtG-130", 0.6565*1.34, ROOT.kRed-1, "WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph" )
wg_130.label = "#gammaW p_{T}>130GeV"

# k-factor from http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2016/078 EXO-16-014
wg = Dataset( "WGToLNuG-madgraphMLM_PtG-0to130", 405.271*1.34, ROOT.kRed-3, "WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" )
wg = wg + wg_130
wg.label = "#gammaW"

znunu100 = Dataset( "ZJetsToNuNu_HT-100To200", 280.47, ROOT.kMagenta-2 , "ZJetsToNuNu_HT-100To200_13TeV-madgraph" )
znunu200 = Dataset( "ZJetsToNuNu_HT-200To400_ext", 77.67, ROOT.kMagenta-1 , "ZJetsToNuNu_HT-200To400_13TeV-madgraph" )
znunu400 = Dataset( "ZJetsToNuNu_HT-400To600", 10.73, ROOT.kMagenta+4 , "ZJetsToNuNu_HT-400To600_13TeV-madgraph" )
znunu600 = Dataset( "ZJetsToNuNu_HT-600To800", 2.559, ROOT.kMagenta+3 , "ZJetsToNuNu_HT-600To800_13TeV-madgraph" )
znunu800 = Dataset( "ZJetsToNuNu_HT-800To1200", 1.1796, ROOT.kMagenta+2 , "ZJetsToNuNu_HT-800To1200_13TeV-madgraph" )
znunu1200 = Dataset( "ZJetsToNuNu_HT-1200To2500", 0.28833, ROOT.kMagenta+1 , "ZJetsToNuNu_HT-1200To2500_13TeV-madgraph" )
znunu2500 = Dataset( "ZJetsToNuNu_HT-2500ToInf", 0.006945, ROOT.kMagenta+0 , "ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph" )
znunu600Inf = Dataset( "ZJetsToNuNu_HT-600ToInf", 4.116, ROOT.kMagenta+0 , "ZJetsToNuNu_HT-600ToInf_13TeV-madgraph" )
znunuSamples = znunu100, znunu200, znunu400, znunu600, znunu800, znunu1200, znunu2500
for ds in znunuSamples: ds.xsecs = [ ds.xsecs[0] * 1.23 ]
znunu = sum( znunuSamples )
znunu.label = "Z#rightarrow#nu#nu"


# k-factor from http://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2016/078 EXO-16-014
zg40 = Dataset( "ZNuNuGJets_MonoPhoton_PtG-40to130", 2.789*1.39, ROOT.kMagenta+1, "ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph" )
zg130 = Dataset( "ZNuNuGJets_MonoPhoton_PtG-130", 0.1832*1.39, ROOT.kMagenta+2, "ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph" )
zg = zg40 + zg130
zg.label = "#gammaZ(#nu#nu)"

zgll = Dataset( "ZGTo2LG", 117.864, ROOT.kMagenta+3, "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
zgll.label = "#gammaZ(ll)"


# cross section from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#DY_Z
dy = Dataset( "DYJetsToLL_M-50", 6025.2, ROOT.kRed+3, "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8" )
dy.label = "Z#rightarrowll"


# signal samples

t5wg_1600_100 = Dataset("SMS-T5Wg_1600_100", 0.00810078, ROOT.kRed, "")
t5wg_1600_100.label = "T5Wg 1600 100"
t5wg_1600_1500 = Dataset("SMS-T5Wg_1600_1500", 0.00810078, ROOT.kRed+4, "")
t5wg_1600_1500.label = "T5Wg 1600 1500"

import collections
class SampleCollection(collections.MutableMapping):
    """Dictionary used to store signal datasets.
    The first time a dataset is requested, it is created"""

    colors = [ROOT.kGreen+4, ROOT.kGreen+1] + [ ROOT.kGreen-i for i in range(5) ] + range(1000)

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        if key not in self.store:
            self.store[key] = Dataset( key, col=self.colors[len(self)] )
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

if __name__ == "__main__":
    # print information
    for i in gjets, qcd, ttjets, ttg, wjets, wg_mg, wg_130, znunu,zg: print i.getLatexTableLine(),
    print

    print gjets.getLatexTableHeader()
    for i in gjets40, gjets100, gjets200, gjets400, gjets600, qcd100, qcd200, qcd300, qcd500, qcd700, qcd1000, qcd1500, qcd2000, ttjets, wjets100, wjets200, wjets400, wjets600, wjets800, wjets1200, wjets2500, ttg, wg_mg, wg_130, znunu100, znunu200, znunu400, znunu600, zg:
        print i.getLatexTableLine(),
    print
