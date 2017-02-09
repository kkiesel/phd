#!/usr/bin/env python2

from include import *
import DatacardParser
import multiprocessing
import glob


def getPointFromDir(name):
    m = re.match("(.*)_(.*)_(.*)", name)
    combi, m1, m2 = m.groups()
    m1,m2 = int(m1), int(m2)
    return combi, m1, m2

def readDict( filename ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )
    d = {}
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()
        obj = ROOT.gROOT.CloneObject( obj )
        d[element.GetName()] = obj
    return d

def scaleObsWithXsec(gr2d):
    out = gr2d.Clone(aux.randomName())
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for ip, (x,y,z) in enumerate(points):
        xsec = aux.getXsecSMSglu(x)
        out.SetPoint(ip, x, y, z*xsec)
    return out

def getXsecLimitHist( gr2d, h ):
    h.SetDirectory(0) # or the next line will overwrite the hist?
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for x,y,z in points:
        xsec = aux.getXsecSMSglu(x)
        h.SetBinContent(h.FindBin(x,y), z*xsec )
    return h

def getXsecLimitHistDelaunay(gr):
    grScaled = scaleObsWithXsec(gr)
    grScaled.SetNpx(500)
    grScaled.SetNpy(500)
    h = grScaled.GetHistogram()
    h = h.Clone(aux.randomName())
    return h

def getObsUncertainty(gr2d, xsecFile):
    gr2dup = gr2d.Clone(aux.randomName())
    gr2ddn = gr2d.Clone(aux.randomName())
    gr2dup.SetDirectory(0)
    gr2ddn.SetDirectory(0)
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for ip, (x,y,z) in enumerate(points):
        xsec, unc = aux.getXsecInfoSMS(x, xsecFile)
        gr2dup.SetPoint(ip, x, y, z*(1-unc/100))
        gr2ddn.SetPoint(ip, x, y, z*(1+unc/100))
    obsUp = limitTools.getContour(gr2dup)
    obsUp.SetName("obs1up")
    obsDn = limitTools.getContour(gr2ddn)
    obsDn.SetName("obs1dn")
    return {"obs1up": obsUp, "obs1dn": obsDn}

def writeDict( d, filename ):
    f = ROOT.TFile( filename, "recreate")
    for name, ob in d.iteritems():
        if ob:
            ob.Write(name)
    f.Close()

def getHistForModel( model ):
    h = ROOT.TH2F()
    if "T5" in model: h = ROOT.TH2F("","", 54, 775, 2125, 420, 0, 2100)
    elif "T6" in model: h = ROOT.TH2F("","", 17, 975, 1825, 210, 0, 2100)
    else: print "Not specified model", model
    h.SetMinimum(0)
    return h

def writeSMSLimitConfig(infile, configName):
    text = """
HISTOGRAM {0} obs_hist
EXPECTED {0} exp exp1up exp1dn kRed kOrange
OBSERVED {0} obs obs1up obs1dn kBlack kGray
PRELIMINARY PrivateWork
LUMI {1:.2f}
ENERGY 13
""".format(infile,aux.intLumi/1e3)
    #t2 = "OBSERVED ../../master/singlePhoton/PlotsSMS/config/SUS14004/2015-01-09-limits/SMS_T5wg/ROOT/SMS_T5wg_gluino_chi1_Exclusion_witXsecLimit.root Expected_limit Expected_limit_up Expected_limit_dn kBlack kGray"
    #text = "\n".join([l for l in text.split("\n") if "OBSERVED" not in l]+[t2, ""])
    with open(configName, "w+") as f:
        f.write(text)


def writeDataCards(outputDir, inputData, inputSignal, combi="", xsecFile=""):
    f = ROOT.TFile(inputSignal)
    dirs = [k.GetName() for k in f.GetListOfKeys() if k.GetName().startswith(combi)]
    #dirs = ["Wg_1600_100"] # cross check

    xBins = [350, 450, 600]

    bInfo = {
            "binlowEMHT_24": ("lowEMHT",1),
            "binlowEMHT_25": ("lowEMHT",2),
            "binlowEMHT_26": ("lowEMHT",3),
            "binhighEMHT_24": ("highEMHT",1),
            "binhighEMHT_25": ("highEMHT",2),
            "binhighEMHT_26": ("highEMHT",3),
        }

    dc = limitTools.MyDatacard(inputData)
    for d in dirs:
        combi2, m1, m2 = getPointFromDir(d)
        xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        hists = {
            "lowEMHT": f.Get(d+"/signal_lowEMHT/met"),
            "highEMHT": f.Get(d+"/signal_highEMHT/met"),
        }
        for name, h in hists.iteritems():
            h.Scale(aux.intLumi*xsec)
            hists[name] = aux.rebinX(h, xBins)

        res = {}
        for bname, (name,bin) in bInfo.iteritems():
            c = hists[name].GetBinContent(bin)
            err = hists[name].GetBinError(bin)
            res[bname] = (c, 1.+err/c) if err else (0, 1)
        dc.newSignal(res)
        dcName = "{}/{}.txt".format(outputDir, d)
        dc.write(dcName)

def callMultiCombine(outputDir, combi):
    files = glob.glob("{}/{}_*.txt".format(outputDir, combi))
    p = multiprocessing.Pool()
    p.map(limitTools.callCombine, files)

def build2dGraphs(outputDir, combi):
    files = glob.glob("{}/{}_*.txt.limit".format(outputDir, combi))

    defaultGr = ROOT.TGraph2D(len(files))
    graphs = dict( (x,defaultGr.Clone(x)) for x in ["obs","exp","exp1up","exp1dn","exp2up","exp2dn"] )
    for g in graphs.values(): g.SetDirectory(0)
    for ifile, _file in enumerate(files):
        m = re.match(".*{}_(\d+)_(\d+).txt.limit".format(combi), _file)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        with open(_file) as f:
            rInfo = limitTools.infoFromOut(f.read())
        for name, gr in graphs.iteritems():
            graphs[name].SetPoint(ifile, m1, m2, rInfo[name] )
    writeDict(graphs, outputDir+"/Graphs2d.root")
    return graphs

def signalScan(name, combi, inputData, inputSignal):
    outputDir = "limitCalculations/"+name
    if not os.path.isdir(outputDir): os.mkdir(outputDir)
    writeDataCards(outputDir, inputData, inputSignal, combi)
    callMultiCombine(outputDir, combi)
    build2dGraphs(outputDir, combi)
    graphs = readDict(outputDir+"/Graphs2d.root")
    toDraw = dict( [(name,limitTools.getContour(gr)) for name,gr in graphs.iteritems() ] )
    toDraw.update(getObsUncertainty(graphs["obs"]))
    toDraw["obs_hist"] = getXsecLimitHistDelaunay(graphs["obs"])
    toDraw["obs_hist"] = getXsecLimitHist( graphs["obs"], getHistForModel("T5Wg") )
    writeDict(toDraw, outputDir+"/Graphs1d.root")

    scanName = "T5Wg"
    writeSMSLimitConfig(outputDir+"/Graphs1d.root", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName)
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName, "plots/%s_limits_"%scanName])

def checkHistogramConsistency(h1, h2, relTolerance=1e-6):
    n1 = h1.GetEntries()
    n2 = h2.GetEntries()
    if n1 != n2: print "Number of entries not the same", n1, n2
    nx1 = h1.GetNbinsX()
    nx2 = h2.GetNbinsX()
    ny1 = h1.GetNbinsY()
    ny2 = h2.GetNbinsY()
    if nx1 != nx2: print "Number of x bins is not the same:", nx1, nx2
    if ny1 != ny2: print "Number of y bins is not the same:", ny1, ny2
    for x,y in aux.loopH(h1):
        c1 = h1.GetBinContent(x,y)
        c2 = h2.GetBinContent(x,y)
        if c1+c2 and abs(c1-c2)/(c1+c2)/2 > relTolerance: print "Not same bin content in bin {}:{}:".format(x,y), c1, c2


def checkConsistency(datacardFile, signalScan, treeFile):
    hScan = aux.getFromFile(signalScan, "Wg_1600_100/met_vs_emht")
    hScan.Scale(aux.intLumi*aux.getXsecSMSglu(1600))
    hPlot = aux.getFromFile(treeFile, "tr/met_vs_emht")
    hPlot = t5wg_1600_100.getHist("tr/met_vs_emht")
    print "Plot vs scan differences"
    sameHists = checkHistogramConsistency(hScan, hPlot, 1e-2)

    tree = ROOT.TChain("tr/simpleTree")
    tree.AddFile(treeFile)
    hTree = hPlot.Clone("hTree")
    hTree.Reset("ICEM")
    tree.Draw("emht:met>>hTree", "weight", "goff")
    print "scan versus tree differences"
    sameHists2 = checkHistogramConsistency(hScan, hTree, 1e-2)

    print "check rebinned hists"
    xBins = [350, 450, 600]
    yBins = [0, 2000]

    bInfo = {
            "binemht2000_24": (1,1),
            "binemht2000_25": (2,1),
            "binemht2000_26": (3,1),
            "bin2000emht_24": (1,2),
            "bin2000emht_25": (2,2),
            "bin2000emht_26": (3,2),
        }

    hPlotRebinned = aux.rebinX(hPlot, xBins, yBins)
    hScanRebinned = aux.rebinX(hScan, xBins, yBins)
    sameHists = checkHistogramConsistency(hScanRebinned, hPlotRebinned, 1e-5)


if __name__ == "__main__":
    checkConsistency("limitCalculations/observation_v2.txt", "../histogramProducer/SMS-T5Wg_signalScan.root", "../histogramProducer/SMS-T5Wg_1600_100_hists.root")

    #signalScan("T5Wg_v4", "Wg", "limitCalculations/observation_v1.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T5Wg_v5", "Wg", "limitCalculations/observation_v2.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T5Wg_v5", "gg", "limitCalculations/observation_v2.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T5Wg_v5", "WW", "limitCalculations/observation_v2.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T5Wg_v5", "Wg", "limitCalculations/observation_v2.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")
    #signalScan("T5Wg_v5", "gg", "limitCalculations/observation_v2.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")
    #signalScan("T5Wg_v5", "WW", "limitCalculations/observation_v2.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")

