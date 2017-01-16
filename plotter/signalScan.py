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

def getXsecLimitHist( gr2d, h ):
    h.SetDirectory(0) # or the next line will overwrite the hist?
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for x,y,z in points:
        xsec = aux.getXsecSMSglu(x)
        h.SetBinContent(h.FindBin(x,y), z*xsec )
    return h

def writeDict( d, filename ):
    f = ROOT.TFile( filename, "recreate")
    for name, ob in d.iteritems():
        if ob:
            ob.Write(name)
    f.Close()

def getHistForModel( model ):
    h = ROOT.TH2F()
    if model == "T5gg": h = ROOT.TH2F("","", 21, 975, 2025, 18, 50, 1850 )
    elif model == "T5Wg": h = ROOT.TH2F("","", 16, 775, 1575, 31, -25, 1525 )
    else: print "Not specified model", model
    h.SetMinimum(0)
    return h

def writeSMSLimitConfig(infile, configName):
    text = """
HISTOGRAM {0} obs_hist
EXPECTED {0} exp exp1up exp1dn kRed kOrange
OBSERVED ../../master/singlePhoton/PlotsSMS/config/SUS14004/2015-01-09-limits/SMS_T5wg/ROOT/SMS_T5wg_gluino_chi1_Exclusion_witXsecLimit.root Expected_limit Expected_limit_up Expected_limit_dn kBlack kGray
PRELIMINARY PrivateWork
LUMI {1:.2f}
ENERGY 13
""".format(infile,aux.intLumi/1e3)
    with open(configName, "w+") as f:
        f.write(text)


def writeDataCards(outputDir,inputData, inputSignal, combi=""):
    f = ROOT.TFile(inputSignal)
    dirs = [k.GetName() for k in f.GetListOfKeys() if k.GetName().startswith(combi)]
    #dirs = ["Wg_1600_100"] # cross check

    xBins = [350, 600, 700]
    yBins = [0, 1500, 2000]

    bInfo = {
            "bin25": (1,3),
            "bin26": (2,3),
            "bin27": (3,3)
        }

    options, b = DatacardParser.addDatacardParserOptions(optparse.OptionParser())
    dc = limitTools.MyDatacard(DatacardParser.parseCard(file(inputData), options))

    for d in dirs:
        h = f.Get(d+"/met_vs_emht")
        h = aux.rebinX(h, xBins, yBins)
        combi2, m1, m2 = getPointFromDir(d)
        xsec = aux.getXsecSMSglu(m1)
        h.Scale(aux.intLumi*xsec)
        res = {}
        for bname, (x,y) in bInfo.iteritems():
            c = h.GetBinContent(x,y)
            err = h.GetBinError(x,y)
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
    toDraw["obs_hist"] = getXsecLimitHist( graphs["obs"], getHistForModel("T5Wg") )
    toDraw["obs_hist"] = graphs["obs"].GetHistogram()
    if False:
        toDraw["obs_hist"] = interpolateH2( toDraw["obs_hist"] )
        toDraw["obs_hist"] = interpolateH2( toDraw["obs_hist"] )
    writeDict(toDraw, outputDir+"/Graphs1d.root")

    scanName = "T5Wg"
    writeSMSLimitConfig(outputDir+"/Graphs1d.root", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName)
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName, "plots/%s_limits_"%scanName])


if __name__ == "__main__":
    signalScan("T5Wg_v4", "Wg", "limitCalculations/observation_v1.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")

