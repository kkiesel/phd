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

def getObsUncertainty(gr2d):
    gr2dup = gr2d.Clone(aux.randomName())
    gr2ddn = gr2d.Clone(aux.randomName())
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for ip, (x,y,z) in enumerate(points):
        xsec, unc = aux.getXsecInfoSMS(x, "data/xSec_SMS_Gluino_13TeV.pkl")
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
    t2 = "OBSERVED ../../master/singlePhoton/PlotsSMS/config/SUS14004/2015-01-09-limits/SMS_T5wg/ROOT/SMS_T5wg_gluino_chi1_Exclusion_witXsecLimit.root Expected_limit Expected_limit_up Expected_limit_dn kBlack kGray"
    text = "\n".join([l for l in text.split("\n") if "OBSERVED" not in l]+[t2, ""])
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
    toDraw.update(getObsUncertainty(graphs["obs"]))
    toDraw["obs_hist"] = getXsecLimitHistDelaunay(graphs["obs"])
    toDraw["obs_hist"] = getXsecLimitHist( graphs["obs"], getHistForModel("T5Wg") )
    writeDict(toDraw, outputDir+"/Graphs1d.root")

    scanName = "T5Wg"
    writeSMSLimitConfig(outputDir+"/Graphs1d.root", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName)
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName, "plots/%s_limits_"%scanName])


if __name__ == "__main__":
    signalScan("T5Wg_v4", "Wg", "limitCalculations/observation_v1.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")

