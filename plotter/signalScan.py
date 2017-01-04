#!/usr/bin/env python2

from include import *
import DatacardParser

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

def getRgraphs(outputDir, combi, inputData, inputSignal):

    f = ROOT.TFile(inputSignal)
    dirs = [k.GetName() for k in f.GetListOfKeys() if k.GetName().startswith(combi)]
    dirs = [d for d in dirs if "_100" in d or "_200" in d]

    xBins = [350, 600, 700]
    yBins = [0, 1500, 2000]

    bInfo = {
            "bin25": (1,1),
            "bin26": (1,2),
            "bin27": (1,3)
        }

    options, b = DatacardParser.addDatacardParserOptions(optparse.OptionParser())
    dc = limitTools.MyDatacard(DatacardParser.parseCard(file(inputData), options))

    defaultGr = ROOT.TGraph2D(len(dirs))
    graphs = dict( (x,defaultGr.Clone(x)) for x in ["obs","exp","exp1up","exp1dn","exp2up","exp2dn"] )
    for g in graphs.values(): g.SetDirectory(0)
    for idir, d in enumerate(dirs):
        h = f.Get(d+"/met_vs_emht")
        h = aux.rebinX(h, xBins, yBins)
        combi2, m1, m2 = getPointFromDir(d)
        xsec = aux.getXsecSMSglu(m1)
        res = {}
        for bname, (x,y) in bInfo.iteritems():
            acc = h.GetBinContent(x,y)
            accErr = h.GetBinError(x,y)
            res[bname] = (acc*aux.intLumi*xsec, accErr*aux.intLumi*xsec)
        dc.newSignal(res)
        dcName = "{}/{}.txt".format(outputDir, d)
        dc.write(dcName)
        #rInfo = limitTools.infosFromDatacard(dcName)
        rInfo = dc.limitFast()
        for name, gr in graphs.iteritems():
            graphs[name].SetPoint(idir, m1, m2, rInfo[name] )
    return graphs

def getXsecLimitHist( gr2d, h ):
    h.SetDirectory(0) # or the next line will overwrite the hist?
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for x,y,z in points:
        xsec = aux.getXsecSMSglu( x )
        h.SetBinContent(h.FindBin(x,y), z*xsec )
    return h

def writeDict( d, filename ):
    f = ROOT.TFile( filename, "recreate")
    for name, ob in d.iteritems():
        ob.Write(name)
    f.Close()

def getHistForModel( model ):
    if model == "T5gg": return ROOT.TH2F("","", 21, 975, 2025, 18, 50, 1850 )
    if model == "T5Wg": return ROOT.TH2F("","", 16, 775, 1575, 31, -25, 1525 )
    print "Not specified model", model


def calculateLimit(name, combi, inputData, inputSignal):
    outputDir = "limitCalculations/"+name
    if not os.path.isdir(outputDir):
        os.mkdir(outputDir)
        graphs = getRgraphs(outputDir, combi, inputData, inputSignal)
        writeDict(graphs, outputDir+"/Graphs2d.root")

    graphs = readDict(outputDir+"/Graphs2d.root")
    toDraw = dict( [(name,limitTools.getContour(gr)) for name,gr in graphs.iteritems() ] )
    toDraw["obs_hist"] = getXsecLimitHist( graphs["obs"], getHistForModel("T5Wg") )
    if False:
        toDraw["obs_hist"] = interpolateH2( toDraw["obs_hist"] )
        toDraw["obs_hist"] = interpolateH2( toDraw["obs_hist"] )
    print toDraw
    writeDict(toDraw, outputDir+"/Graphs1d.root")

    scanName = "T5Wg"
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName, "plots/%s_limits_"%scanName])




if __name__ == "__main__":
    calculateLimit("T5Wg_v1", "Wg", "limitCalculations/observation_v1.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")

