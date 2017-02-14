#!/usr/bin/env python2

from include import *
import DatacardParser
import multiprocessing
import glob

###############################################################################
# check consistency
###############################################################################
def checkHistogramConsistency2d(h1, h2, relTolerance=1e-6):
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

def checkHistogramConsistency1d(h1, h2, relTolerance=1e-6):
    n1 = h1.GetEntries()
    n2 = h2.GetEntries()
    if n1 != n2: print "Number of entries not the same", n1, n2
    nx1 = h1.GetNbinsX()
    nx2 = h2.GetNbinsX()
    if nx1 != nx2: print "Number of x bins is not the same:", nx1, nx2
    for x in aux.loopH(h1):
        c1 = h1.GetBinContent(x)
        c2 = h2.GetBinContent(x)
        if c1+c2 and abs(c1-c2)/(c1+c2)/2 > relTolerance: print "Not same bin content in bin {}:".format(x), c1, c2


def checkConsistency(datacardFile, signalScan, treeFile):
    hScan = aux.getFromFile(signalScan, "Wg_1600_100/signal_lowEMHT/met")
    hScan.Scale(aux.intLumi*aux.getXsecSMSglu(1600))
    #hPlot = aux.getFromFile(treeFile, "signal_lowEMHT/met")
    hPlot = t5wg_1600_100.getHist("signal_lowEMHT/met")
    print "Plot vs scan differences"
    sameHists = checkHistogramConsistency1d(hScan, hPlot, 1e-3)

    xBins = [350, 450, 600]
    hScanRebinned = aux.rebinX(hScan, xBins)
    dc = limitTools.MyDatacard(datacardFile)
    print "Datacard versus scan:"
    print dc.exp["binlowEMHT_24"]["signal"], hScanRebinned.GetBinContent(1)
    print dc.exp["binlowEMHT_25"]["signal"], hScanRebinned.GetBinContent(2)
    print dc.exp["binlowEMHT_26"]["signal"], hScanRebinned.GetBinContent(3)

    return

###############################################################################
# end check consistency
###############################################################################

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

def interpolateAlongY(h2):
    for xbin, ybin in aux.loopH(h2):
        if not h2.GetBinContent(xbin,ybin):
            ybinUp, cUp = 0, 0
            for ybinUp in range(ybin,h2.GetNbinsY()):
                cUp = h2.GetBinContent(xbin,ybinUp)
                if cUp: break
            ybinDn, cDn = 0, 0
            for ybinDn in range(ybin,0, -1):
                cDn = h2.GetBinContent(xbin,ybinDn)
                if cDn: break
            if not cUp or not cDn: continue
            h2.SetBinContent(xbin,ybin, ((cUp-cDn)*ybin+(cDn*ybinUp-cUp*ybinDn))/(ybinUp-ybinDn)) # linear interpolation

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

def getLimit2dHist(inputSignal):
    f = ROOT.TFile(inputSignal)
    dirs = [k.GetName() for k in f.GetListOfKeys() if "Tree" not in k.GetName()]
    xValues = set()
    yValues = set()
    for d in dirs:
        c, m1, m2 = getPointFromDir(d)
        xValues.add(m1)
        yValues.add(m2)
    xValues = sorted(xValues)
    yValues = sorted(yValues)
    print xValues
    print yValues

def getHistForModel( model ):
    h = ROOT.TH2F()
    if "T5" in model: h = ROOT.TH2F("","", 54, 775, 2125, 420, 0, 2100)
    elif "T6" in model: h = ROOT.TH2F("","", 17, 975, 1825, 210, 0, 2100)
    else: print "Not specified model", model
    h.SetMinimum(0)
    return h

def proceedWithWeakScan(outputDir, scanName, xsecFile):
    scanRes = {}
    for fname in glob.glob("{}/*.txt.limit".format(outputDir)):
        m = re.match(".*_(\d+)_(\d+).txt.limit", fname)
        with open(fname) as f: scanRes[int(m.group(1))] = limitTools.infoFromOut(f.read())

    defaultGr = ROOT.TGraph(len(scanRes))
    graphs = dict( (x,defaultGr.Clone(x)) for x in ["obs","exp","exp1up","exp1dn","exp2up","exp2dn"] )
    obsGr = ROOT.TGraph()
    xsecGr = ROOT.TGraph()
    xsecGrUp = ROOT.TGraph()
    xsecGrDn = ROOT.TGraph()
    exp1sigma = ROOT.TGraphAsymmErrors()
    exp2sigma = ROOT.TGraphAsymmErrors()
    for i, m in enumerate(sorted(scanRes)):
        xsec, xsec_unc = aux.getXsecInfoSMS(m, xsecFile)
        xsecGr.SetPoint(i, m, xsec)
        xsecGrUp.SetPoint(i, m, xsec*(1+xsec_unc))
        xsecGrDn.SetPoint(i, m, xsec*(1-xsec_unc))
        obsGr.SetPoint(i, m, xsec*scanRes[m]["obs"])
        expR = scanRes[m]["exp"]
        exp1sigma.SetPoint(i, m, xsec*expR)
        exp2sigma.SetPoint(i, m, xsec*expR)
        exp1sigma.SetPointEYhigh(i, xsec*(scanRes[m]["exp1up"]-expR))
        exp2sigma.SetPointEYhigh(i, xsec*(scanRes[m]["exp2up"]-expR))
        exp1sigma.SetPointEYlow(i, xsec*(expR-scanRes[m]["exp1dn"]))
        exp2sigma.SetPointEYlow(i, xsec*(expR-scanRes[m]["exp2dn"]))
        for name in graphs:
            graphs[name].SetPoint(i, m, xsec*scanRes[m][name] )
    writeDict(graphs, outputDir+"/Graphs2d.root")

    # beautify
    obsGr.SetLineWidth(2)

    for g in xsecGr, xsecGrUp, xsecGrDn:
        g.SetLineColor(ROOT.kBlue)
    xsecGrUp.SetLineStyle(2)
    xsecGrDn.SetLineStyle(2)

    exp2sigma.SetFillColor(ROOT.kOrange)
    exp1sigma.SetFillColor(ROOT.kGreen+1)
    exp1sigma.SetLineColor(2)
    exp1sigma.SetLineStyle(2)
    exp1sigma.SetLineWidth(2)

    exp2sigma.SetMaximum(2)
    exp2sigma.SetTitle(";m_{#tilde{#chi}_{1}} (GeV);95% CL cross section upper limit (pb)")
    exp2sigma.GetXaxis().SetRangeUser(400,1400)

    # draw
    exp2sigma.Draw("ap3")
    exp1sigma.Draw("3 same")
    exp1sigma.Draw("xc")

    xsecGr.Draw("xc")
    xsecGrUp.Draw("xc")
    xsecGrDn.Draw("xc")

    obsGr.Draw("xcp")

    # legend
    leg = ROOT.TLegend(.56,.59,.94,.92)
    leg.SetFillColor( ROOT.kWhite )
    leg.SetFillStyle(0)
    leg.AddEntry(exp1sigma, "Expected limit", "l")
    leg.AddEntry(exp1sigma, "#pm1 s.d.", "f")
    leg.AddEntry(exp2sigma, "#pm 2 s.d.", "f")
    leg.AddEntry(obsGr, "Observed limit", "l")
    leg.AddEntry(xsecGr, "Signal cross section", "l")
    leg.Draw()

    aux.Label(info="TChiWG")
    ROOT.gPad.SetLogy()
    aux.save("TChiWG_limit", log=False)
    exit()


def writeSMSLimitConfig(infile, configName):
    text = """
HISTOGRAM {0} obs_hist
EXPECTED {0} exp exp1up exp1dn kRed kOrange
OBSERVED {0} obs obs1up obs1dn kBlack kGray
PRELIMINARY
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

def recalculateLimits(outputDir):
    for filename in sorted(glob.glob(outputDir+"/*txt")):
        l = limitTools.Limit(filename)
        l.getInfo()
        if l.error: print filename, l.rMinNLL, l.error


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

def getScanName(inputSignal, combi):
    m = re.match(".*/SMS-(..)Wg.*", inputSignal)
    if m: return m.group(1) + combi
    if "TChiWG" in inputSignal: return "TChiWG"
    print "Do not know what to do"


def signalScan(dirName, combi, inputData, inputSignal):
    outputDir = "limitCalculations/"+dirName
    scanName = getScanName(inputSignal, combi)
    if not os.path.isdir(outputDir): os.mkdir(outputDir)
    if "T5" in inputSignal:
        xsecFile = "data/xSec_SMS_Gluino_13TeV.pkl"
    elif "T6" in inputSignal:
        xsecFile = "data/xSec_SMS_Squark_13TeV.pkl"
    elif "TChiWG" in inputSignal:
        xsecFile = "data/xSec_SMS_N2C1_13TeV.pkl"
    else:
        print "Do not know if squark or gluino scan"
        return

    #writeDataCards(outputDir, inputData, inputSignal, combi, xsecFile)
    #callMultiCombine(outputDir, combi)
    if "TChiWG" in inputSignal: proceedWithWeakScan(outputDir, scanName)
    build2dGraphs(outputDir, combi)
    graphs = readDict(outputDir+"/Graphs2d.root")
    toDraw = dict( [(name,limitTools.getContour(gr)) for name,gr in graphs.iteritems() ] )
    toDraw.update(getObsUncertainty(graphs["obs"], xsecFile))
    #toDraw["obs_hist"] = getXsecLimitHistDelaunay(graphs["obs"])
    toDraw["obs_hist"] = getXsecLimitHist( graphs["obs"], getHistForModel(scanName) )
    interpolateAlongY(toDraw["obs_hist"])
    writeDict(toDraw, outputDir+"/Graphs1d.root")

    writeSMSLimitConfig(outputDir+"/Graphs1d.root", "smsPlotter/config/SUS16047/%s_SUS16047.cfg"%scanName)
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS16047/%s_SUS16047.cfg"%scanName, "plots/%s_limits_"%scanName])

if __name__ == "__main__":
    #checkConsistency("limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T5Wg_signalScan.root", "../histogramProducer/SMS-T5Wg_1600_100_hists.root")
    #checkConsistency("testDatacard.txt", "../histogramProducer/SMS-T5Wg_signalScan.root", "../histogramProducer/SMS-T5Wg_1600_100_hists.root")

    #signalScan("T5Wg_v7", "Wg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T5gg_v7", "gg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T6Wg_v7", "Wg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")
    #signalScan("T6gg_v7", "gg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")
    signalScan("TChiWg_v7", "Wg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-TChiWG_signalScan.root")

    #getLimit2dHist("../histogramProducer/SMS-T6Wg_signalScan.root")


