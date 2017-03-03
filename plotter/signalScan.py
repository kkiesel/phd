#!/usr/bin/env python2

from include import *
import DatacardParser
import multiprocessing
import glob
'''
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

def scaleObsWithXsec(gr2d):
    out = gr2d.Clone(aux.randomName())
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for ip, (x,y,z) in enumerate(points):
        xsec = aux.getXsecSMSglu(x)
        out.SetPoint(ip, x, y, z*xsec)
    return out

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

def recalculateLimits(outputDir):
    for filename in sorted(glob.glob(outputDir+"/*txt")):
        l = limitTools.Limit(filename)
        l.getInfo()
        if l.error: print filename, l.rMinNLL, l.error




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
    if "TChiWG" in inputSignal: proceedWithWeakScan(outputDir, scanName, xsecFile)
    #build2dGraphs(outputDir, combi)
    #graphs = readDict(outputDir+"/Graphs2d.root")
    #toDraw = dict( [(name,limitTools.getContour(gr)) for name,gr in graphs.iteritems() ] )
    #toDraw.update(getObsUncertainty(graphs["obs"], xsecFile))
    #toDraw["obs_hist"] = getXsecLimitHistDelaunay(graphs["obs"])
    #toDraw["obs_hist"] = getXsecLimitHist( graphs["obs"], getHistForModel(scanName), xsecFile )
    #interpolateAlongY(toDraw["obs_hist"])
    #writeDict(toDraw, outputDir+"/Graphs1d.root")

    writeSMSLimitConfig(outputDir+"/Graphs1d.root", "smsPlotter/config/SUS16047/%s_SUS16047.cfg"%scanName)
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS16047/%s_SUS16047.cfg"%scanName, "plots/%s_limits_"%scanName])

def significanceScan(dirName):
    gr = ROOT.TGraph2D()
    for ipoint, f in enumerate(glob.glob("limitCalculations/{}/*.txt".format(dirName))):
        out = limitTools.callCombineSignificance(f)
        for line in out.split("\n"):
            if line.startswith("Limit:"):
                r = float(line.split()[3])
        m = re.match(".*_(\d+)_(\d+).txt", f)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        gr.SetPoint(ipoint, m1, m2, r)
        print out, r


#if __name__ == "__main__":
    #checkConsistency("limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T5Wg_signalScan.root", "../histogramProducer/SMS-T5Wg_1600_100_hists.root")
    #checkConsistency("testDatacard.txt", "../histogramProducer/SMS-T5Wg_signalScan.root", "../histogramProducer/SMS-T5Wg_1600_100_hists.root")

    #signalScan("T5Wg_v7", "Wg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T5gg_v7", "gg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("T6Wg_v7", "Wg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")
    #signalScan("T6gg_v7", "gg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")
    #signalScan("TChiWg_v7", "Wg", "limitCalculations/observation_v4.txt", "../histogramProducer/SMS-TChiWG_signalScan.root")

    #significanceScan("T6gg_v7")

    #getLimit2dHist("../histogramProducer/SMS-T6Wg_signalScan.root")
    #getLimit2dHist("../histogramProducer/SMS-T5Wg_signalScan.root")

###############################################################################
'''
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
    exp2sigma.SetLineColor(exp2sigma.GetFillColor())
    exp1sigma.SetFillColor(ROOT.kGreen+1)
    exp1sigma.SetLineColor(2)
    exp1sigma.SetLineStyle(2)
    exp1sigma.SetLineWidth(2)

    exp2sigma.SetMaximum(2)
    exp2sigma.SetTitle(";m_{#tilde{#chi}_{1}} (GeV);95% CL cross section upper limit (pb)")
    exp2sigma.GetXaxis().SetLimits(300,1300)

    # draw
    exp2sigma.Draw("ap3")
    exp1sigma.Draw("3 same")
    exp1sigma.Draw("xc")

    xsecGr.Draw("xc")
    xsecGrUp.Draw("xc")
    xsecGrDn.Draw("xc")

    obsGr.Draw("xcp")

    # legend
    exp1sigmaClone = exp1sigma.Clone()
    exp1sigmaClone.SetLineColor(exp1sigma.GetFillColor())
    exp1sigmaClone.SetLineStyle(1)
    leg = ROOT.TLegend(.45,.59,.94,.92)
    leg.SetFillStyle(0)
    leg.AddEntry(obsGr, "Observed limit", "l")
    leg.AddEntry(exp2sigma, "Expected limit #pm 1(2) s.d._{experiment}", "f")
    leg.AddEntry(xsecGr, "Signal cross section #pm s.d._{theory}", "l")
    leg.Draw()

    leg2 = ROOT.TLegend(.45,.66,.94,.85)
    leg2.SetFillStyle(0)
    leg2.AddEntry(None, "", "")
    leg2.AddEntry(exp1sigmaClone, "", "f")
    leg2.AddEntry(None, "", "")
    leg2.Draw()

    leg3 = leg.Clone()
    leg3.Clear()
    leg3.AddEntry(None, "", "")
    leg3.AddEntry(exp1sigma, "", "l")
    leg3.AddEntry(None, "", "")
    leg3.Draw()

    leg4 = ROOT.TLegend(.45,.63,.94,.662)
    leg4.SetFillStyle(0)
    leg4.AddEntry(xsecGrUp, "", "l")
    leg4.AddEntry(xsecGrUp, "", "l")
    leg4.Draw()

    lsp_s1 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
    lsp_s2 = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"

    aux.Label(info="pp#rightarrow%s(#tilde{G}#gamma)%s(#tilde{G}W^{#pm})"%(lsp_s1,lsp_s2))
    ROOT.gPad.SetLogy()
    aux.save("TChiWG_limit", log=False)
    exit()

def getPointFromDir(name):
    m = re.match("(.*)_(.*)_(.*)", name)
    combi, m1, m2 = m.groups()
    m1,m2 = int(m1), int(m2)
    return combi, m1, m2

def writeDict(d, filename):
    f = ROOT.TFile(filename, "recreate")
    for name, ob in d.iteritems():
        if ob: ob.Write(name)
    f.Close()

def readDict( filename ):
    f = ROOT.TFile( filename )
    tmpDir = f.GetDirectory( path )
    d = {}
    for element in tmpDir.GetListOfKeys():
        obj = element.ReadObj()
        obj = ROOT.gROOT.CloneObject( obj )
        d[element.GetName()] = obj
    return d

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

def getXsecFile(name):
    xsecFile = ""
    if "T5" in name: xsecFile = "data/xSec_SMS_Gluino_13TeV.pkl"
    elif "T6" in name: xsecFile = "data/xSec_SMS_Squark_13TeV.pkl"
    elif "TChiWG" in name: xsecFile = "data/xSec_SMS_N2C1_13TeV.pkl"
    else: print "Do not know which cross section belongs to", name
    return xsecFile

def getMultiScanName(inputSignal):
    return os.path.basename(inputSignal).split("_")[0][4:]

def getScanName(inputSignal, combi):
    return getMultiScanName(inputSignal).replace("Wg", combi)

def getSignalUncertainties(inputSignal, dirname):
    xBins = [350, 450, 600]
    out = {}
    f = ROOT.TFile(inputSignal)
    hNominal = f.Get(dirname+"/met")
    hNominal = aux.rebinX(hNominal, xBins)

    hGenMet = f.Get(dirname+"/metGen")
    hGenMet = aux.rebinX(hGenMet, xBins)
    out["genMet"] = [1.+abs(hNominal.GetBinContent(b)-hGenMet.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hPuUp = f.Get(dirname+"/met_nopuUp")
    hPuDn = f.Get(dirname+"/met_nopuDn")
    hPuUp.Scale(hNominal.GetEntries()/hPuUp.GetEntries())
    hPuDn.Scale(hNominal.GetEntries()/hPuDn.GetEntries())
    hPuUp = aux.rebinX(hPuUp, xBins)
    hPuDn = aux.rebinX(hPuDn, xBins)
    out["pu"] = [1.+abs(hPuUp.GetBinContent(b)-hPuDn.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hJesUp = f.Get(dirname+"/met_jesUp")
    hJesDn = f.Get(dirname+"/met_jesDn")
    hJesUp = aux.rebinX(hJesUp, xBins)
    hJesDn = aux.rebinX(hJesDn, xBins)
    hJerUp = f.Get(dirname+"/met_jerUp")
    hJerDn = f.Get(dirname+"/met_jerDn")
    hJerUp = aux.rebinX(hJerUp, xBins)
    hJerDn = aux.rebinX(hJerDn, xBins)
    out["jes"] = [1.+math.sqrt((abs(hJesUp.GetBinContent(b)-hJesDn.GetBinContent(b))/2/hNominal.GetBinContent(b))**2 \
        + (abs(hJerUp.GetBinContent(b)-hJerDn.GetBinContent(b))/2/hNominal.GetBinContent(b))**2) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hIsr = f.Get(dirname+"/met_isr")
    hIsrUp = f.Get(dirname+"/met_isrUp")
    hIsrDn = f.Get(dirname+"/met_isrDn")
    hIsr.Scale(hNominal.Integral(0,-1)/hIsr.Integral(0,-1))
    hIsrUp.Scale(hNominal.Integral(0,-1)/hIsrUp.Integral(0,-1))
    hIsrDn.Scale(hNominal.Integral(0,-1)/hIsrDn.Integral(0,-1))
    hIsr = aux.rebinX(hIsr, xBins)
    hIsrUp = aux.rebinX(hIsrUp, xBins)
    hIsrDn = aux.rebinX(hIsrDn, xBins)
    out["isr"] = [1.+abs(hIsrUp.GetBinContent(b)-hIsrDn.GetBinContent(b))/2/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    scaleHists = dict([(i,f.Get(dirname+"/met_weight{}".format(i))) for i in range(1,9)])
    for hname, h in scaleHists.iteritems():
        h.Scale(hNominal.Integral(0,-1)/h.Integral(0,-1))
        scaleHists[hname] = aux.rebin(h, xBins)
    out["scale"] = [1.+(max([h.GetBinContent(b) for h in scaleHists.values()])-min([h.GetBinContent(b) for h in scaleHists.values()]))/2./hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    out["stat"] = [1.+hNominal.GetBinError(b)/hNominal.GetBinContent(b) if hNominal.GetBinContent(b) else 1 for b in range(1,4)]

    hCont = f.Get(dirname+"/met_contamination")
    hCont = aux.rebinX(hCont, xBins)
    cont = [hCont.GetBinContent(b) for b in range(1,4)]
    acc = [hNominal.GetBinContent(b) for b in range(1,4)]

    return acc, cont, out



def writeDataCards(outputDir, inputData, inputSignal, combi="", xsecFile=""):
    f = ROOT.TFile(inputSignal)
    dirs = [k.GetName() for k in f.GetListOfKeys() if k.GetName().startswith(combi)]
    #dirs = ["Wg_1600_100"] # cross check
    binNames = ["binlowEMHT_24", "binlowEMHT_25", "binlowEMHT_26", "binhighEMHT_24", "binhighEMHT_25", "binhighEMHT_26"]

    dc = limitTools.MyDatacard(inputData)
    for d in dirs:
        combi2, m1, m2 = getPointFromDir(d)
        xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        acc, cont, syst = getSignalUncertainties(inputSignal, d+"/signal_lowEMHT")
        acc2, cont2, syst2 = getSignalUncertainties(inputSignal, d+"/signal_highEMHT")
        acc.extend(acc2)
        cont.extend(cont2)
        for a,b in syst.iteritems():
            syst[a] = syst[a]+syst2[a]

        acc = [a*aux.intLumi*xsec for a in acc]
        cont = [a*aux.intLumi*xsec for a in cont]

        acc = [a-b for a,b in zip(acc,cont)]

        systs = {}
        for sName, valueList in syst.iteritems():
            if sName == "stat":
                for binName, unc in zip(binNames, valueList):
                    systs["signalStat_{}".format(binName)] = {binName:unc}
            else:
                systs[sName] = dict(zip(binNames, valueList))

        dc.newSignal( dict(zip(binNames, acc)), systs)
        dcName = "{}/{}.txt".format(outputDir, d)
        dc.write(dcName)

def callMultiCombine(outputDir):
    files = glob.glob("{}/*.txt".format(outputDir))
    p = multiprocessing.Pool()
    p.map(limitTools.callCombine, files)

def callMultiSignificance(outputDir):
    files = glob.glob("{}/*.txt".format(outputDir))
    p = multiprocessing.Pool()
    p.map(limitTools.callCombineSignificance, files)

def clearWrongCombineOutputs(outputDir):
    files = glob.glob("{}/*.txt.limit".format(outputDir))
    for fname in files:
        with open(fname) as f:
            rInfo = limitTools.infoFromOut(f.read())
        if rInfo["obs"] == 0:
            os.remove(fname)
    files = glob.glob("{}/*.txt.significance".format(outputDir))
    for fname in files:
        with open(fname) as f:
            rInfo = limitTools.significanceFromOut(f.read())
        if not rInfo:
            os.remove(fname)

def build2dGraphsLimit(outputDir):
    files = glob.glob("{}/*.txt.limit".format(outputDir))
    defaultGr = ROOT.TGraph2D(len(files))
    graphs = dict( (x,defaultGr.Clone(x)) for x in ["obs","exp","exp1up","exp1dn","exp2up","exp2dn"] )
    for g in graphs.values(): g.SetDirectory(0)
    for ifile, _file in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt.limit", _file)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        with open(_file) as f:
            rInfo = limitTools.infoFromOut(f.read())
        for name, gr in graphs.iteritems():
            graphs[name].SetPoint(ifile, m1, m2, rInfo[name])
    writeDict(graphs, outputDir+"/saved_graphs2d_limit.root")
    return graphs

def build2dGraphsSignificance(outputDir):
    files = glob.glob("{}/*.txt.significance".format(outputDir))
    defaultGr = ROOT.TGraph2D(len(files))
    defaultGr.SetDirectory(0)
    for ifile, _file in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt.significance", _file)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        with open(_file) as f:
            r = limitTools.significanceFromOut(f.read())
        defaultGr.SetPoint(ifile, m1, m2, r)
    writeDict({"significance": defaultGr}, outputDir+"/saved_graphs2d_significance.root")

def latexScanName(scanName):
    if scanName == "T6gg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
    elif scanName == "T6Wg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{q}#tilde{q}, #tilde{q} #rightarrow q%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
    elif scanName == "T5gg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma#tilde{G}"%(lsp_s,lsp_s)
    elif scanName == "T5Wg":
        lsp_s = "#lower[-0.12]{#tilde{#chi}}#lower[0.2]{#scale[0.85]{^{0/#pm}}}#kern[-1.3]{#scale[0.85]{_{1}}}"
        return "pp #rightarrow #tilde{g}#tilde{g}, #tilde{g} #rightarrow qq%s, %s #rightarrow #gamma/W^{#pm}#tilde{G}"%(lsp_s,lsp_s)
    return scanName

def drawSignalScanHist(h, scanName, saveName):
    style.style2d()
    h.SetTitle("")
    if "T5" in scanName: h.GetXaxis().SetTitle("m_{#tilde{g}} (GeV)")
    if "T6" in scanName: h.GetXaxis().SetTitle("m_{#tilde{q}} (GeV)")
    if "gg" in scanName: h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0}_{1}} (GeV)")
    if "Wg" in scanName: h.GetYaxis().SetTitle("m_{#tilde{#chi}^{0/#pm}_{1}} (GeV)")
    hname = h.GetName()
    if hname == "significance":
        h.SetMinimum(-2)
        h.SetMaximum(2)
        style.setPaletteRWB()
        h.GetZaxis().SetTitle("Significance (s.d.)")
    elif hname == "xsec":
        h.Scale(1000)
        h.GetZaxis().SetTitle("signal cross section (fb)")
    else:
        uncerts = {"isr": "ISR uncert", "jes": "Jet energy uncert", "pu": "Pile-up", "scale": "Renorm.+Factor. scale uncert", "genMet": "FastSim p_{T}^{miss}"}
        var, emhtSel, metSel = h.GetName().split("_")
        if var in uncerts.keys():
            h.SetMinimum(0)
            h.SetMaximum(0.1)
        emhtSel = emhtSel.replace("binlowEMHT", "700-2000").replace("binhighEMHT", "2000-#infty")
        metSel = metSel.replace("24", "350-450").replace("25", "450-600").replace("26", "600-#infty")
        var = uncerts[var] if var in uncerts else var
        h.GetZaxis().SetTitle("{} ({},{})".format(var, emhtSel, metSel))
    h.GetZaxis().SetTitleOffset(1.42)
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset()*.95)

    c = ROOT.TCanvas()
    h.Draw("colz")
    aux.Label2D(info="#scale[.76]{{{}}}".format(latexScanName(scanName)))
    aux.save("{}_{}".format(scanName,saveName))
    style.defaultStyle()

def drawLimitInput(outputDir, scanName, xsecFile):
    files = glob.glob("{}/*.txt".format(outputDir))
    dc = limitTools.MyDatacard(files[0])
    defaultGr = ROOT.TGraph2D(len(files))
    graphs = dict((x,defaultGr.Clone(x)) for x in ["Acceptance_"+b for b in dc.bins]+["xsec"])
    uncerts = "isr", "jes", "pu", "scale", "genMet"
    graphs.update(dict([(x+"_"+y,defaultGr.Clone(x+"_"+y)) for x in uncerts for y in dc.bins]))
    for ifile, f in enumerate(files):
        m = re.match(".*_(\d+)_(\d+).txt", f)
        m1 = int(m.group(1))
        m2 = int(m.group(2))
        xsec = aux.getXsecInfoSMS(m1, xsecFile)[0]
        graphs["xsec"].SetPoint(ifile, m1, m2, xsec)
        dc = limitTools.MyDatacard(f)
        systDict = dict([(l[0],l) for l in dc.systs])
        for b in dc.bins:
            graphs["Acceptance_"+b].SetPoint(ifile, m1, m2, dc.exp[b]["signal"]/(xsec*aux.intLumi))
            for unc in uncerts:
                graphs[unc+"_"+b].SetPoint(ifile, m1, m2, systDict[unc][4][b]["signal"]-1)
    for name, gr in graphs.iteritems():
        h = gr.GetHistogram()
        drawSignalScanHist(h, scanName, name)


def getXsecLimitHistDelaunay(gr):
    grScaled = scaleObsWithXsec(gr)
    grScaled.SetNpx(500)
    grScaled.SetNpy(500)
    h = grScaled.GetHistogram()
    h = h.Clone(aux.randomName())
    return h

def build2dGraphs(outputDir, xsecFile):
    build2dGraphsLimit(outputDir)
    build2dGraphsSignificance(outputDir)

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

def interpolateHoles(h2):
    for xbin, ybin in aux.loopH(h2):
        c = h2.GetBinContent(xbin, ybin)
        if abs(c)>1e-5: continue
        cNord = h2.GetBinContent(xbin, ybin+1)
        cSout = h2.GetBinContent(xbin, ybin-1)
        cWest = h2.GetBinContent(xbin+1, ybin)
        cEast = h2.GetBinContent(xbin-1, ybin)
        if cNord and cSout and cWest and cEast:
            h2.SetBinContent(xbin, ybin, sum([cNord,cSout,cWest,cEast])/4.)
    return h2



def getXsecLimitHist( gr2d, h, xsecFile ):
    h.SetDirectory(0) # or the next line will overwrite the hist?
    points = [ (gr2d.GetX()[i], gr2d.GetY()[i], gr2d.GetZ()[i]) for i in range(gr2d.GetN()) ]
    for x,y,z in points:
        xsec = aux.getXsecInfoSMS(x, xsecFile)[0]
        h.SetBinContent(h.FindBin(x,y), z*xsec )
    return h

def getHistForModel( model ):
    h = ROOT.TH2F()
    if "T5" in model: h = ROOT.TH2F("","", 35, 775, 2525, 500, 0, 2500)
    elif "T6" in model: h = ROOT.TH2F("","", 17, 975, 1825, 210, 0, 2100)
    else: print "Not specified model", model
    h.SetMinimum(0)
    return h

def smoothContour(gr, neighbors=5, sigma=.5):
    fgaus = ROOT.TF1("fgaus", "gaus", -10, 10)
    fgaus.SetParameters(1,0,1)
    weights = [fgaus.Eval(i*sigma) for i in range(neighbors)]
    out = gr.Clone(aux.randomName())
    n = gr.GetN()
    Xs = [gr.GetX()[i] for i in range(n)]
    Ys = [gr.GetY()[i] for i in range(n)]
    for i, (x, y) in enumerate(zip(Xs,Ys)):
        pNeigh = min(neighbors, i+1, n-i)
        newX, ws, newY = 0, 0, 0
        for j in range(pNeigh):
            if j:
                newX += (Xs[i-j]+Xs[i+j])*weights[j]
                newY += (Ys[i-j]+Ys[i+j])*weights[j]
                ws += 2*weights[j]
            else:
                newX += x*weights[0]
                newY += y*weights[0]
                ws += weights[0]
        out.SetPoint(i, newX/ws, newY/ws)
    return out

def build1dGraphs(outputDir, xsecFile, scanName):
    graphs = readDict(outputDir+"/saved_graphs2d_limit.root")
    toDraw = dict( [(name,limitTools.getContour(gr)) for name,gr in graphs.iteritems() ] )
    toDraw.update(getObsUncertainty(graphs["obs"], xsecFile))
    toDraw = dict([(n,smoothContour(gr)) for n,gr in toDraw.iteritems()])
    #toDraw["obs_hist"] = getXsecLimitHistDelaunay(graphs["obs"])
    toDraw["obs_hist"] = getXsecLimitHist( graphs["obs"], getHistForModel(scanName), xsecFile )
    interpolateAlongY(toDraw["obs_hist"])
    writeDict(toDraw, outputDir+"/saved_graphs1d_limit.root")

def drawSignificance(outputDir, scanName):
    graphs = readDict(outputDir+"/saved_graphs2d_significance.root")
    h = graphs["significance"].GetHistogram()
    h.SetName("significance")
    interpolateHoles(h)
    drawSignalScanHist(h, scanName, "significance")

def signalScan(combi, version, inputData, inputSignal):
    scanName = getScanName(inputSignal, combi)
    outputDir = "limitCalculations/{}_{}".format(scanName, version)
    if not os.path.isdir(outputDir): os.mkdir(outputDir)
    xsecFile = getXsecFile(inputSignal)
    #writeDataCards(outputDir, inputData, inputSignal, combi, xsecFile)
    #drawLimitInput(outputDir, scanName, xsecFile)
    #callMultiCombine(outputDir)
    #callMultiSignificance(outputDir)
    ##clearWrongCombineOutputs(outputDir)
    if "TChiWG" in inputSignal: proceedWithWeakScan(outputDir, scanName, xsecFile)
    #build2dGraphs(outputDir, xsecFile)
    build1dGraphs(outputDir, xsecFile, scanName)
    #drawSignificance(outputDir, scanName)
    writeSMSLimitConfig(outputDir+"/saved_graphs1d_limit.root", "smsPlotter/config/SUS16047/%s_SUS16047.cfg"%scanName)
    subprocess.call(["python2", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS16047/%s_SUS16047.cfg"%scanName, "plots/%s_limits_"%scanName])


if __name__ == "__main__":
    #signalScan("Wg", "v9", "limitCalculations/observation_v6.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    #signalScan("gg", "v9", "limitCalculations/observation_v6.txt", "../histogramProducer/SMS-T5Wg_signalScan.root")
    signalScan("Wg", "v9", "limitCalculations/observation_v6.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")
    #signalScan("gg", "v9", "limitCalculations/observation_v6.txt", "../histogramProducer/SMS-T6Wg_signalScan.root")

    #signalScan("Wg", "v9", "limitCalculations/observation_v6.txt", "../histogramProducer/SMS-T5Wg_mGo2150To2500_signalScan.root")
    #signalScan("gg", "v9", "limitCalculations/observation_v6.txt", "../histogramProducer/SMS-T5Wg_mGo2150To2500_signalScan.root")


