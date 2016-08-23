#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import re
if sys.version_info[:2] == (2,6):
    print "Initialize correct python version first!"
    sys.exit()

# private libs
import ratio
import style
import multiplot
from rwthColors import rwth
import auxiliary as aux
from datasets import *

folder2tex = {
    "tr": "#gamma(EB)",
    "tr_tight": "#gamma (tight)",
    "tr_hltHT600": "#gamma #it{H}_{T}-trigger",
    "tr_0met100": "#gamma #it{E}_{T}^{miss}<100GeV",
    "tr_100met": "#gamma #it{E}_{T}^{miss}>100GeV",
    "tr_ee": "#gamma(EE)",
    "tr_jControl": "jet",
    "tr_jControl_0met100": "jet #it{E}_{T}^{miss}<100GeV",
    "tr_jControl_100met": "jet #it{E}_{T}^{miss}>100GeV",
    "tr_jControl_ee": "jet(EE)",
}

def getNames(dataset=data, folder="tr", objects=[ROOT.TH1F]):
    return aux.getObjectNames(dataset.files[0], folder, objects)

def getBinningAndNames(names=None):
    if not names: names = getNames()
    return [ (n, bn, b) for n in names for bn, b in aux.getBinningsFromName(n).iteritems() ]

def getStyledHist(dataset, name, binning, directory="tr"):
    h = dataset.getHist("{}/{}".format(directory,name))
    if not h: return
    if binning: h = aux.rebin(h, binning)
    aux.appendFlowBin(h)
    if style.divideByBinWidth: h.Scale(1, "width")
    h.SetYTitle(aux.getYAxisTitle(h))
    return h

def integral(h):
    binOption = "width" if style.divideByBinWidth else ""
    return h.Integral(0, -1, binOption)

def compareDirs(sampleName, dir1, dir2, dataset1, dataset2=None):
    if not dataset2: dataset2 = dataset1
    for name, binningName, binning in getBinningAndNames():
        h1 = getStyledHist(dataset1, name, binning, dir1)
        h2 = getStyledHist(dataset2, name, binning, dir2)
        i1 = integral(h1)
        i2 = integral(h2)
        if not i1 or not i2: continue
        aux.drawOpt(h1, "data")
        aux.drawOpt(h2, "pre")
        h2.Scale(i1/i2)

        h2.drawOption_="hist e"

        can = ROOT.TCanvas()
        m = multiplot.Multiplot()
        m.add(h1, "{} {}".format(folder2tex[dir1], aux.metricPrefix(i1)))
        m.add(h2, "{} {}".format(folder2tex[dir2], aux.metricPrefix(i2)))

        additional = {}
        #if name == "met": additional = {"metUp":None, "metDn":None}
        #if name == "metPar": additional = {"metParUp":None, "metParDn":None}
        for n in additional.keys():
            h = getStyledHist(dataset1, n, binning, dir2)
            aux.drawOpt(h, "pre")
            if not integral(h): continue
            h.Scale(i1/integral(h))
            h.drawOption_="hist"
            #m.add(h, n)
            additional[n] = h
        syst = aux.getSystFromDifference(additional.values()[0], additional.values()[1]) if additional else None
        if syst: m.add(syst, "#pm #sigma_{p_{T}}")
        m.Draw()
        l = aux.Label(info=dataset1.label, sim=dataset1 is not data)
        r = ratio.Ratio("{}/{}".format(folder2tex[dir1], folder2tex[dir2]), h1, h2, syst)
        r.draw(.5, 1.5)

        funcName = sys._getframe().f_code.co_name
        if binningName: binningName = "_" + binningName
        saveName = "{}_{}_{}_vs_{}_{}{}".format(funcName, sampleName, dir1, dir2, name, binningName)
        aux.save(saveName)

def drawSameHistogram(sampleName, folder="tr", stack=[], additional=[]):
    stack.reverse()
    for name, binningName, binning in getBinningAndNames():

        can = ROOT.TCanvas()
        m = multiplot.Multiplot()

        for d in stack:
            h = getStyledHist(d, name, binning, folder)
            m.addStack(h, d.label)

        dataHist = None
        for d in additional:
            h = getStyledHist(d, name, binning, folder)
            if d == data:
                aux.drawOpt(h, "data")
                dataHist = h
            else: aux.drawOpt(h, "signal")
            m.add(h, d.label)

        m.sortStackByIntegral()
        m.Draw()
        #hsm = m.hists[0].GetStack().Last()

        if binningName: binningName = "_"+binningName
        saveName = "sameHistograms_{}_{}__{}{}".format(sampleName, folder, name, binningName)
        aux.save(saveName)

def weighting2d(h2num, h2den, xCut=100, savename=""):
    xCutBin = h2num.GetXaxis().FindBin(xCut-1e-6)
    h1numY = h2num.ProjectionY(aux.randomName(), 0, xCutBin)
    h1denY = h2den.ProjectionY(aux.randomName(), 0, xCutBin)
    h1numY.Divide(h1denY)

    if savename:
        aux.write2File(h1numY, savename, "weights.root")

    h2denW = h2den.Clone(aux.randomName())
    h2denWsys = h2den.Clone(aux.randomName())
    for xbin, ybin in aux.loopH2(h2den):
        w = h1numY.GetBinContent(ybin)
        we = h1numY.GetBinError(ybin)
        c = h2den.GetBinContent(xbin,ybin)
        h2denW.SetBinContent(xbin, ybin, c*w )
        h2denW.SetBinError(xbin, ybin, h2den.GetBinError(xbin, ybin) * w )
        h2denWsys.SetBinContent(xbin, ybin, c*w )
        h2denWsys.SetBinError(xbin, ybin, c*we )
    return h2denW, h2denWsys

def finalPrediction(allSets, simpleScale=False, suffix=""):
    jetSet = dataHt if allSets == data else allSets

    hName = "metPer_vs_emht"
    axVars = dict(zip(["x","y"],re.match("(.*)_vs_(.*)",hName).groups()))

    metBinning = aux.getBinningsFromName("met")["3"]
    emhtBinning = aux.getBinningsFromName("emht")["2"]


    h2data = allSets.getHist("tr{}/{}".format(suffix,hName))
    h2jetControl = jetSet.getHist("tr_jControl{}/{}".format(suffix,hName))
    h2eControl = allSets.getHist("tr_eControl{}/{}".format(suffix,hName))
    if True:
        h2eControl = allSets.getHist("tr{}/{}".format(suffix,hName))
        h2eControl.Scale(0)
    h2ttg = ttg.getHist("tr{}/{}".format(suffix,hName))
    h2tg = tg.getHist("tr{}/{}".format(suffix,hName))
    h2wg = wg_mg.getHist("tr{}/{}".format(suffix,hName))
    h2zg = zg_130.getHist("tr{}/{}".format(suffix,hName))
    h2dy = zgll.getHist("tr{}/{}".format(suffix,hName))
    h2s1 = signal["T5Wg_1550_100"].getHist("tr{}/{}".format(suffix,hName))
    h2s2 = signal["T5Wg_1550_1500"].getHist("tr{}/{}".format(suffix,hName))

    h2data = aux.rebin2d(h2data, metBinning, emhtBinning)
    h2jetControl = aux.rebin2d(h2jetControl, metBinning, emhtBinning)
    h2eControl = aux.rebin2d(h2eControl, metBinning, emhtBinning)
    h2ttg = aux.rebin2d(h2ttg, metBinning, emhtBinning)
    h2tg = aux.rebin2d(h2tg, metBinning, emhtBinning)
    h2wg = aux.rebin2d(h2wg, metBinning, emhtBinning)
    h2zg = aux.rebin2d(h2zg, metBinning, emhtBinning)
    h2dy = aux.rebin2d(h2dy, metBinning, emhtBinning)
    h2s1 = aux.rebin2d(h2s1, metBinning, emhtBinning)
    h2s2 = aux.rebin2d(h2s2, metBinning, emhtBinning)

    # subtract mc from numerator (small influence on denominator)
    num = h2data.Clone(aux.randomName())
    num.Add(h2ttg, -1)
    num.Add(h2tg, -1)
    num.Add(h2wg, -1)
    num.Add(h2zg, -1)
    num.Add(h2dy, -1)
    num.Add(h2eControl, -0.0197)

    if simpleScale:
        h2QcdW, h2QcdWsys = qcdPrediction2dsimple(num, h2jetControl, 100)
    else:
        h2QcdW, h2QcdWsys = weighting2d(num, h2jetControl, 100)

    if allSets == data:
        h2eControl.Scale(0.0197) # fakeRate Data
    else:
        h2eControl.Scale(0.0107) # fakeRate MC
    for h in h2ttg, h2tg, h2wg, h2zg, h2dy, h2s1, h2s2:
        h.Scale(0.983) # trigger efficiency, uncertainty: 0.002

    cuts = [
        ("x", 0, 1e6),
        ("x", 0, 800),
        ("x", 0, 1000),
        ("x", 800, 1000),
        ("x", 0, 2000),
        ("x", 0, 1500),
        ("x", 700, 1e6),
        ("x", 800, 1e6),
        ("x", 900, 1e6),
        ("x", 1000, 1e6),
        ("x", 1100, 1e6),
        ("x", 1200, 1e6),
        ("x", 1300, 1e6),
        ("x", 1400, 1e6),
        ("x", 1500, 1e6),
        ("x", 2000, 1e6),
        ("y", 0, 1e6),
        ("y", 0, 100),
        ("y", 100, 1e6),
        ]
    #if emhtBinning: cuts = [ ("x", low, high) for low, high in zip(emhtBinning,emhtBinning[1:]+[1e6]) ]

    if axVars["y"] == "emht" and False:
        cuts = [("x", 0, 1e6), ("x", 0, 800), ("x", 0, 1000), ("x", 800, 1000), ("x", 1000, 1e6), ("x", 2000, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "mht":
        cuts = [("x", 0, 10), ("x", 40, 100), ("x", 100, 200), ("x", 200, 300), ("x", 300, 400), ("x", 400, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "gPt":
        cuts = [("x", 100, 300), ("x", 300, 500), ("x", 500, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "jetPt":
        cuts = [("x", 100, 400), ("x", 400, 600), ("x", 600, 1e6), ("y", 0, 1e6)]
    elif axVars["y"] == "njet":
        cuts = [("x", 0, 2), ("x", 3, 4), ("x", 5, 6), ("x", 7, 1e6), ("y", 0, 1e6)]

    for dir, cut1, cut2 in cuts:
        parDir = "xy".replace(dir, "") # get the other axis
        cut1Bin = aux.getAxis(h2data,parDir).FindBin(cut1)
        cut2Bin = aux.getAxis(h2data,parDir).FindBin(cut2-1e-6)
        h1data = aux.getProjection(h2data, dir, cut1Bin, cut2Bin)
        h1QcdW = aux.getProjection(h2QcdW, dir, cut1Bin, cut2Bin)
        h1QcdWsys = aux.getProjection(h2QcdWsys, dir, cut1Bin, cut2Bin)
        h1e = aux.getProjection(h2eControl, dir, cut1Bin, cut2Bin)
        h1ttg = aux.getProjection(h2ttg, dir, cut1Bin, cut2Bin)
        h1tg = aux.getProjection(h2tg, dir, cut1Bin, cut2Bin)
        h1wg = aux.getProjection(h2wg, dir, cut1Bin, cut2Bin)
        h1zg = aux.getProjection(h2zg, dir, cut1Bin, cut2Bin)
        h1dy = aux.getProjection(h2dy, dir, cut1Bin, cut2Bin)
        h1s1 = aux.getProjection(h2s1, dir, cut1Bin, cut2Bin)
        h1s2 = aux.getProjection(h2s2, dir, cut1Bin, cut2Bin)

        for h in h1data, h1QcdW, h1QcdWsys, h1e, h1ttg, h1tg, h1wg, h1zg, h1dy, h1s1, h1s2:
            h.SetYTitle("Events/Bin")
            aux.appendFlowBin(h)
            if style.divideByBinWidth: h.Scale(1, "width")
            if dir == "y": h.SetTitleOffset(1)
        h1QcdW_integral = h1QcdW.Integral()
        if h1QcdW_integral:
            h1QcdW.Scale(h1data.Integral()/h1QcdW_integral)
            h1QcdWsys.Scale(h1data.Integral()/h1QcdW_integral)

        if allSets == data:
            if dir=="y": maxi=2000e5
            if dir=="x": maxi=160
            for bin in range(h1data.FindBin(maxi),h1data.GetNbinsX()+2):
                h1data.SetBinContent(bin,0)

        h1eSys = aux.getSysHisto(h1e,.3)
        h1ttgSys = aux.getSysHisto(h1ttg,.3)
        h1tgSys = aux.getSysHisto(h1tg,.3)
        h1wgSys = aux.getSysHisto(h1wg,.3)
        h1zgSys = aux.getSysHisto(h1zg,.5)
        h1dySys = aux.getSysHisto(h1dy,.5)

        aux.drawOpt(h1data, "data")
        h1QcdW.SetLineColor(rwth.blue75)
        h1e.SetLineColor(rwth.green)
        h1ttg.SetLineColor(rwth.bordeaux)
        h1tg.SetLineColor(rwth.bordeaux75)
        h1wg.SetLineColor(rwth.red75)
        h1zg.SetLineColor(rwth.yellow)
        h1dy.SetLineColor(rwth.yellow50)
        h1s1.SetLineWidth(3)
        h1s1.drawOption_ = "hist"
        h1s2.SetLineWidth(3)
        h1s2.drawOption_ = "hist"
        h1s2.SetLineStyle(2)

        dataLabel = "(Pseudo)Data"
        if allSets == data: dataLabel = "Data"
        if allSets == gjets+qcd: dataLabel = "#gammaJet+QCD"
        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        m.add(h1data, dataLabel)
        if allSets == data:
            m.addStack(h1tg, "t#gamma")
            m.addStack(h1dy, "DY")
            m.addStack(h1e, "e#rightarrow#gamma")
            m.addStack(h1ttg, "t#bar{t}#gamma")
            m.addStack(h1zg, "Z#gamma")
            m.addStack(h1wg, "W#gamma")
            m.add(h1s1, signal["T5Wg_1550_100"].label)
            m.add(h1s2, signal["T5Wg_1550_1500"].label)
        m.addStack(h1QcdW, "#gammaJet+QCD pred.")

        # systematics
        sysStack = ROOT.THStack()
        if allSets == data:
            systHists = h1QcdWsys, h1eSys, h1ttgSys, h1wgSys, h1zgSys, h1dy
        else:
            systHists = [h1QcdWsys]
        for h in systHists:
            sysStack.Add(h)
        sysHist = sysStack.GetStack().Last()
        aux.drawOpt(sysHist, "sys")

        info = aux.getAxis(h2data, parDir).GetTitle().replace(" (GeV)","")
        if cut1: info = str(cut1)+"<"+info
        if cut2<1e5: info = info+"<"+str(cut2)
        if "<" not in info and ">" not in info: info=""
        m.leg.SetHeader(info)
        m.Draw()
        sysHist.Draw("same e2")

        l = aux.Label(sim=allSets is not data)

        r = ratio.Ratio("Data/SM", h1data, m.hists[0].GetStack().Last(), sysHist)
        r.draw(0.5,1.5)

        appendix = ""
        if allSets == data: appendix = "_data"
        if any([ x.startswith("T5") for x in allSets.names ]): appendix = "_sigCont"

        title = "finalPlot" if allSets == data else "closure"
        simpleString = "_simple" if simpleScale else ""
        aux.save("{}{}{}{}_{}_{}{}{}".format(title,appendix,simpleString,axVars[dir],suffix,int(cut1),axVars[parDir],int(cut2)), normal=False)

def getAxis(hname):
    m = re.match("(.*)_vs_(.*)", hname)
    if m:
        return m.groups()
    else:
        print "Could not determine x and y variable from", hname


def projectAll2d(sampleName, hname, dir1, dir2, dataset1, dataset2=None):
    if not dataset2: dataset2 = dataset1

    varX, varY = getAxis(hname)
    dim1Binning = aux.getBinningsFromName(varX)["3"]
    dim2Binning = aux.getBinningsFromName("emht")["2"]
    #dim2Binning = range(0,1000, 10)
    dim2Binning = [-5000,5000]
    #dim2Binning = None

    h2Raw1 = dataset1.getHist(dir1+"/"+hname)
    h2Raw2 = dataset2.getHist(dir2+"/"+hname)

    h2_1 = aux.rebin2d(h2Raw1, dim1Binning, dim2Binning)
    h2_2 = aux.rebin2d(h2Raw2, dim1Binning, dim2Binning)

    h2_2_w, h2_2_wSys = weighting2d(h2_1, h2_2)
    #h2_2_w, h2_2_wSys = h2_2, h2_2
    #h2_2_w.Scale(h2_1.Integral()/h2_2.Integral())
    #h2_2_wSys.Scale(h2_1.Integral()/h2_2.Integral())

    cuts = [
        ("xy", 0, 1e6),
        ("yx", 0, 100),
        ("yx", 100, 1e6),
        ]
    for axOrder, cut1, cut2 in cuts:
        ax, parAx = axOrder
        cut1Bin = aux.getAxis(h2_1, parAx).FindBin(cut1)
        cut2Bin = aux.getAxis(h2_1, parAx).FindBin(cut2-1e-6)
        h1_1 = aux.getProjection(h2_1, ax, cut1Bin, cut2Bin)
        h1_2w = aux.getProjection(h2_2_w, ax, cut1Bin, cut2Bin)
        h1_2wSys = aux.getProjection(h2_2_wSys, ax, cut1Bin, cut2Bin)

        for h in h1_1, h1_2w, h1_2wSys:
            h.Scale(1, "width" if style.divideByBinWidth else "")
            h.SetYTitle("Events/Bin") # TODO: fix if divide by bin
            aux.appendFlowBin( h )

        aux.drawOpt(h1_1, "data")
        aux.drawOpt(h1_2w, "pre")
        aux.drawOpt(h1_2wSys, "sys")

        c = ROOT.TCanvas()
        m = multiplot.Multiplot()
        m.add(h1_1, folder2tex[dir1])
        m.add(h1_2w, folder2tex[dir2])
        m.add(h1_2wSys, "#sigma_{weight}")

        info = ""
        if cut1: info += "{}<".format(cut1)
        if cut1 or cut2 < 1e6: info += aux.getAxis(h2_1, parAx).GetTitle().replace(" (GeV)","")
        if cut2>0 and cut2 < 1e6: info += "<{}".format(cut2)
        if info: m.leg.SetHeader(info)
        m.Draw()

        l = aux.Label(sim=dataset1 is not data, info=dataset1.label)

        r = ratio.Ratio("{}/{}".format(folder2tex[dir1],folder2tex[dir2]), h1_1, h1_2w, h1_2wSys )
        r.draw(0.5,1.5)

        funcName = sys._getframe().f_code.co_name
        aux.save("{}_{}_{}_{}_{}_{}to{:.0f}".format(funcName, sampleName, hname, axOrder[0], dir1, cut1, cut2), normal=False)

def getStyledHists(dataset1, name, binning, dirs):
    out = None
    for d in dirs:
        h = getStyledHist(dataset1, name, binning, d)
        if out:
            out.Add(h)
        else:
            out = h
    return out

def beautifyDir(s):
    m = re.match("(.*)_(\d)_(\d)", s)
    label, pos, nPart = m.groups()
    out = ["j"]*int(nPart)
    obj = "j" if "jControl" in label else "#gamma"
    out[int(pos)] = "#font[22]{%s}"%obj
    return ''.join(out)


def compareAddPositions(sampleName, dataset1, selections, dataset2=None):
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange]
    if not dataset2: dataset2 = dataset1
    for name, binningName, binning in getBinningAndNames():
        can = ROOT.TCanvas()
        m = multiplot.Multiplot()
        normInt = 0
        hists = []
        for iS, s in enumerate(selections):
            if "jControl" in s:
                h = getStyledHist(dataset2, name, binning, s)
            else:
                h = getStyledHist(dataset1, name, binning, s)
            thisInt = h.Integral()
            if thisInt and not normInt: normInt = thisInt
            if thisInt: h.Scale(normInt/thisInt)
            h.SetLineColor(colors[iS])
            m.add(h, beautifyDir(s))
            hists.append(h)
        additionals = {
#            "met": ("metUp5","metDn5"),
#            "metPar": ("metParUp5", "metParDn5"),
        }
        if name in additionals:
            additionalHists = []
            for varHist in additionals[name]:
                hUp = getStyledHist(dataset2, varHist, binning, s)
                thisInt = hUp.Integral()
                if thisInt: hUp.Scale(normInt/thisInt)
                hUp.SetLineColor(ROOT.kOrange)
                hUp.SetMarkerColor(ROOT.kOrange)
                m.add(hUp, "")
                additionalHists.append(hUp)
#        else: continue
        syst = aux.getSystFromDifference(additionalHists[0], additionalHists[1]) if name in additionals else None

        m.Draw()
        l = aux.Label(info=dataset1.label, sim=dataset1 is not data)
        if len(selections)>1:
            r = ratio.Ratio("#frac{{{}}}{{{}}}".format(selections[0], selections[1]), hists[0], hists[1], syst)
            r.draw(.5, 1.5)

        funcName = sys._getframe().f_code.co_name
        if binningName: binningName = "_" + binningName
        dirs = "_vs_".join(selections)
        saveName = "{}_{}_{}_{}{}".format(funcName, sampleName, dirs, name, binningName)
        aux.save(saveName)

def getMetParInfoForDifferentSelections(dataset):
    dirs = [ x for x in getNames(dataset,"", objects=[ROOT.TDirectory]) if x.startswith("tr_")]
    for d in dirs:
        for hname in "metParRaw", "metPar":
            h = dataset.getHist("{}/{}".format(d,hname))
            print d, hname, h.GetMean(), h.GetRMS()


if __name__ == "__main__":
    #compareDirs("gqcd", "tr", "tr_jControl", gjets+qcd)
    """
    projectAll2d("gqcd", "met_vs_emht", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_150", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_200", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_250", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_300", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_350", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_400", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_njet3", "tr", "tr_jControl", gjets+qcd)
    projectAll2d("gqcd", "metRaw_vs_emht_njet5", "tr", "tr_jControl", gjets+qcd)
    """
    """
    projectAll2d("data", "metRaw_vs_emht", "tr", "tr_jControl", data)
    projectAll2d("data", "metRaw_vs_emht_110", "tr", "tr_jControl", data)
    projectAll2d("data", "metRaw_vs_emht_120", "tr", "tr_jControl", data)
    projectAll2d("data", "metRaw_vs_emht_130", "tr", "tr_jControl", data)
    """


    compareAddPositions("gqcd", gjets+qcd, ["tr_0_1","tr_0_2", "tr_0_3", "tr_0_4"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_0_2","tr_1_2"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_0_3","tr_1_3", "tr_2_3"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_0_4","tr_1_4", "tr_2_4", "tr_3_4"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_0_1","tr_jControl_0_1"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_0_2","tr_jControl_0_2"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_1_2","tr_jControl_1_2"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_0_3","tr_jControl_0_3"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_1_3","tr_jControl_1_3"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_2_3","tr_jControl_2_3"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_0_4","tr_jControl_0_4"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_1_4","tr_jControl_1_4"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_2_4","tr_jControl_2_4"])
    compareAddPositions("gqcd", gjets+qcd, ["tr_3_4","tr_jControl_3_4"])


    #finalPrediction(gjets+qcd, suffix="_3_4")

    """
    compareAddPositions("data", data, ["tr_0_1","tr_jControl_0_1"], dataHt)
    compareAddPositions("data", data, ["tr_0_2","tr_jControl_0_2"], dataHt)
    compareAddPositions("data", data, ["tr_1_2","tr_jControl_1_2"], dataHt)
    compareAddPositions("data", data, ["tr_0_3","tr_jControl_0_3"], dataHt)
    compareAddPositions("data", data, ["tr_1_3","tr_jControl_1_3"], dataHt)
    compareAddPositions("data", data, ["tr_2_3","tr_jControl_2_3"], dataHt)
    compareAddPositions("data", data, ["tr_0_4","tr_jControl_0_4"], dataHt)
    compareAddPositions("data", data, ["tr_1_4","tr_jControl_1_4"], dataHt)
    compareAddPositions("data", data, ["tr_2_4","tr_jControl_2_4"], dataHt)
    compareAddPositions("data", data, ["tr_3_4","tr_jControl_3_4"], dataHt)
    """
