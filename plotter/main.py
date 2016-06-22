#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
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
        aux.drawOpt(h1, "data")
        aux.drawOpt(h2, "pre")
        h2.Scale(i1/i2)

        h2.drawOption_="hist e"

        can = ROOT.TCanvas()
        m = multiplot.Multiplot()
        m.add(h1, "{} {}".format(folder2tex[dir1], aux.metricPrefix(i1)))
        m.add(h2, "{} {}".format(folder2tex[dir2], aux.metricPrefix(i2)))
        m.Draw()
        l = aux.Label(info=dataset1.label, sim=dataset1 is not data)
        r = ratio.Ratio("{}/{}".format(folder2tex[dir1], folder2tex[dir2]), h1, h2)
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

def projectAll2d(sampleName, hname, dir1, dir2, dataset1, dataset2=None):
    if not dataset2: dataset2 = dataset1

    dim1Binning = aux.getBinningsFromName("met")["3"]
    dim2Binning = aux.getBinningsFromName("emht")["2"]

    h2Raw1 = dataset1.getHist(dir1+"/"+hname)
    h2Raw2 = dataset2.getHist(dir2+"/"+hname)

    h2_1 = aux.rebin2d(h2Raw1, dim1Binning, dim2Binning)
    h2_2 = aux.rebin2d(h2Raw2, dim1Binning, dim2Binning)

    h2_2_w, h2_2_wSys = weighting2d(h2_1, h2_2)

    cuts = [
        ("xy", 0, 1e6),
        ("xy", 0, 2000),
        ("xy", 0, 800),
        ("xy", 0, 1000),
        ("xy", 1000, 2000),
        ("xy", 2000, 1e6),
        ("xy", 700, 1e6),
        ("xy", 800, 1e6),
        ("xy", 900, 1e6),
        ("xy", 1000, 1e6),
        ("xy", 1100, 1e6),
        ("xy", 1200, 1e6),
        ("xy", 1300, 1e6),
        ("xy", 1400, 1e6),
        ("xy", 1500, 1e6),
        ("xy", 2000, 1e6),
        ("yx", 0, 1e6),
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
        if cut1 or cut2: info += aux.getAxis(h2_1, parAx).GetTitle().replace(" (GeV)","")
        if cut2>0: info += "<{}".format(cut2)
        m.leg.SetHeader(info)
        m.Draw()

        l = aux.Label(sim=True, info=dataset1.label)

        r = ratio.Ratio("{}/{}".format(folder2tex[dir1],folder2tex[dir2]), h1_1, h1_2w, h1_2wSys )
        r.draw(0.5,1.5)

        funcName = sys._getframe().f_code.co_name
        aux.save("{}_{}_{}_{}_{}to{:.0f}".format(funcName, sampleName, hname, dir1, cut1, cut2))

def allSame():
    drawSameHistogram("mc", "tr", [gjets, qcd, ttjets, ttg, wjets, wg_mg, zg_130, znunu], [signal[s] for s in "T5Wg_1550_100","T5Wg_1550_1500"])
    drawSameHistogram("mc", "tr_ee", [gjets, qcd, ttjets, ttg, wjets, wg_mg, zg_130, znunu], [signal[s] for s in "T5Wg_1550_100","T5Wg_1550_1500"])

def allCompare():
    compareDirs("gqcd", "tr", "tr_ee", gjets+qcd)
    compareDirs("gqcd", "tr_tight", "tr_ee", gjets+qcd)
    compareDirs("data", "tr_tight", "tr_ee", data)
    compareDirs("ttg", "tr_tight", "tr_ee", ttg)
    compareDirs("wg", "tr_tight", "tr_ee", wg_mg)
    compareDirs("zg", "tr_tight", "tr_ee", zg_130)

def allWeighed():
    projectAll2d("gqcd", "met_vs_emht", "tr_tight", "tr_ee", gjets+qcd)

if __name__ == "__main__":
    allSame()
    allCompare()
