#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

selectionLabels = {
    "tr_ee": "EE",
    "tr_jControl_ee": "EE",
    "tr_eControl_ee": "EE"
}

def sampleCompositionFromTree(name, samples, dirName, variable, cut, binning):
    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    for s in samples:
        h = aux.createHistoFromDatasetTree(s, variable, "weight*({})".format(cut), binning, dirName+"/simpleTree")
        m.addStack(h, s.label)
    m.sortStackByIntegral()
    m.Draw()
    ratio.clearXaxisCurrentPad()
    p = ratio.createBottomPad()
    x = aux.drawContributions(m.getStack())
    l = aux.Label(sim=True)
    saveName = "sameHistogramsFromTree_{}_{}_{}_{}".format(name, dirName, variable, len(binning))
    aux.save( saveName )

def drawSameHistogram(sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scale=False):
    if "jControl" in name and data in additional:
        additional.remove(data)
        additional.append(dataHt)
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    dataIntegral = 0
    for d in additional:
        if d == data or d == dataHt:
            h = d.getHist(name)
            dataIntegral = h.Integral(0,-1)
    bkgHist = None
    for d in bkg:
        h = d.getHist(name)
        if bkgHist: bkgHist.Add(h)
        else: bkgHist = h
    bkgIntegral = bkgHist.Integral(0,-1)

    for d in bkg[-1::-1]:
        h = d.getHist(name)
        if not h: continue
        if not h.Integral(): continue
        if scale: h.Scale(dataIntegral/bkgIntegral)
        if binning: h = aux.rebin( h, binning )

        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        m.addStack( h, d.label )

    dataHist = None
    for d in additional:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        if binning: h = aux.rebin( h, binning )
        aux.appendFlowBin( h )
        if d == data and name == "tr/met":
            for bin in range(h.GetNbinsX()+2):
                if h.GetBinCenter(bin)>200:
                    h.SetBinContent(bin,0);h.SetBinError(bin,0)

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            # disable errors for data, so that ErrorOption is working
            # TODO: kPoisson also for rebinned and scaled histograms
            if not binning: h.Sumw2(False)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
            dataHist = h
        else:
            h.drawOption_ = "hist"
            h.SetLineWidth(3)

        m.add( h, d.label )

    m.sortStackByIntegral()
    if name.split("/")[0] in selectionLabels:
        m.leg.SetHeader(selectionLabels[name.split("/")[0]])
    if not m.Draw(): return

    # ratio
    hsm = m.hists[0].GetStack().Last()
    if dataHist:
        r = ratio.Ratio( "Data/SM", dataHist, hsm )
        rMean = dataHist.Integral()/hsm.Integral()
        if rMean > 2 or rMean < 0.25: # for jcontrol with prescaled data
            r.draw(rMean/2,1.5*rMean)
        else:
            x = r.draw(0,1.5, m.getStack())
    else:
        ratio.clearXaxisCurrentPad()
        p = ratio.createBottomPad()
        x = aux.drawContributions(m.getStack())


    l = aux.Label(sim=data not in additional)

    if binningName: binningName = "_"+binningName
    name = name.replace("/","__")
    saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
    aux.save( saveName )

    if "emht" in name and dataHist:
        aux.write2File(r.ratio, saveName, "emhtWeights.root")


def drawSameHistograms( sampleNames="test", stack=[], additional=[] ):
    file = stack[0].files[0] if stack else additional[0].files[0]
    names = aux.getObjectNames( file, "tr", [ROOT.TH1F] )
    dirs = [d for d in aux.getDirNames(file) if "tr_" in d and not "gen" in d and not "true" in d]

    if data in additional:
        if "genMatch" in names: names.remove("genMatch")
        if "genHt" in names: names.remove("genHt")

    names = ["met", "emht"]
    dirs = "tr", "tr_ee", "tr_jControl", "tr_eControl","tr_eControl_ee", "tr_0met100"
    #dirs = ["tr_jControl", "tr_jControl_noLep"]
    #dirs = ["tr_genWZ11","tr_genWZ11",  "tr_jControl", "tr_eControl", "tr"]
    #names = ["g_eta"]
    #dirs = ["tr", "tr_ee"]

    for name in names:
        for binningName, binning in aux.getBinningsFromName( name ).iteritems():
            #if binningName != "1": continue
            for directory in dirs:
                if "jControl" in directory and data in additional:
                    drawSameHistogram( sampleNames, directory+"/"+name, stack, [dataHt]+[x for x in additional if x != data], binning, binningName )
                else:
                    drawSameHistogram( sampleNames, directory+"/"+name, stack, additional, binning, binningName )


if __name__ == "__main__":
    #drawSameHistograms( "gqcd_data", [gjets, qcd], additional=[data])
    #drawSameHistograms( "mc_data", [gjets, qcd, ttjets_ht, ttg, wjets, wg, zg, znunu], additional=[data])
    #drawSameHistograms( "mc", [gjets, qcd, ttjets, ttg, wjets, wg, zg, znunu], additional=[t5wg_1600_100, t5wg_1600_1500])
    #drawSameHistograms( "gqcd_splitted", [gjets200, gjets400, gjets600, qcd500, qcd700, qcd1000, qcd1500, qcd2000])
    #sampleCompositionFromTree("gqcd_2000", [gjets200, gjets400, gjets600, qcd500, qcd700, qcd1000, qcd1500, qcd2000], "tr", "met", "2000<emht", range(0,200,10)+[200, 300, 400, 500, 600])
    #sampleCompositionFromTree("gqcd_ee_2000", [gjets200, gjets400, gjets600, qcd500, qcd700, qcd1000, qcd1500, qcd2000], "tr_ee", "met", "2000<emht", range(0,200,10)+[200, 300, 400, 500, 600])
    #sampleCompositionFromTree("gqcd_2000", [gjets200, gjets400, gjets600, qcd500, qcd700, qcd1000, qcd1500, qcd2000], "tr_jControl", "met", "2000<emht", range(0,200,10)+[200, 300, 400, 500, 600])
    #sampleCompositionFromTree("gqcd_ee_2000", [gjets200, gjets400, gjets600, qcd500, qcd700, qcd1000, qcd1500, qcd2000], "tr_ee", "met", "2000<emht", range(0,200,10)+[200, 300, 400, 500, 600])
    #sampleCompositionFromTree("gqcd", [gjets200, gjets400, gjets600, qcd500, qcd700, qcd1000, qcd1500, qcd2000], "tr", "met", "1", range(0,200,10)+[200, 300, 400, 500, 600])
    #sampleCompositionFromTree("gqcd_jControl", [gjets200dr, gjets400dr, gjets600dr, qcd500, qcd700, qcd1000, qcd1500, qcd2000], "tr_jControl", "met", "emht>2000", range(0,200,10)+[200, 250, 300, 350, 450, 600, 700])
    #sampleCompositionFromTree("gqcd_jControl", [gjets200dr, gjets400dr, gjets600dr, qcd500, qcd700, qcd1000, qcd1500, qcd2000], "tr_jControl", "met", "emht>2000", range(0,200,10)+[200, 250, 300, 350, 450, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500])

    nBins = range(0,200,10)+[200, 250, 300, 350, 450, 600, 700]
    style.divideByBinWidth = True
    #drawSameHistogram("sampleComposition_qcdClosure_highEMHT", "signal_highEMHT/met", [qcd500, qcd700, qcd1000, qcd1500, qcd2000,gjets40dr, gjets100dr, gjets200dr, gjets400dr, gjets600dr], [], nBins, "")
    #drawSameHistogram("sampleComposition_electronClosure_highEMHT", "signal_highEMHT/met", [wjets200, wjets400, wjets600, wjets800, wjets1200, wjets2500, ttjets0, ttjets600, ttjets800, ttjets1200, ttjets2500], [], nBins, "")
    #drawSameHistogram("sampleComposition_electronClosure_highEMHT", "signal_highEMHT/met", [wjets1200, wjets2500, ttjets0, ttjets600, ttjets800, ttjets1200, ttjets2500], [], nBins, "2")
    #drawSameHistogram("sampleComposition_electronClosure_highEMHT", "signal_highEMHT/met", [wjets1200, wjets2500, ttjets600, ttjets800, ttjets1200, ttjets2500], [], nBins, "3")
    #drawSameHistogram("sampleComposition_electronClosure_highEMHT", "signal_highEMHT/met", [wjets1200, wjets2500, ttjets800, ttjets1200, ttjets2500], [], nBins, "4")
    #drawSameHistogram("qcd_wztjets", "tr_jControl/met", [qcd, wjets, ttjets_ht, znunu], [], nBins, "")
    style.divideByBinWidth = False


