#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def createHistoFromTree(tree, variable, weight="", nBins=20, firstBin=None, lastBin=None ):
    """
    tree: tree to create histo from
    variable: variable to plot (must be a branch of the tree)
    weight: weights to apply (e.g. "var1*(var2 > 15)" will use weights from var1 and cut on var2 > 15
    nBins, firstBin, lastBin: number of bins, first bin and last bin (same as in TH1F constructor)
    nBins: if nBins is a list, and to a int, a user binned plot will be generated
    returns: histogram
    """
    from ROOT import TH1F
    name = aux.randomName()
    if isinstance( nBins, list ):
        import array
        xBins = array.array('d', nBins )
        result = TH1F(name, variable, len(nBins)-1, xBins)
        result.Sumw2()
        tree.Draw("%s>>%s"%(variable, name), weight, "goff")
    elif firstBin==None and lastBin==None:
        import ROOT
        tree.Draw("%s>>%s(%s,,)"%(variable,name,nBins), weight, "goff")
        result = ROOT.gDirectory.Get( name )
        if isinstance( result, ROOT.TTree ):
            print "Warning, no entries"
            return ROOT.TH1F()
        result.Sumw2() # applying the errors here is perhaps not entirely correct
    else:
        result = TH1F(name, variable, nBins, firstBin, lastBin)
        result.Sumw2()
        tree.Draw("%s>>%s"%(variable, name), weight, "goff")
    aux.appendFlowBin(result)
    result.SetTitle("")
    if variable.startswith("met"):
        result.SetXTitle("#it{E}_{T}^{miss} (GeV)")
    aux.setYAxisTitle(result)
    return result


def getMetMean(tree, variable, weight):
    import bisect
    nBins = 10000
    h = createHistoFromTree(tree, variable, weight, nBins, 0, 100)
    cum, cumW, means, medians = [], [], [], []
    meanEnd, cutEnd, tmp = 0, 0, 0
    for bin in range(nBins):
        x = h.GetBinCenter(bin+1)
        c = h.GetBinContent(bin+1)
        oldC = cum[-1] if cum else 0
        oldX = cumW[-1] if cumW else 0
        cum.append( oldC+c )
        cumW.append( oldX+x*c )
        mean = cumW[-1]/cum[-1] if cum[-1] else 0
        means.append(mean)
        median = h.GetBinLowEdge(bisect.bisect(cum, cum[-1]/2))
        medians.append(median)
        tmp2 = mean-median
        if tmp2 * tmp < 0:
            meanEnd = mean
            cutEnd = x
        tmp = tmp2
    return meanEnd, cutEnd

def getMetMeanRMS(h, mini=0,maxi=100):
    h.GetXaxis().SetRangeUser(mini,maxi)
    m, mE = h.GetMean(), h.GetMeanError()
    rms, rmsE = h.GetRMS(), h.GetRMSError()
    h.GetXaxis().SetRange()
    return m, mE, rms, rmsE


def getGJetPrediction(dirTree, preTree, variable, weight, nBins):

    dirHist = createHistoFromTree(dirTree, variable, weight, nBins, 0, 200)
    preHistUnw = createHistoFromTree(preTree, variable, weight, nBins, 0, 200)

    m1, cut1 = getMetMean(dirTree, variable, weight)
    m2, cut2 = getMetMean(preTree, variable, weight)
    m11, m1E, rms1, rms1E = getMetMeanRMS(dirHist, 0, cut1)
    m22, m2E, rms2, rms2E = getMetMeanRMS(preHistUnw, 0, cut2)
    #m1, m1E, rms1, rms1E = getMetMeanRMS(dirHist, 0, 100)
    #m2, m2E, rms2, rms2E = getMetMeanRMS(preHistUnw, 0, 100)

    rmsDiff = aux.sqrt( (m1-m2)**2 + (rms1-rms2)**2 + m1E**2 + m2E**2 + rms1E**2 +rms2E**2 )
    corrDn, corr, corrUp = m1/(m2-rmsDiff), m1/m2, m1/(m2+rmsDiff)

    print
    print "Using events with MET < {} ({})".format(cut1, cut2)
    print "Dir: μ = {}({}) ± {}  σ = {} ± {}".format(m1, m11, m1E, rms1, rms1E)
    print "Pre: μ = {}({}) ± {}  σ = {} ± {}".format(m2, m22, m2E, rms2, rms2E)
    print "Scales-1:", corrDn-1, corr-1, corrUp-1
    print "Uncertainties:", abs(m1-m2), abs(rms1-rms2), m1E, m2E, rms1E, rms2E

    preHist = createHistoFromTree(preTree, "{}*{}".format(variable,corr), weight, nBins, 0, 200)
    preHistUp = createHistoFromTree(preTree, "{}*{}".format(variable,corrUp), weight, nBins, 0, 200)
    preHistDn = createHistoFromTree(preTree, "{}*{}".format(variable,corrDn), weight, nBins, 0, 200)
    syst = aux.getSystFromDifference(preHistDn, preHistUp)

    # Scale
    dirInt = dirHist.Integral(0, dirHist.FindBin(100))
    preInt = preHist.Integral(0, preHist.FindBin(100))
    for h in preHist, syst:
        h.Scale(dirInt/preInt)

    # Force symmetric systematic uncertainty
    for bin in range(syst.GetNbinsX()+2):
        syst.SetBinContent(bin, preHist.GetBinContent(bin))
    return preHist, syst


if __name__ == "__main__":
    filename = "../toyStudies/toy_tree.root"
    dirTree = ROOT.TChain("gam")
    preTree = ROOT.TChain("jet")
    sTree = ROOT.TChain("gamjet")
    for t in dirTree, preTree, sTree:
        t.AddFile(filename)

    variable = "met"
    weight = "1"
    nBins = range(0,200,5)

    dirHist = createHistoFromTree(dirTree, variable, weight, nBins)
    preHist, gjetSyst = getGJetPrediction(dirTree, preTree, variable, weight, nBins)
    totUnc = aux.addHistUncert(preHist, gjetSyst)

    gjetHistUnw = createHistoFromTree(preTree, variable, weight, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100))/gjetHistUnw.Integral(0,dirHist.FindBin(100)))
    gjetHistUnw.SetLineColor(ROOT.kBlue)

    smearedHist = createHistoFromTree(sTree, variable, weight, nBins)

    # beautify
    aux.drawOpt(dirHist, "data")
    aux.drawOpt(preHist, "pre")
    aux.drawOpt(totUnc, "totUnc")
    aux.drawOpt(smearedHist, "data")
    smearedHist.SetMarkerColor(ROOT.kGray)

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add(dirHist, "Toy Data")
    #m.add(smearedHist, "Smeared toys")
    m.add(preHist, "")
    m.add(totUnc, "Tot. uncert.")
    m.add(gjetHistUnw, "Unweighted")
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, preHist, gjetSyst)
    r.draw(0.5, 1.5)
    gjetHistUnwClone = gjetHistUnw.Clone(aux.randomName())
    gjetHistUnwClone.Divide(preHist)
    gjetHistUnwClone.Draw("same hist")
    smearedHistClone = smearedHist.Clone(aux.randomName())
    smearedHistClone.Divide(preHist)
    #smearedHistClone.Draw("same hist")

    leg = ROOT.TLegend(.17,.12,.39,.18)
    leg.SetFillStyle(0)
    leg.AddEntry(gjetHistUnw, "Unweighted/Weighted", "l")
    leg.Draw()
    info = ROOT.TLatex( .02, .95, "Toys (p_{T},#phi): (100,#frac{#pi}{4}), (200,#frac{#pi}{4}), (200,#frac{3#pi}{2}), (100,#frac{5#pi}{4})")
    info.SetNDC()
    info.Draw()
    aux.save("toyStudies")

