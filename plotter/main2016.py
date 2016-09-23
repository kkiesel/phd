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

def getMetMeanRMS(h, mini=0,maxi=100):
    h.GetXaxis().SetRangeUser(mini,maxi)
    m, mE = h.GetMean(), h.GetMeanError()
    rms, rmsE = h.GetRMS(), h.GetRMSError()
    h.GetXaxis().SetRange()
    return m, mE, rms, rmsE

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
    #hRes = h.Clone()
    #for bin in range(nBins):
    #    hRes.SetBinContent(bin+1, means[bin]-medians[bin])
    #hRes.Draw("hist")
    #ROOT.gPad.SaveAs("test.pdf")
    return meanEnd, cutEnd


def getGJetPrediction(dirTree, preTree, variable, weight, nBins):

    dirHist = createHistoFromTree(dirTree, variable, weight, nBins, 0, 200)
    preHistUnw = createHistoFromTree(preTree, variable, weight, nBins, 0, 200)

    m1, cut1 = getMetMean(dirTree, variable, weight)
    m2, cut2 = getMetMean(preTree, variable, weight)
    m11, m1E, rms1, rms1E = getMetMeanRMS(dirHist, 0, cut1)
    m22, m2E, rms2, rms2E = getMetMeanRMS(preHistUnw, 0, cut2)

    rmsDiff = aux.sqrt( (m1-m2)**2 + (rms1-rms2)**2 + m1E**2 + m2E**2 + rms1E**2 +rms2E**2 )
    corrDn, corr, corrUp = m1/(m2-rmsDiff), m1/m2, m1/(m2+rmsDiff)

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

def createHistoFromDatasetTree(dset, variable, weight, nBins):
    tree = ROOT.TChain("tr/simpleTree")
    for f in dset.files: tree.Add(f)
    h = createHistoFromTree(tree, variable, weight, nBins)
    h.SetLineColor(dset.color)
    return h

def fitPrediction(dirTree, preTree, variable, weight, isData=True):
    binning = range(0,100,1)
    scales = [.85+.01*i for i in range(20)]
    gr = ROOT.TGraph(len(scales))
    dirHist = createHistoFromTree(dirTree, variable, weight, binning)
    for iscale, scale in enumerate(scales):
        preHist = createHistoFromTree(preTree, "{}*{}".format(variable,scale), weight, binning)
        preHist.Scale(dirHist.Integral(0,100)/preHist.Integral(0,100))
        chi2 = dirHist.Chi2Test(preHist, "OF UF CHI2 WW NORM")
        c = ROOT.TCanvas()
        dirHist.Draw()
        preHist.SetLineColor(2)
        preHist.Draw("same")
        r = ratio.Ratio("g/p", dirHist,preHist)
        r.draw(.5,1.5)
        l = aux.Label(info="scale={}  chi2={}".format(scale, chi2))
        aux.save("fitPrediction_scale{}percent".format(scale*100), log=False)
        gr.SetPoint(iscale, scale, chi2)
        print scale, chi2
    c = ROOT.TCanvas()
    gr.Draw()
    gr.Fit("pol2", "W")
    aux.save("fitPredictionTest", log=False)

def qcdClosure(name, dirSet, preSet=None, additionalSets=[], cut="1"):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain("tr/simpleTree")
    preTree = ROOT.TChain("tr_jControl/simpleTree")

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,250,10)+range(250,500,20)+[600,700,800,900,910]
    #fitPrediction(dirTree, preTree, variable, weight, isData=False)
    #return

    dirHist = createHistoFromTree(dirTree, variable, weight, nBins)
    preHist, gjetSyst = getGJetPrediction(dirTree, preTree, variable, weight, nBins)
    totUnc = aux.addHistUncert(preHist, gjetSyst)

    # beautify
    aux.drawOpt(dirHist, "data")
    aux.drawOpt(preHist, "pre")
    aux.drawOpt(totUnc, "totUnc")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add(dirHist, "(Pseudo)Data")
    m.add(preHist, "")
    m.add(totUnc, "Tot. uncert.")
    if cut is not "1": m.leg.SetHeader(cut)
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, preHist, gjetSyst)
    r.draw(0.5, 1.5)

    l = aux.Label(sim= not dirSet==data, info=dirSet.label)
    aux.save("qcdClosure_{}".format(name))

def finalDistribution(name, dirSet, preSet=None, additionalSets=[], cut="1"):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain("tr/simpleTree")
    preTree = ROOT.TChain("tr_jControl/simpleTree")
    eTree = ROOT.TChain("tr_eControl/simpleTree")

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,250,10)+range(250,500,20)+[600,700,800,900,910]

    dirHist = createHistoFromTree(dirTree, variable, weight, nBins)
    gjetHist, gjetSyst = getGJetPrediction(dirTree, preTree, variable, weight, nBins)
    gjetHist = createHistoFromTree(preTree, variable, weight, nBins)
    gjetHist.Scale(dirHist.Integral(0,100)/gjetHist.Integral(0,100))
    gjetHist.SetLineColor(ROOT.kCyan-1)

    eHist = createHistoFromTree(dirTree, variable, weight, nBins)
    eHist.Scale( 0.0197 if dirSet==data else 0.0107 )
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    zgHist = createHistoFromDatasetTree(zg, variable, weight, nBins)
    wgHist = createHistoFromDatasetTree(wg_130, variable, weight, nBins)
    ttgHist = createHistoFromDatasetTree(ttg, variable, weight, nBins)
    zgSyst = aux.getSysHisto(zgHist, 0.3)
    wgSyst = aux.getSysHisto(wgHist, 0.3)
    ttgSyst = aux.getSysHisto(ttgHist, 0.3)

    totStat = aux.addHists(gjetHist, eHist, zgHist, wgHist, ttgHist)
    totSyst = aux.addHists(gjetSyst, eSyst, zgSyst, wgSyst, ttgSyst)
    totUnc = aux.addHistUncert(totStat, totSyst)

    # beautify
    aux.drawOpt(dirHist, "data")
    aux.drawOpt(totUnc, "totUnc")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add(dirHist, "(Pseudo)Data")
    m.addStack(eHist, "e#rightarrow#gamma")
    m.addStack(ttgHist, "t#bar{t}#gamma")
    m.addStack(zgHist, "Z#gamma")
    m.addStack(wgHist, "W#gamma")
    m.addStack(gjetHist, "#gamma + Jet")
    m.add(totUnc, "Tot. uncert.")
    if cut is not "1": m.leg.SetHeader(cut)
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, totStat, totSyst)
    r.draw(0.5, 1.5)

    l = aux.Label(sim= not dirSet==data, info=dirSet.label)
    aux.save("finalDistribution_{}".format(name))




if __name__ == "__main__":
    finalDistribution("data", data, dataHt)
    #qcdClosure("data", data, dataHt)
    #qcdClosure("gqcd", gjets+qcd)
    qcdClosure("gqcd_emht1000", gjets+qcd, cut="emht<1000")
    qcdClosure("gqcd_2000emht", gjets+qcd, cut="emht>2000")

