#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins):
    weightName = weight.replace("(","").replace(")","").replace("*","").replace(">","").replace("<","").replace("&&","").replace(" ", "")
    saveName = "__".join([name, str(len(dirSet.names)), treename.replace("/","_"), str(len(preSet.names)), weightName, variable, str(len([x for x in nBins if x<=100]))])
    saveNameRoot = "savedFitPredictions/{}.root".format(saveName)
    if os.path.isfile(saveNameRoot):
        print "Using saved events from", saveName
        f = ROOT.TFile(saveNameRoot)
        dirHist = f.Get("dir")
        preHists = {}
        for key in f.GetListOfKeys():
            name = key.GetName()
            if name == "dir": continue
            preHists[float(name)] = key.ReadObj()
    else:
        print "Calculating new", saveName
        dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
        dirHist.SetName("dir")
        scales = [.85+.01*i for i in range(30)]
        preHists = {}
        for iscale, scale in enumerate(scales):
            preHist = aux.createHistoFromTree(preTree, "{}*{}".format(variable,scale), weight, nBins)
            preHist.SetName(str(scale))
            preHists[scale] = preHist
        # write histograms to file
        f = ROOT.TFile(saveNameRoot, "recreate")
        dirHist.Write()
        for k, h in preHists.iteritems(): h.Write()
        f.Close()

    maxBin = dirHist.FindFixBin(100)-1 # do not take met=100 bin
    dirInt, dirIntErr = aux.integralAndError(dirHist, 0, maxBin)
    for bin in range(maxBin+1, dirHist.GetNbinsX()+2):
        dirHist.SetBinContent(bin, 0)
        dirHist.SetBinError(bin, 0)
    gr = ROOT.TGraph()
    preHists = collections.OrderedDict(sorted(preHists.items(), key=lambda t: t[0]))
    for iscale, (scale, preHist) in enumerate(preHists.iteritems()):
        for bin in range(maxBin+1, dirHist.GetNbinsX()+2):
            preHist.SetBinContent(bin, 0)
            preHist.SetBinError(bin, 0)

        preHist.Scale(dirInt/preHist.Integral(0,maxBin))
        chi2 = dirHist.Chi2Test(preHist, "OF UF CHI2 WW")
        c = ROOT.TCanvas()
        dirHist.GetXaxis().SetRangeUser(0,100)
        preHist.GetXaxis().SetRangeUser(0,100)
        dirHist.Draw()
        preHist.SetLineColor(2)
        preHist.Draw("same")
        r = ratio.Ratio("#gamma/Jet", dirHist,preHist)
        r.draw(.95,1.05)
        l = aux.Label(info="Scale={}  #chi^{{2}}={}".format(scale, chi2))
        #aux.save(saveName+"_scale{}percent".format(int(scale*100)), "savedFitPredictions/",log=False)
        gr.SetPoint(iscale, scale, chi2)
    c = ROOT.TCanvas()
    points = [(i, gr.GetX()[i],gr.GetY()[i]) for i in range(gr.GetN())]
    ys = [gr.GetY()[i] for i in range(gr.GetN())]
    minIndex = ys.index(min(ys))
    sidePoints = 3
    gr2 = ROOT.TGraph()
    for iNew,iOld in enumerate(range(max(0,minIndex-sidePoints),min(gr.GetN(),minIndex+sidePoints+1))):
        gr2.SetPoint(iNew, points[iOld][1], points[iOld][2])
    gr = gr2
    gr.SetTitle(";Scale;#chi^{2}")
    gr.SetMarkerStyle(20)
    gr.Draw("ap")
    gr.Fit("pol2", "Q")
    fitFunc = gr.GetFunction("pol2")
    fitFunc.SetLineColor(ROOT.kRed)
    parameters = fitFunc.GetParameters()
    a1, a2 = parameters[1], parameters[2]
    fitScale = -a1/(2*a2)
    deltaChi2 = 1 # change chi2 by this value for the uncertainty
    fitErr = aux.sqrt(deltaChi2/a2) if a2>0 else 0
    errFunc = fitFunc.Clone()
    errFunc.SetRange(fitScale-fitErr, fitScale+fitErr)
    errFunc.SetFillColorAlpha(errFunc.GetLineColor(), .5)
    errFunc.SetFillStyle(1001)

    errFunc.Draw("FC same")
    aux.save(saveName+"_fit", "savedFitPredictions/",log=False)

    err = aux.sqrt(fitErr**2 + (1-fitScale)**2)
    fitScaleUp = fitScale + err
    fitScaleDn = fitScale - err
    preHist = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScale), weight, nBins)
    preHistUp = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScaleUp), weight, nBins)
    preHistDn = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScaleDn), weight, nBins)
    syst = aux.getSystFromDifference(preHistDn, preHistUp)

    # Scale
    preInt, preIntErr = aux.integralAndError(preHist, 0, maxBin)
    preIntUp, preIntErrUp = aux.integralAndError(preHistUp, 0, maxBin)
    preIntDn, preIntErrDn = aux.integralAndError(preHistDn, 0, maxBin)
    preIntErr2 = abs(preIntDn-preIntUp)/2.
    relScaleUncert = aux.sqrt( (dirIntErr/dirInt)**2 + (preIntErr/preInt)**2 + (preIntErr2/preInt)**2 )

    for h in preHist, syst:
        h.Scale(dirInt/preInt)
        h.SetDirectory(0)

    # Force symmetric systematic uncertainty
    for bin in range(syst.GetNbinsX()+2):
        syst.SetBinContent(bin, preHist.GetBinContent(bin))
        syst.SetBinError(bin, aux.sqrt(syst.GetBinError(bin)**2 + (relScaleUncert*syst.GetBinContent(bin))**2))
    return preHist, syst, fitScale, fitErr, dirInt/preInt

def qcdClosure(name, dirSet, treename="tr/simpleTree", preSet=None, additionalSets=[], cut="1", noScale=False):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain(treename.replace("tr", "tr_jControl"))

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+[130, 170, 230, 300, 400, 500, 600]
    #range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200, 300, 400, 500, 600]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data:
        for bin in range(dirHist.FindBin(100),dirHist.GetNbinsX()+2):
            dirHist.SetBinContent(bin, 0)
            dirHist.SetBinError(bin, 0)
    preHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins)
    #preHist, gjetSyst = getGJetPrediction(dirTree, preTree, variable, weight, nBins)
    totUnc = aux.addHistUncert(preHist, gjetSyst)

    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weight, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

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
    gjetHistUnw.Divide(preHist)
    gjetHistUnw.SetLineColor(ROOT.kBlue)
    gjetHistUnw.Draw("same hist")
    leg = ROOT.TLegend(.17,.15,.39,.23)
    leg.SetFillStyle(0)
    leg.AddEntry(gjetHistUnw, "Unweighted/Weighted", "l")
    leg.Draw()

    l = aux.Label(sim= not dirSet==data, info=dirSet.label)
    aux.save("qcdClosure_{}".format(name))

def finalDistribution(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1"):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain(treename.replace("tr", "tr_jControl"))
    eTree = ROOT.TChain(treename.replace("tr", "tr_eControl"))

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 450, 500, 550, 600, 610]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 160 and False:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)
    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins)
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weight, nBins, treename.replace("tr", "tr_jControl"))
    mcPreHist.Scale(norm) # scale also with trigger prescale?
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        #gjetHist.SetBinContent(bin, cOld - subT)
        #gjetSyst.SetBinContent(bin, cOld - subT)

    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weight, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.createHistoFromTree(eTree, variable, weight, nBins)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    #eHist.Scale( 0.0 ) # disable data driven bkg
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    # warning: scaling systematic errors are displayed as statistical ones
    zgHist = aux.createHistoFromDatasetTree(zg, variable, weight, nBins, treename)
    wgHist = aux.createHistoFromDatasetTree(wg, variable, weight, nBins, treename)
    ttgHist = aux.createHistoFromDatasetTree(ttg, variable, weight, nBins, treename)
    zHist = aux.createHistoFromDatasetTree(znunu, variable, weight, nBins, treename)
    wHist = aux.createHistoFromDatasetTree(wjets, variable, weight, nBins, treename)
    wHist.Add(aux.createHistoFromDatasetTree(wjets, variable, weight, nBins, treename.replace("/","_genE/")), -1)
    ttHist = aux.createHistoFromDatasetTree(ttjets_ht, variable, weight, nBins, treename)
    ttHist.Add(aux.createHistoFromDatasetTree(ttjets_ht, variable, weight, nBins, treename.replace("/","_genE/")), -1)
    mcSystUncert = 0.3
    zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    wgSyst = aux.getSysHisto(wgHist, mcSystUncert)
    ttgSyst = aux.getSysHisto(ttgHist, mcSystUncert)
    zSyst = aux.getSysHisto(zHist, mcSystUncert)
    wSyst = aux.getSysHisto(wHist, mcSystUncert)
    ttSyst = aux.getSysHisto(ttHist, mcSystUncert)

    totStat = aux.addHists(gjetHist, eHist, zgHist, wgHist, ttgHist, wHist, zHist, ttHist)
    totSyst = aux.addHists(gjetSyst, eSyst, zgSyst, wgSyst, ttgSyst, wSyst, zSyst, ttSyst)
    totUnc = aux.addHistUncert(totStat, totSyst)

    signal = aux.createHistoFromDatasetTree(t5wg_1600_100, variable, weight, nBins, treename)
    aux.drawOpt(signal, "signal")
    #signal.Add(totStat)
    signal.SetLineColor(ROOT.kAzure)

    # beautify
    aux.drawOpt(dirHist, "data")
    aux.drawOpt(totUnc, "totUnc")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add(dirHist, "Data" if dirSet == data else "Pseudodata")
    m.add(signal, "T5Wg 1600 100")
    m.addStack(eHist, "e#rightarrow#gamma")
    m.addStack(zHist, "Z#rightarrow#nu#nu")
    m.addStack(wHist, "W")
    m.addStack(ttHist, "t#bar{t}")
    m.addStack(ttgHist, "t#bar{t}#gamma")
    m.addStack(zgHist, "Z#gamma")
    m.addStack(wgHist, "W#gamma")
    m.addStack(gjetHist, "#gamma + Jet")
    m.add(totUnc, "Tot. uncert.")
    legHeader = ""
    if cut is not "1": legHeader += cut
    if "_ee" in treename: legHeader += " EE"
    if "dPhi" in treename: legHeader += " dPhi"
    if legHeader: m.leg.SetHeader(legHeader)
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, totStat, totSyst)
    x = r.draw(0., 1.5, m.getStack())
    gjetHistUnw.Add(eHist)
    gjetHistUnw.Add(zgHist)
    gjetHistUnw.Add(wgHist)
    gjetHistUnw.Add(ttgHist)
    gjetHistUnw.Add(wHist)
    gjetHistUnw.Add(zHist)
    gjetHistUnw.Add(ttHist)
    gjetHistUnw.Divide(totSyst)
    gjetHistUnw.SetLineColor(ROOT.kBlue)
    gjetHistUnw.Draw("same hist")
    leg = ROOT.TLegend(.17,.12,.39,.18)
    leg.SetFillStyle(0)
    leg.AddEntry(gjetHistUnw, "Unweighted/Weighted", "l")
    leg.Draw()

    #x = aux.drawContributions(m.getStack())

    l = aux.Label(sim= not dirSet==data, info=dirSet.label if dirSet != data else "")
    aux.save("finalDistribution_{}".format(name), normal=False)

    for bin in range(dirHist.GetNbinsX()-5, dirHist.GetNbinsX()+1):
        print
        print "Bin: {} to {}".format(dirHist.GetBinLowEdge(bin),dirHist.GetBinLowEdge(bin+1))
        print "Data:", dirHist.GetBinContent(bin)
        print "Bkg :", totStat.GetBinContent(bin)
        print "gJet:", gjetHist.GetBinContent(bin)
        print "e->g:", eHist.GetBinContent(bin)
        print "wg:  ", wgHist.GetBinContent(bin)
        print "zg:  ", zgHist.GetBinContent(bin)
        print "ttg: ", ttgHist.GetBinContent(bin)
        print "w:   ", wHist.GetBinContent(bin)
        print "z:   ", zHist.GetBinContent(bin)
        print "tt:  ", ttHist.GetBinContent(bin)


def finalDistributionEasy(name, dirSet, preSet=None):
    if not preSet: preSet = dirSet
    hname = "met"
    folder = "tr"
    binning = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    dirHist = aux.stdHist(dirSet, "{}/{}".format(folder,hname), binning)
    if dirSet == data and False:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 160:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)
    gjetHist = aux.stdHist(preSet, "{}/{}".format("tr_jControl",hname), binning)
    maxScaleBin = dirHist.FindFixBin(100)-1
    gjetHist.Scale(dirHist.Integral(0,maxScaleBin)/gjetHist.Integral(0,maxScaleBin))
    gjetSyst = aux.getSysHisto(gjetHist, 0.2)
    gjetHist.SetLineColor(ROOT.kCyan)

    eHist = aux.stdHist(dirSet, "{}/{}".format("tr_eControl",hname), binning)
    eHist.Scale( 0.0197 if dirSet==data else 0.0107 )
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    zgHist = aux.stdHist(zg+znunu, "{}/{}".format(folder,hname), binning)
    wgHist = aux.stdHist(wg+wjets, "{}/{}".format(folder,hname), binning)
    ttgHist = aux.stdHist(ttg+ttjets, "{}/{}".format(folder,hname), binning)
    #for h in zgHist, wgHist, ttgHist:
    #    h.Scale(1.5)
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
    m.add(dirHist, "Data" if dirSet == data else "Pseudodata")
    m.addStack(eHist, "e#rightarrow#gamma")
    m.addStack(ttgHist, "t#bar{t}#gamma")
    m.addStack(zgHist, "Z#gamma")
    m.addStack(wgHist, "W#gamma")
    m.addStack(gjetHist, "#gamma + Jet")
    m.add(totUnc, "Tot. uncert.")
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, totStat, totSyst)
    r.draw(0.5, 1.5)

    l = aux.Label(sim= not dirSet==data, info=dirSet.label if dirSet != data else "")
    aux.save("finalDistributionEasy_{}".format(name))

def emhtReweighting(treename="tr/simpleTree"):
    dataSet = data
    if "jControl" in treename: dataSet = dataHt
    mcSet = gjets #+ qcd + zg + wg_130 + ttg + wjets + ttjets
    mcSet = qcd

    dataTree = ROOT.TChain(treename)
    mcTree = ROOT.TChain(treename)

    for f in dataSet.files: dataTree.Add(f)
    for f in mcSet.files: mcTree.Add(f)

    cut = "1"
    variable = "emht"
    weight = "weight*({})".format(cut)
    #nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    #nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,1300,10)+range(1300,2200,50)+range(2200,3000,100)

    dataHist = aux.createHistoFromTree(dataTree, variable, weight, nBins)
    mcHist = aux.createHistoFromTree(mcTree, variable, weight, nBins)
    mcHist.Scale(dataHist.Integral()/mcHist.Integral())

    dataHist.Divide(mcHist)
    dataHist.Draw()
    f1 = ROOT.TF1("f1","pol2", 0, 3000)
    dataHist.Fit("f1")
    aux.save("test", endings=[".pdf", ".root"], log=False)
    return f1

def plotOverlayedPredictions(filename):
    f = ROOT.TFile(filename)
    hdir = f.Get("dir")
    hdir.GetXaxis().SetRangeUser(0,100)
    m = multiplot.Multiplot()
    m.drawOption_ = "ep"
    m.add(hdir, "Pseudo(Data)")
    for key in sorted([k for k in f.GetListOfKeys()]):
        name = key.GetName()
        if name == "dir": continue
        h = key.ReadObj()
        h.GetXaxis().SetRangeUser(0,100)
        h.Scale(hdir.Integral()/h.Integral())
        h.SetLineColor(ROOT.kRed)
        h.SetMarkerColor(ROOT.kRed)
        if name == "1.0":
            h.SetLineColor(ROOT.kBlue)
            h.SetMarkerColor(ROOT.kBlue)
        m.add(h, "Jet E_{T}^{miss} scaled by "+name)
    m.Draw()
    aux.save(filename.replace(".root","_overlay").replace("/","_"))


if __name__ == "__main__":
    allMC = gjets+qcd+zg+wg+ttg+wjets+ttjets_ht+znunu
    allMC.label = "MC mix"
    #qcdClosure("gqcd", gjets+qcd)
    #finalDistribution("mc", allMC, allMC)
    #finalDistribution("mc_2000emht_dPhi3", allMC, allMC, cut="2000<emht&&dPhi>.3")
    #finalDistribution("mc_2000emht", allMC, allMC, cut="emht>2000")
    #finalDistribution("mc_ee", allMC, allMC, treename="tr_ee/simpleTree")
    #finalDistribution("mc_ee_2000emht", allMC, allMC, treename="tr_ee/simpleTree", cut="emht>2000")

    #finalDistribution("data", data, dataHt)
    #finalDistribution("data_ee", data, dataHt, "tr_ee/simpleTree")
    #finalDistribution("data_2000emht", data, dataHt, cut="2000<emht")
    finalDistribution("data_ee_2000emht", data, dataHt, "tr_ee/simpleTree", "2000<emht")
    #finalDistribution("data_dPhi3", data, dataHt, cut="dPhi<.3")
    #finalDistribution("data_ee_dPhi3", data, dataHt, "tr_ee/simpleTree", "dPhi<.3")
    #finalDistribution("data_3dPhi", data, dataHt, cut=".3<dPhi")
    #finalDistribution("data_ee_3dPhi", data, dataHt, "tr_ee/simpleTree", ".3<dPhi")

    #finalDistribution("data_ee", data, dataHt, "tr_ee/simpleTree")
    #finalDistribution("data_ee_dPhi3", data, dataHt, "tr_ee/simpleTree", cut="dPhi<.3")
    #finalDistribution("data_ee_2000emht_dPhi3", data, dataHt, "tr_ee/simpleTree", cut="2000<emht&&dPhi<.3")
    #finalDistribution("data", data, dataHt)
    #finalDistribution("data_dphi1", data, dataHt, cut="dPhi>.1")
    #finalDistribution("data_dphi2", data, dataHt, cut="dPhi>.2")
    #finalDistribution("data_dphi3", data, dataHt, cut="dPhi>.3")
    #finalDistribution("data_dphi4", data, dataHt, cut="dPhi>.4")
    #finalDistribution("data_emht1000", data, dataHt, cut="emht<1000")
    #finalDistribution("data_1500emht", data, dataHt, cut="emht>1500")
    #finalDistribution("data_2000emht_dphi3", data, dataHt, cut="2000<emht&&dPhi>.3")
    #qcdClosure("data", data, preSet=dataHt)
    #qcdClosure("data_2000emht", data, preSet=dataHt, cut="emht>2000")
    #qcdClosure("gqcd_emht1000", gjets+qcd, cut="emht<1000")
    #qcdClosure("gqcd_1000emht2000", gjets+qcd, cut="1000<emht&&emht<2000")
    #qcdClosure("gqcd_1000emht", gjets+qcd, cut="emht>1000")
    #qcdClosure("gqcd_1500emht", gjets+qcd, cut="emht>1500")
    #qcdClosure("gqcd_1700emht", gjets+qcd, cut="emht>1700")
    #qcdClosure("gqcd_1800emht", gjets+qcd, cut="emht>1800")
    #qcdClosure("gqcd_1900emht", gjets+qcd, cut="emht>1900")
    #qcdClosure("gqcd_2000emht", gjets+qcd, cut="emht>2000")
    #qcdClosure("gqcd_2100emht", gjets+qcd, cut="emht>2100")
    #qcdClosure("gqcd_2200emht", gjets+qcd, cut="emht>2200")
    #qcdClosure("gqcd_2300emht", gjets+qcd, cut="emht>2300")
    #qcdClosure("gqcd_2500emht", gjets+qcd, cut="emht>2500")
    #finalDistributionEasy("data", data, dataHt)
    #plotOverlayedPredictions("savedFitPredictions/gqcd__8__tr_simpleTree__8__weight1__met__24.root")
    #qcdClosure("gqcd_ee_2000emht", gjets+qcd, treename="tr_ee/simpleTree", cut="emht>2000")

    """
    qcdClosure("gqcd_700emht800", gjets+qcd, cut="700<emht && emht<800")
    qcdClosure("gqcd_800emht900", gjets+qcd, cut="800<emht && emht<900")
    qcdClosure("gqcd_900emht1000", gjets+qcd, cut="900<emht && emht<1000")
    qcdClosure("gqcd_1000emht1100", gjets+qcd, cut="1000<emht && emht<1100")
    qcdClosure("gqcd_1100emht1200", gjets+qcd, cut="1100<emht && emht<1200")
    qcdClosure("gqcd_1200emht1300", gjets+qcd, cut="1200<emht && emht<1300")
    qcdClosure("gqcd_1300emht1400", gjets+qcd, cut="1300<emht && emht<1400")
    qcdClosure("gqcd_1400emht1500", gjets+qcd, cut="1400<emht && emht<1500")
    qcdClosure("gqcd_1500emht1600", gjets+qcd, cut="1500<emht && emht<1600")
    qcdClosure("gqcd_1600emht1700", gjets+qcd, cut="1600<emht && emht<1700")
    qcdClosure("gqcd_1700emht1800", gjets+qcd, cut="1700<emht && emht<1800")
    qcdClosure("gqcd_1800emht1900", gjets+qcd, cut="1800<emht && emht<1900")
    qcdClosure("gqcd_1900emht2000", gjets+qcd, cut="1900<emht && emht<2000")
    qcdClosure("gqcd_2000emht2100", gjets+qcd, cut="2000<emht && emht<2100")
    qcdClosure("gqcd_2100emht2200", gjets+qcd, cut="2100<emht && emht<2200")
    qcdClosure("gqcd_2200emht", gjets+qcd, cut="2200<emht")
    qcdClosure("data_700emht800", data, preSet=dataHt, cut="700<emht && emht<800")
    qcdClosure("data_800emht900", data, preSet=dataHt, cut="800<emht && emht<900")
    qcdClosure("data_900emht1000", data, preSet=dataHt, cut="900<emht && emht<1000")
    qcdClosure("data_1000emht1100", data, preSet=dataHt, cut="1000<emht && emht<1100")
    qcdClosure("data_1100emht1200", data, preSet=dataHt, cut="1100<emht && emht<1200")
    qcdClosure("data_1200emht1300", data, preSet=dataHt, cut="1200<emht && emht<1300")
    qcdClosure("data_1300emht1400", data, preSet=dataHt, cut="1300<emht && emht<1400")
    qcdClosure("data_1400emht1500", data, preSet=dataHt, cut="1400<emht && emht<1500")
    qcdClosure("data_1500emht1600", data, preSet=dataHt, cut="1500<emht && emht<1600")
    qcdClosure("data_1600emht1700", data, preSet=dataHt, cut="1600<emht && emht<1700")
    qcdClosure("data_1700emht1800", data, preSet=dataHt, cut="1700<emht && emht<1800")
    qcdClosure("data_1800emht1900", data, preSet=dataHt, cut="1800<emht && emht<1900")
    qcdClosure("data_1900emht2000", data, preSet=dataHt, cut="1900<emht && emht<2000")
    qcdClosure("data_2000emht2100", data, preSet=dataHt, cut="2000<emht && emht<2100")
    qcdClosure("data_2200emht", data, preSet=dataHt, cut="2200<emht")
    """
