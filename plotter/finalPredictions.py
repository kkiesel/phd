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

def qcdClosure(name, dirSet, treename="tr/simpleTree", preSet=None, additionalSets=[], cut="1", noScale=False, binning=None):
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
    if binning: nBins = binning

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data:
        for bin in range(dirHist.FindBin(100),dirHist.GetNbinsX()+2):
            dirHist.SetBinContent(bin, 0)
            dirHist.SetBinError(bin, 0)
    preHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins)
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
    legHeader = ""
    if cut is not "1": legHeader += aux.beautifyCutString(cut)
    if "_ee" in treename: legHeader += " EE"
    if legHeader: m.leg.SetHeader(legHeader)
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, preHist, gjetSyst)
    r.draw(0.5, 1.5)
    gjetHistUnw.Divide(preHist)
    gjetHistUnw.SetLineColor(ROOT.kGreen)
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
    nBins = range(0,200,10)+[200, 250, 300, 350, 450, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 100 and "_ee" not in name and False:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins)
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weight, nBins, treename.replace("tr", "tr_jControl"))
    mcPreHist.Scale(norm)
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        gjetHist.SetBinContent(bin, cOld - subT)
        gjetSyst.SetBinContent(bin, cOld - subT)
    """
    emhtRange = [700,750,800,900,1000,1500,2000,25000]
    gjetHist = None
    gjetSyst = None
    for emhtBin in range(len(emhtRange)-1):
        emhtSliceCut = "&&{}<emht&&emht<{}".format(emhtRange[emhtBin],emhtRange[emhtBin+1])
        gjetHist_i, gjetSyst_i, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight.replace(")", emhtSliceCut+")"), variable, nBins)
        gjetHist_i.SetLineColor(ROOT.kCyan)
        # correct for other backgrounds
        mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weight.replace(")", emhtSliceCut+")"), nBins, treename.replace("tr", "tr_jControl"))
        mcPreHist.Scale(norm)
        mcPreHist.Scale(5.93369235121/36.459205898)
        for bin in range(gjetHist_i.GetNbinsX()+2):
            cOld = gjetHist_i.GetBinContent(bin)
            subT = mcPreHist.GetBinContent(bin)
            gjetHist_i.SetBinContent(bin, cOld - subT)
            gjetSyst_i.SetBinContent(bin, cOld - subT)
        if gjetHist: gjetHist.Add(gjetHist_i)
        else: gjetHist = gjetHist_i
        if gjetSyst:
#            gjetSyst.Add(gjetSyst_i)
            for bin in range(gjetSyst.GetNbinsX()+2):
                gjetSyst.SetBinContent(bin, gjetSyst.GetBinContent(bin)+gjetSyst_i.GetBinContent(bin))
                gjetSyst.SetBinError(bin, gjetSyst.GetBinError(bin)+gjetSyst_i.GetBinError(bin))
        else: gjetSyst = gjetSyst_i
    """
    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weight, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.createHistoFromTree(eTree, variable, weight, nBins)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
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
    m.addStack(gjetHist, "#gamma + Jet (datadriven)")
    m.add(totUnc, "Tot. uncert.")
    legHeader = ""
    if cut is not "1": legHeader += aux.beautifyCutString(cut)
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
    gjetHistUnw.SetLineColor(ROOT.kGreen)
    gjetHistUnw.Draw("same hist")
    leg = ROOT.TLegend(.17,.12,.39,.18)
    leg.SetFillStyle(0)
    leg.AddEntry(gjetHistUnw, "Unweighted/Weighted", "l")
    leg.Draw()

    #x = aux.drawContributions(m.getStack())

    l = aux.Label(sim= not dirSet==data, info=dirSet.label if dirSet != data else "")
    aux.save("finalDistribution_{}".format(name), normal=False)
    return
    dc = limitTools.MyDatacard("photonHt")
    for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
        dc.addBin(
            "bin"+str(bin), dirHist.GetBinContent(bin),
            gjetHist.GetBinContent(bin), 1.+gjetHist.GetBinError(bin)/gjetHist.GetBinContent(bin), 1.+gjetSyst.GetBinError(bin)/gjetHist.GetBinContent(bin),
            eHist.GetBinContent(bin), 1.+eHist.GetBinError(bin)/eHist.GetBinContent(bin), 1.+eSyst.GetBinError(bin)/eHist.GetBinContent(bin),
            zgHist.GetBinContent(bin)+zHist.GetBinContent(bin), 1.+aux.sqrt(zHist.GetBinError(bin)**2+zgHist.GetBinError(bin)**2)/(zgHist.GetBinContent(bin)+zHist.GetBinContent(bin)), 1.+aux.sqrt(zSyst.GetBinError(bin)**2+zgSyst.GetBinError(bin)**2)/(zgHist.GetBinContent(bin)+zHist.GetBinContent(bin)),
            wgHist.GetBinContent(bin)+wHist.GetBinContent(bin), 1.+aux.sqrt(wHist.GetBinError(bin)**2+wgHist.GetBinError(bin)**2)/(wgHist.GetBinContent(bin)+wHist.GetBinContent(bin)), 1.+aux.sqrt(wSyst.GetBinError(bin)**2+wgSyst.GetBinError(bin)**2)/(wgHist.GetBinContent(bin)+wHist.GetBinContent(bin)),
            ttgHist.GetBinContent(bin)+ttHist.GetBinContent(bin), 1.+aux.sqrt(ttHist.GetBinError(bin)**2+ttgHist.GetBinError(bin)**2)/(ttgHist.GetBinContent(bin)+ttHist.GetBinContent(bin)), 1.+aux.sqrt(ttSyst.GetBinError(bin)**2+ttgSyst.GetBinError(bin)**2)/(ttgHist.GetBinContent(bin)+ttHist.GetBinContent(bin)),
            signal.GetBinContent(bin), 1.+signal.GetBinError(bin)/signal.GetBinContent(bin), 1.3)
    #print dc.limit()

    return
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

def plotOverlayedPredictions(filename):
    c = ROOT.TCanvas()
    f = ROOT.TFile(filename)
    hdir = f.Get("dir")
    hdir.GetXaxis().SetRangeUser(0,100)
    aux.drawOpt(hdir, "data")
    m = multiplot.Multiplot()
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
            m.add(h, "Jet selection")
        elif name == "1.1":
            m.add(h, "Jet selection, shifted #it{E}_{T}^{miss}")
        else:
            m.add(h, "")
    m.Draw()
    l = aux.Label(sim="data" not in filename)
    aux.save(filename.replace(".root","_overlay").replace("/","_"))

def gjetPrediction(dirHist, preSet, subSet, variable, nBins, weight="weight", saveName=""):
    dirHist = dirHist.Clone("dir")
    saveName = "savedFitPredictions/weight_{}_{}".format(aux.modifySaveName(saveName),aux.modifySaveName(weight))
    saveNameRoot = saveName + ".root"
    if os.path.isfile(saveNameRoot):
        print "Using saved events from {}".format(saveNameRoot)
        f = ROOT.TFile(saveNameRoot)
        dirHist = f.Get("dir")
        preHists = {}
        for key in f.GetListOfKeys():
            name = key.GetName()
            if name == "dir": continue
            preHists[float(name)] = key.ReadObj()
    else:
        print "Calculating new", saveNameRoot
        scales = [.85+.01*i for i in range(30)]
        preHists = {}
        for iscale, scale in enumerate(scales):
            preHist   = aux.createHistoFromDatasetTree(preSet, "{}*{}".format(variable,scale), weight, nBins, "tr_jControl/simpleTree")
            mcPreHist = aux.createHistoFromDatasetTree(subSet, "{}*{}".format(variable,scale), weight, nBins, "tr_jControl/simpleTree")
            preHist.Add(mcPreHist, -1)
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
    preHist = aux.createHistoFromDatasetTree(preSet, "{}*{}".format(variable,fitScale), weight, nBins, "tr_jControl/simpleTree")
    preHistUp = aux.createHistoFromDatasetTree(preSet, "{}*{}".format(variable,fitScaleUp), weight, nBins, "tr_jControl/simpleTree")
    preHistDn = aux.createHistoFromDatasetTree(preSet, "{}*{}".format(variable,fitScaleDn), weight, nBins, "tr_jControl/simpleTree")
    mcPreHist = aux.createHistoFromDatasetTree(subSet, "{}*{}".format(variable,fitScale), weight, nBins, "tr_jControl/simpleTree")
    mcPreHistUp = aux.createHistoFromDatasetTree(subSet, "{}*{}".format(variable,fitScaleUp), weight, nBins, "tr_jControl/simpleTree")
    mcPreHistDn = aux.createHistoFromDatasetTree(subSet, "{}*{}".format(variable,fitScaleDn), weight, nBins, "tr_jControl/simpleTree")
    preHist.Add(mcPreHist, -1)
    preHistUp.Add(mcPreHistUp, -1)
    preHistDn.Add(mcPreHistDn, -1)

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
    return preHist, syst

def finalDistribution1dHist(name, dirSet, preSet):
    if not preSet: preSet = dirSet

    nBins = range(0,200,10)+[200, 250, 300, 350, 450, 600, 700]
    dirHist = aux.stdHist(dirSet, name, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst = gjetPrediction(dirHist, preSet, zg+wg+ttg+wjets+ttjets_nlo+znunu, "met", nBins, "weight", name)
    gjetHist.SetLineColor(rwth.myLightBlue)

    #gjetHistUnw = aux.stdHist(preSet, "tr_jControl/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    #gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    #eHist = aux.stdHist(dirSet, "tr_eControl/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    #eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    #eHist.SetLineColor(rwth.myYellow)
    #eSyst = aux.getSysHisto(eHist, 0.3)

    zgHist = aux.stdHist(zg+znunu, name, nBins)
    wgHist = aux.stdHist(wg+wjets, name, nBins)
    tgHist = aux.stdHist(ttjets_nlo+ttg+tg, name, nBins)
    zgHist.SetLineColor(rwth.myRed)
    wgHist.SetLineColor(rwth.myOrange)
    tgHist.SetLineColor(rwth.myBlue)
    # todo: suptracet electron stuff
    mcSystUncert = 0.3
    zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    wgSyst = aux.getSysHisto(wgHist, mcSystUncert)
    tgSyst = aux.getSysHisto(tgHist, mcSystUncert)

    totStat = aux.addHists(gjetHist, zgHist, wgHist, tgHist)
    totSyst = aux.addHists(gjetSyst, zgSyst, wgSyst, tgSyst)
    totUnc = aux.addHistUncert(totStat, totSyst)

    signal = aux.stdHist(t5wg_1600_100, name, nBins)

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
    m.addStack(zgHist, "#gammaZ")
    m.addStack(tgHist, "#gammat#bar{t}")
    m.addStack(wgHist, "#gammaW")
    m.addStack(gjetHist, "#gamma + Jet (datadriven)")
    m.add(totUnc, "Tot. uncert.")
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, totStat, totSyst)
    x = r.draw(0., 1.5, m.getStack())
    l = aux.Label(sim= not dirSet==data, info=dirSet.label if dirSet != data else "")
    aux.save("finalDistribution1dHist_{}".format(name), normal=False)


def finalDistributionHist(name, dirSet, preSet, minEmht=0, maxEmht=1e8, directory="tr"):
    if not preSet: preSet = dirSet

    nBins = range(0,200,10)+[200, 250, 300, 350, 450, 600, 700]
    dirHist = aux.stdHist(dirSet, "tr/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    weight = "weight*({}<emht&&emht<{})".format(minEmht, maxEmht)
    gjetHist, gjetSyst = gjetPrediction(dirHist, preSet, zg+wg+ttg+wjets+ttjets_nlo+znunu, "met", nBins, weight)
    gjetHist.SetLineColor(rwth.myLightBlue)

    gjetHistUnw = aux.stdHist(preSet, "tr_jControl/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.stdHist(dirSet, "tr_eControl/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    eHist.SetLineColor(rwth.myYellow)
    eSyst = aux.getSysHisto(eHist, 0.3)

    zgHist = aux.stdHist(zg+znunu, "tr/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    wgHist = aux.stdHist(wg+wjets, "tr/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    tgHist = aux.stdHist(ttjets_nlo+ttg+tg, "tr/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)
    zgHist.SetLineColor(rwth.myRed)
    wgHist.SetLineColor(rwth.myOrange)
    tgHist.SetLineColor(rwth.myBlue)
    # todo: suptracet electron stuff
    mcSystUncert = 0.3
    zgSyst = aux.getSysHisto(zgHist, mcSystUncert)
    wgSyst = aux.getSysHisto(wgHist, mcSystUncert)
    tgSyst = aux.getSysHisto(tgHist, mcSystUncert)

    totStat = aux.addHists(gjetHist, eHist, zgHist, wgHist, tgHist)
    totSyst = aux.addHists(gjetSyst, eSyst, zgSyst, wgSyst, tgSyst)
    totUnc = aux.addHistUncert(totStat, totSyst)

    signal = aux.stdHist(t5wg_1600_100, "tr/met_vs_emht", nBins, False, minEmht, maxEmht-1e-6)

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
    m.addStack(eHist, "e#rightarrow#gamma (datadriven)")
    m.addStack(zgHist, "#gammaZ")
    m.addStack(tgHist, "#gammat#bar{t}")
    m.addStack(wgHist, "#gammaW")
    m.addStack(gjetHist, "#gamma + Jet (datadriven)")
    m.add(totUnc, "Tot. uncert.")
    m.Draw()

    r = ratio.Ratio("Data/Pred", dirHist, totStat, totSyst)
    x = r.draw(0., 1.5, m.getStack())
    gjetHistUnw.Add(eHist)
    gjetHistUnw.Add(zgHist)
    gjetHistUnw.Add(wgHist)
    gjetHistUnw.Add(tgHist)
    gjetHistUnw.Divide(totSyst)
    gjetHistUnw.SetLineColor(ROOT.kGreen)
    gjetHistUnw.Draw("same hist")
    leg = ROOT.TLegend(.17,.12,.39,.18)
    leg.SetFillStyle(0)
    leg.AddEntry(gjetHistUnw, "Unweighted/Weighted", "l")
    leg.Draw()
    l = aux.Label(sim= not dirSet==data, info=dirSet.label if dirSet != data else "")
    aux.save("finalDistribution_{}".format(name), normal=False)

    if name == "data_emht2000": dc = limitTools.MyDatacard("photonHt")
    elif name == "data_2000emht": dc = limitTools.MyDatacard("test.txt")
    for bin in range(dirHist.GetNbinsX()-2, dirHist.GetNbinsX()+1):
        dc.addBin(
            name.split("_")[1]+"_"+str(bin), dirHist.GetBinContent(bin),
            gjetHist.GetBinContent(bin), 1.+gjetHist.GetBinError(bin)/gjetHist.GetBinContent(bin), 1.+gjetSyst.GetBinError(bin)/gjetHist.GetBinContent(bin),
            eHist.GetBinContent(bin), 1.+eHist.GetBinError(bin)/eHist.GetBinContent(bin), 1.+eSyst.GetBinError(bin)/eHist.GetBinContent(bin),
            zgHist.GetBinContent(bin), 1.+zgHist.GetBinError(bin)/zgHist.GetBinContent(bin), 1.+zgSyst.GetBinError(bin)/zgHist.GetBinContent(bin),
            wgHist.GetBinContent(bin), 1.+wgHist.GetBinError(bin)/wgHist.GetBinContent(bin), 1.+wgSyst.GetBinError(bin)/wgHist.GetBinContent(bin),
            tgHist.GetBinContent(bin), 1.+tgHist.GetBinError(bin)/tgHist.GetBinContent(bin), 1.+tgSyst.GetBinError(bin)/tgHist.GetBinContent(bin),
            signal.GetBinContent(bin), 1.+signal.GetBinError(bin)/signal.GetBinContent(bin), 1.3)
    dc.write("testDatacard.txt")


if __name__ == "__main__":
    allMC = gjets+qcd+zg+wg+ttg+wjets+ttjets_ht+znunu
    allMC.label = "MC mix"
    #qcdClosure("gqcd", gjets+qcd)
    #finalDistribution("mc", allMC, allMC)
    #finalDistribution("mc_2000emht_dPhi3", allMC, allMC, cut="2000<emht&&dPhi>.3")
    #finalDistribution("mc_2000emht", allMC, allMC, cut="emht>2000")
    #finalDistribution("mc_ee", allMC, allMC, treename="tr_ee/simpleTree")
    #finalDistribution("mc_ee_2000emht", allMC, allMC, treename="tr_ee/simpleTree", cut="emht>2000")
    #dataB= Dataset("SinglePhoton_Run2016B-23Sep2016-v3", 0, ROOT.kBlack )
    #dataB.label = "Data RunB: 5.9/fb"
    #dataHtB = Dataset("JetHT_Run2016B-23Sep2016-v3", 0, ROOT.kBlack )

    #finalDistribution("dataB", dataB, dataHtB)
    #finalDistribution("dataB_emhtbins_syst", dataB, dataHtB)
    #finalDistribution("dataB_ee", dataB, dataHtB, "tr_ee/simpleTree")
    #finalDistribution("dataB_2000emht", dataB, dataHtB, cut="2000<emht")
    #finalDistribution("dataB_emht800", dataB, dataHtB, cut="emht<800")
    #finalDistribution("dataB_800emht900", dataB, dataHtB, cut="800<emht&&emht<900")
    #finalDistribution("dataB_900emht1000", dataB, dataHtB, cut="900<emht&&emht<1000")
    #finalDistribution("data", data, dataHt)
    #finalDistribution("data_ee", data, dataHt, "tr_ee/simpleTree")
    #finalDistribution("data_emht2000", data, dataHt, cut="emht<2000")
    #finalDistribution("data_2000emht", data, dataHt, cut="2000<emht")
    #finalDistribution("data_ee_2000emht", data, dataHt, "tr_ee/simpleTree", "2000<emht")
    #finalDistribution("data_emht1000", data, dataHt, cut="emht<1000")
    #finalDistribution("data_ee_emht1000", data, dataHt, "tr_ee/simpleTree", "emht<1000")
    #finalDistribution("data_1000emht2000", data, dataHt, cut="1000<emht&&emht<2000")
    #finalDistribution("data_ee_1000emht2000", data, dataHt, "tr_ee/simpleTree", "1000<emht&&emht<2000")

    #finalDistribution("data_dPhi3", data, dataHt, cut="dPhi<.3")
    #finalDistribution("data_ee_dPhi3", data, dataHt, "tr_ee/simpleTree", "dPhi<.3")
    #finalDistribution("data_3dPhi", data, dataHt, cut=".3<dPhi")
    #finalDistribution("data_ee_3dPhi", data, dataHt, "tr_ee/simpleTree", ".3<dPhi")
    #finalDistribution("data_dPhi2", data, dataHt, cut="dPhi<.2")
    #finalDistribution("data_ee_dPhi2", data, dataHt, "tr_ee/simpleTree", "dPhi<.2")
    #finalDistribution("data_2dPhi", data, dataHt, cut=".2<dPhi")
    #finalDistribution("data_ee_2dPhi", data, dataHt, "tr_ee/simpleTree", ".2<dPhi")
    #finalDistribution("data_dPhi1", data, dataHt, cut="dPhi<.1")
    #finalDistribution("data_ee_dPhi1", data, dataHt, "tr_ee/simpleTree", "dPhi<.1")
    #finalDistribution("data_1dPhi", data, dataHt, cut=".1<dPhi")
    #finalDistribution("data_ee_1dPhi", data, dataHt, "tr_ee/simpleTree", ".1<dPhi")

    #finalDistribution("data_700emht800", data, dataHt, cut="700<emht && emht<800")
    #finalDistribution("data_800emht900", data, dataHt, cut="800<emht && emht<900")
    #finalDistribution("data_900emht1000", data, dataHt, cut="900<emht && emht<1000")
    #finalDistribution("data_1000emht1100", data, dataHt, cut="1000<emht && emht<1100")
    #finalDistribution("data_1100emht1200", data, dataHt, cut="1100<emht && emht<1200")
    #finalDistribution("data_1200emht1300", data, dataHt, cut="1200<emht && emht<1300")
    #finalDistribution("data_1300emht1400", data, dataHt, cut="1300<emht && emht<1400")
    #finalDistribution("data_1400emht1500", data, dataHt, cut="1400<emht && emht<1500")
    #finalDistribution("data_1500emht1600", data, dataHt, cut="1500<emht && emht<1600")
    #finalDistribution("data_1600emht1700", data, dataHt, cut="1600<emht && emht<1700")
    #finalDistribution("data_1700emht1800", data, dataHt, cut="1700<emht && emht<1800")
    #finalDistribution("data_1800emht1900", data, dataHt, cut="1800<emht && emht<1900")
    #finalDistribution("data_1900emht2000", data, dataHt, cut="1900<emht && emht<2000")
    #finalDistribution("data_2000emht2100", data, dataHt, cut="2000<emht && emht<2100")

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
    #qcdClosure("gqcd", gjets+qcd)
    #qcdClosure("gqcd_2000emht", gjets+qcd, cut="2000<emht", binning=range(0,200,10)+[200, 300, 400])
    #qcdClosure("gqcd_emht1000", gjets+qcd, cut="emht<1000")
    #qcdClosure("gqcd_1000emht2000", gjets+qcd, cut="1000<emht&&emht<2000")
    #qcdClosure("gqcd_ee", gjets+qcd, "tr_ee/simpleTree")
    #qcdClosure("gqcd_ee_2000emht", gjets+qcd, "tr_ee/simpleTree", cut="2000<emht", binning=range(0,200,10)+[200, 300, 400])
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
    #plotOverlayedPredictions("savedFitPredictions/data__8__tr_simpleTree__8__weight1__met__11.root")
    #plotOverlayedPredictions("savedFitPredictions/gqcd__8__tr_simpleTree__8__weight1__met__11.root")
    #qcdClosure("gqcd_ee_2000emht", gjets+qcd, treename="tr_ee/simpleTree", cut="emht>2000")
    #qcdClosure("test", gjets600)

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
    #finalDistributionHist("data_emht2000", data, dataHt, 0, 2000)
    #finalDistributionHist("data_2000emht", data, dataHt, 2000)
    #finalDistribution1dHist("tr/met", data, dataHt)
    #finalDistribution1dHist("tr/metCrystalSeedCorrected", data, dataHt)
    #finalDistribution1dHist("tr/metCrystalSeedCorrected2", data, dataHt)
    #finalDistribution1dHist("tr_100mt/met", data, dataHt)
    #finalDistribution1dHist("tr_metp0/met", data, dataHt)
    finalDistribution1dHist("tr_central/met", data, dataHt)
    finalDistribution1dHist("tr_EB_forward/met", data, dataHt)
