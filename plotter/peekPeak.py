#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def getGJetFitPredictionEGscaled(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins, weightPre="weight", variableScaled="met"):
    weightName = weight.replace("(","").replace(")","").replace("*","").replace(">","").replace("<","").replace("&&","").replace(" ", "").replace("/","DIV")
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
        dirHist = aux.createHistoFromTree(dirTree, variableScaled, weight, nBins)
        dirHist.SetName("dir")
        scales = [.85+.01*i for i in range(30)]
        preHists = {}
        for iscale, scale in enumerate(scales):
            preHist = aux.createHistoFromTree(preTree, "{}*{}".format(variable,scale), weightPre, nBins)
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
        if not preHist.Integral(0,maxBin): return 0, 0, 0, 0, 0
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
    #gr = gr2
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
    preHist = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScale), weightPre, nBins)
    preHistUp = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScaleUp), weightPre, nBins)
    preHistDn = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScaleDn), weightPre, nBins)
    syst = aux.getSystFromDifference(preHistDn, preHistUp)

    # Scale
    preInt, preIntErr = aux.integralAndError(preHist, 0, maxBin)
    preIntUp, preIntErrUp = aux.integralAndError(preHistUp, 0, maxBin)
    preIntDn, preIntErrDn = aux.integralAndError(preHistDn, 0, maxBin)
    preIntErr2 = abs(preIntDn-preIntUp)/2.
    if not dirInt or not preInt: return 0, 0, 0, 0, 0
    relScaleUncert = aux.sqrt( (dirIntErr/dirInt)**2 + (preIntErr/preInt)**2 + (preIntErr2/preInt)**2 )

    for h in preHist, syst:
        h.Scale(dirInt/preInt)
        h.SetDirectory(0)

    # Force symmetric systematic uncertainty
    for bin in range(syst.GetNbinsX()+2):
        syst.SetBinContent(bin, preHist.GetBinContent(bin))
        syst.SetBinError(bin, aux.sqrt(syst.GetBinError(bin)**2 + (relScaleUncert*syst.GetBinContent(bin))**2))
    return preHist, syst, fitScale, fitErr, dirInt/preInt


def getGJetFitPredictionTwoCuts(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins, weightPre="weight"):
    weightName = weight.replace("(","").replace(")","").replace("*","").replace(">","").replace("<","").replace("&&","").replace(" ", "").replace("/","DIV")
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
            preHist = aux.createHistoFromTree(preTree, "{}*{}".format(variable,scale), weightPre, nBins)
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
        if not preHist.Integral(0,maxBin): return 0, 0, 0, 0, 0
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
    preHist = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScale), weightPre, nBins)
    preHistUp = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScaleUp), weightPre, nBins)
    preHistDn = aux.createHistoFromTree(preTree, "{}*{}".format(variable,fitScaleDn), weightPre, nBins)
    syst = aux.getSystFromDifference(preHistDn, preHistUp)

    # Scale
    preInt, preIntErr = aux.integralAndError(preHist, 0, maxBin)
    preIntUp, preIntErrUp = aux.integralAndError(preHistUp, 0, maxBin)
    preIntDn, preIntErrDn = aux.integralAndError(preHistDn, 0, maxBin)
    preIntErr2 = abs(preIntDn-preIntUp)/2.
    if not dirInt or not preInt: return 0, 0, 0, 0, 0
    relScaleUncert = aux.sqrt( (dirIntErr/dirInt)**2 + (preIntErr/preInt)**2 + (preIntErr2/preInt)**2 )

    for h in preHist, syst:
        h.Scale(dirInt/preInt)
        h.SetDirectory(0)

    # Force symmetric systematic uncertainty
    for bin in range(syst.GetNbinsX()+2):
        syst.SetBinContent(bin, preHist.GetBinContent(bin))
        syst.SetBinError(bin, aux.sqrt(syst.GetBinError(bin)**2 + (relScaleUncert*syst.GetBinContent(bin))**2))
    return preHist, syst, fitScale, fitErr, dirInt/preInt

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
        preInt = preHist.Integral(0,maxBin)
        if preInt: preHist.Scale(dirInt/preInt)
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

def qcdClosureTwoCuts(name, dirSet, treename="tr/simpleTree", preSet=None, additionalSets=[], cut="1", cut2="1", noScale=False, binning=None):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain(treename.replace("tr", "tr_jControl"))

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    weight2 = "weight*({})".format(cut2)
    nBins = range(0,100,10)+[130, 170, 230, 300, 400, 500, 600]
    #range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200, 300, 400, 500, 600]
    if binning: nBins = binning

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data:
        for bin in range(dirHist.FindBin(100),dirHist.GetNbinsX()+2):
            dirHist.SetBinContent(bin, 0)
            dirHist.SetBinError(bin, 0)
    preHist, gjetSyst, fitScale, err, norm = getGJetFitPredictionTwoCuts(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins, weight2)
    totUnc = aux.addHistUncert(preHist, gjetSyst)

    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weight2, nBins)
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
    aux.save("qcdClosure_twoCuts_{}".format(name))

############################## Original ########################################

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
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
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
            signal.GetBinContent(bin), signal.GetBinError(bin), 1.3)

############################ raw met ##########################################

def finalDistributionRawMet(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1"):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain(treename.replace("tr", "tr_jControl"))
    eTree = ROOT.TChain(treename.replace("tr", "tr_eControl"))

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "metRaw"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
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
    aux.save("finalDistributionRawMet_{}".format(name), normal=False)
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
            signal.GetBinContent(bin), signal.GetBinError(bin), 1.3)

##################### Scale simulation by different lumi ######################

def finalDistributionLowerLumi(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1", lumi=36.45):
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
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins)
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weight, nBins, treename.replace("tr", "tr_jControl"))
    mcPreHist.Scale(norm)
    mcPreHist.Scale(lumi/36.459205898)
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        gjetHist.SetBinContent(bin, cOld - subT)
        gjetSyst.SetBinContent(bin, cOld - subT)
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
    for h in zgHist, wgHist, ttgHist, zHist, wHist, ttHist:
        h.Scale(lumi/36.459205898)
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
    signal.Scale(lumi/36.459205898)

    aux.drawOpt(signal, "signal")
    #signal.Add(totStat)
    signal.SetLineColor(ROOT.kAzure)

    # beautify
    aux.drawOpt(dirHist, "data")
    aux.drawOpt(totUnc, "totUnc")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    m.add(dirHist, "Data")
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
    l.lum = ROOT.TLatex( .63, .95, "%.2f fb^{-1} (13 TeV)"%(lumi) )
    l.draw()

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
            signal.GetBinContent(bin), signal.GetBinError(bin), 1.3)

################# ht stacked ###############################

def finalDistributionHTbins(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1"):
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
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

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
        for bin in range(gjetHist_i.GetNbinsX()+2):
            cOld = gjetHist_i.GetBinContent(bin)
            subT = mcPreHist.GetBinContent(bin)
            gjetHist_i.SetBinContent(bin, cOld - subT)
            gjetSyst_i.SetBinContent(bin, cOld - subT)
        if gjetHist: gjetHist.Add(gjetHist_i)
        else: gjetHist = gjetHist_i
        if gjetSyst:
            gjetSyst.Add(gjetSyst_i)
            #for bin in range(gjetSyst.GetNbinsX()+2):
            #    gjetSyst.SetBinContent(bin, gjetSyst.GetBinContent(bin)+gjetSyst_i.GetBinContent(bin))
            #    gjetSyst.SetBinError(bin, gjetSyst.GetBinError(bin)+gjetSyst_i.GetBinError(bin))
        else: gjetSyst = gjetSyst_i
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

############ two Trees ######################

def finalDistributionTwoTrees(name, dirSet, preSet=None, treename1="tr/simpleTree", treename2="tr_jControl/simpleTree", cut="1", lumi=36.459205898):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename1)
    preTree = ROOT.TChain(treename2)
    eTree = ROOT.TChain(treename1)

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename1, preSet, weight, variable, nBins)
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weight, nBins, treename2)
    mcPreHist.Scale(lumi/36.459205898)
    mcPreHist.Scale(norm)
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        gjetHist.SetBinContent(bin, cOld - subT)
        gjetSyst.SetBinContent(bin, cOld - subT)
    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weight, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.createHistoFromTree(eTree, variable, weight, nBins)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    eHist.Scale( 0.0 ) # disable data driven bkg
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    # warning: scaling systematic errors are displayed as statistical ones
    zgHist = aux.createHistoFromDatasetTree(zg, variable, weight, nBins, treename1)
    wgHist = aux.createHistoFromDatasetTree(wg, variable, weight, nBins, treename1)
    ttgHist = aux.createHistoFromDatasetTree(ttg, variable, weight, nBins, treename1)
    zHist = aux.createHistoFromDatasetTree(znunu, variable, weight, nBins, treename1)
    wHist = aux.createHistoFromDatasetTree(wjets, variable, weight, nBins, treename1)
    ttHist = aux.createHistoFromDatasetTree(ttjets_ht, variable, weight, nBins, treename1)
    for h in zgHist, wgHist, ttgHist, zHist, wHist, ttHist:
        h.Scale(lumi/36.459205898)
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

    signal = aux.createHistoFromDatasetTree(t5wg_1600_100, variable, weight, nBins, treename1)
    signal.Scale(lumi/36.459205898)

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
    legHeader += treename1.split("/")[0]
    if cut is not "1": legHeader += aux.beautifyCutString(cut)
    if "_ee" in treename1: legHeader += " EE"
    if "dPhi" in treename1: legHeader += " dPhi"
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
    l.lum = ROOT.TLatex( .63, .95, "%.2f fb^{-1} (13 TeV)"%(lumi) )
    l.draw()
    aux.save("finalDistribution_{}".format(name), normal=False)

################# several eta regions #########################################
def finalDistributionEta(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1"):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain("tr_jControl/simpleTree")
    eTree = ROOT.TChain(treename)

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins)
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weight, nBins, "tr_jControl/simpleTree")
    mcPreHist.Scale(norm)
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        gjetHist.SetBinContent(bin, cOld - subT)
        gjetSyst.SetBinContent(bin, cOld - subT)
    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weight, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.createHistoFromTree(eTree, variable, weight, nBins)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    eHist.Scale( 0.0 ) # disable data driven bkg
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    # warning: scaling systematic errors are displayed as statistical ones
    zgHist = aux.createHistoFromDatasetTree(zg, variable, weight, nBins, treename)
    wgHist = aux.createHistoFromDatasetTree(wg, variable, weight, nBins, treename)
    ttgHist = aux.createHistoFromDatasetTree(ttg, variable, weight, nBins, treename)
    zHist = aux.createHistoFromDatasetTree(znunu, variable, weight, nBins, treename)
    wHist = aux.createHistoFromDatasetTree(wjets, variable, weight, nBins, treename)
    ttHist = aux.createHistoFromDatasetTree(ttjets_ht, variable, weight, nBins, treename)
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
    legHeader += treename.split("/")[0]
    if cut is not "1": legHeader += beautifyCutString(cut)
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


######## scale variation for eg object

def finalDistributionEGscale(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1", cutPre="1", egScale=1.):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain("tr_jControl/simpleTree")
    eTree = ROOT.TChain(treename)

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "met"
    variableScaled = "sqrt(met**2+(pt*{0})**2+2*met*(pt*{0})*cos(dphi))".format(1-egScale)
    weight = "weight*({})".format(cut)
    weightPre = "weight*({})".format(cutPre)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variableScaled, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPredictionEGscaled(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins, weightPre, variableScaled)
    if not gjetHist: return
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weightPre, nBins, "tr_jControl/simpleTree")
    mcPreHist.Scale(norm)
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        gjetHist.SetBinContent(bin, cOld - subT)
        gjetSyst.SetBinContent(bin, cOld - subT)
    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weightPre, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.createHistoFromTree(eTree, variable, weight, nBins)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    eHist.Scale( 0.0 ) # disable data driven bkg
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    # warning: scaling systematic errors are displayed as statistical ones
    zgHist = aux.createHistoFromDatasetTree(zg, variable, weight, nBins, treename)
    wgHist = aux.createHistoFromDatasetTree(wg, variable, weight, nBins, treename)
    ttgHist = aux.createHistoFromDatasetTree(ttg, variable, weight, nBins, treename)
    zHist = aux.createHistoFromDatasetTree(znunu, variable, weight, nBins, treename)
    wHist = aux.createHistoFromDatasetTree(wjets, variable, weight, nBins, treename)
    ttHist = aux.createHistoFromDatasetTree(ttjets_ht, variable, weight, nBins, treename)
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
    legHeader += treename.split("/")[0].replace("tr_", "").replace("tr","")
    if cut is not "1": legHeader += aux.beautifyCutString(cut)
    if cutPre is not "1": legHeader += "%%%%"+aux.beautifyCutString(cutPre)
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
    selEfficiency = dirHist.GetEntries()/dirTree.GetEntries()

    l = aux.Label(sim= not dirSet==data, info=dirSet.label if dirSet != data else "Selected {:.0%}".format(selEfficiency))
    aux.save("finalDistribution_{}".format(name), normal=False)

# cut on prediction and observation seperately

def finalDistributionTwoCuts(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1", cutPre="1"):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain("tr_jControl/simpleTree")
    eTree = ROOT.TChain(treename)

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    weightPre = "weight*({})".format(cutPre)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPredictionTwoCuts(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins, weightPre)
    if not gjetHist: return
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weightPre, nBins, "tr_jControl/simpleTree")
    mcPreHist.Scale(norm)
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        gjetHist.SetBinContent(bin, cOld - subT)
        gjetSyst.SetBinContent(bin, cOld - subT)
    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weightPre, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.createHistoFromTree(eTree, variable, weight, nBins)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    eHist.Scale( 0.0 ) # disable data driven bkg
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    # warning: scaling systematic errors are displayed as statistical ones
    zgHist = aux.createHistoFromDatasetTree(zg, variable, weight, nBins, treename)
    wgHist = aux.createHistoFromDatasetTree(wg, variable, weight, nBins, treename)
    ttgHist = aux.createHistoFromDatasetTree(ttg, variable, weight, nBins, treename)
    zHist = aux.createHistoFromDatasetTree(znunu, variable, weight, nBins, treename)
    wHist = aux.createHistoFromDatasetTree(wjets, variable, weight, nBins, treename)
    ttHist = aux.createHistoFromDatasetTree(ttjets_ht, variable, weight, nBins, treename)
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
    legHeader += treename.split("/")[0].replace("tr_", "").replace("tr","")
    if cut is not "1": legHeader += aux.beautifyCutString(cut)
    if cutPre is not "1": legHeader += "%%%%"+aux.beautifyCutString(cutPre)
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
    selEfficiency = dirHist.GetEntries()/dirTree.GetEntries()

    l = aux.Label(sim= not dirSet==data, info=dirSet.label if dirSet != data else "Selected {:.0%}".format(selEfficiency))
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
            signal.GetBinContent(bin), signal.GetBinError(bin), 1.3)
    #print dc.limit()


def finalDistributionEta(name, dirSet, preSet=None, treename="tr/simpleTree", cut="1"):
    if not preSet: preSet = dirSet
    dirTree = ROOT.TChain(treename)
    preTree = ROOT.TChain("tr_jControl/simpleTree")
    eTree = ROOT.TChain(treename)

    for f in dirSet.files: dirTree.Add(f)
    for f in preSet.files: preTree.Add(f)
    for f in dirSet.files: eTree.Add(f)

    variable = "met"
    weight = "weight*({})".format(cut)
    nBins = range(0,100,10)+range(100,200,10)+[200, 250, 300, 600]
    nBins = range(0,200,10)+[200,250]+range(300,500,20)+[600,700,800,900,910]
    nBins = range(0,200,10)+[200, 250, 300, 350, 400, 500, 600, 700]
    #nBins = range(0,500,10)+[500, 550, 600, 650, 660]

    dirHist = aux.createHistoFromTree(dirTree, variable, weight, nBins)
    if dirSet == data or "data" in name:
        for bin in range(dirHist.GetNbinsX()+2):
            if dirHist.GetBinCenter(bin) > 350 and "_ee" not in name:
                dirHist.SetBinContent(bin,0)
                dirHist.SetBinError(bin,0)

    gjetHist, gjetSyst, fitScale, err, norm = getGJetFitPrediction(dirTree, preTree, name, dirSet, treename, preSet, weight, variable, nBins)
    gjetHist.SetLineColor(ROOT.kCyan)
    # correct for other backgrounds
    mcPreHist = aux.createHistoFromDatasetTree(zg+wg+ttg+wjets+ttjets_ht+znunu, "{}*{}".format(variable,fitScale), weight, nBins, "tr_jControl/simpleTree")
    mcPreHist.Scale(norm)
    for bin in range(gjetHist.GetNbinsX()+2):
        cOld = gjetHist.GetBinContent(bin)
        subT = mcPreHist.GetBinContent(bin)
        gjetHist.SetBinContent(bin, cOld - subT)
        gjetSyst.SetBinContent(bin, cOld - subT)
    gjetHistUnw = aux.createHistoFromTree(preTree, variable, weight, nBins)
    gjetHistUnw.Scale(dirHist.Integral(0,dirHist.FindBin(100)-1)/gjetHistUnw.Integral(0,dirHist.FindBin(100)-1))

    eHist = aux.createHistoFromTree(eTree, variable, weight, nBins)
    eHist.Scale( 0.0267 if dirSet==data else 0.0154 )
    eHist.Scale( 0.0 ) # disable data driven bkg
    eHist.SetLineColor(ROOT.kGreen)
    eSyst = aux.getSysHisto(eHist, 0.3)

    # warning: scaling systematic errors are displayed as statistical ones
    zgHist = aux.createHistoFromDatasetTree(zg, variable, weight, nBins, treename)
    wgHist = aux.createHistoFromDatasetTree(wg, variable, weight, nBins, treename)
    ttgHist = aux.createHistoFromDatasetTree(ttg, variable, weight, nBins, treename)
    zHist = aux.createHistoFromDatasetTree(znunu, variable, weight, nBins, treename)
    wHist = aux.createHistoFromDatasetTree(wjets, variable, weight, nBins, treename)
    ttHist = aux.createHistoFromDatasetTree(ttjets_ht, variable, weight, nBins, treename)
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
    legHeader += treename.split("/")[0]
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
            signal.GetBinContent(bin), signal.GetBinError(bin), 1.3)
    #print dc.limit()

def comparisons(name="n_vertex"):
    h1 = aux.stdHist(data, "tr/"+name)
    h2 = aux.stdHist(dataHt, "tr_jControl/"+name)
    h2.Scale(h1.Integral()/h2.Integral())
    aux.drawOpt(h1, "data")
    aux.drawOpt(h2, "pre")
    h2.SetLineColor(2)
    c = ROOT.TCanvas()
    h1.Draw(h1.drawOption_)
    h2.Draw(h1.drawOption_+"same")
    if name == "n_jet":
        h3 = h2.Clone(aux.randomName())
        for b in range(h3.GetNbinsX()+2):
            h3.SetBinContent(b, h3.GetBinContent(b+1))
        h3.SetLineStyle(2)
        h3.Draw(h1.drawOption_+"same")

    r = ratio.Ratio("#gamma/Jet", h1, h2)
    if name == "n_jet": r.draw(0, 2)
    elif name == "n_vertex": r.draw(.95,1.05)
    else: r.draw()
    aux.save("simpleTest_"+name)

def photonPtAnalysis(norm=True, treename="tr/simpleTree"):
    name = "photonPtAnalysis_"
    dirname = treename.split("/")[0]
    if dirname != "tr": name += dirname[3:] + "_"
    if norm: name += "norm_"
    t = ROOT.TChain(treename)
    for f in data.files: t.Add(f)
    tSig = ROOT.TChain(treename)
    for f in t5wg_1600_100.files: tSig.Add(f)

    ptBinning = range(100,801,10)

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    for icut, cut in enumerate(["1", "emht>800", "emht>900", "emht>1000", "emht>1100", "emht>1200"]):
        h = aux.createHistoFromTree(t, "pt", cut, ptBinning)
        h.SetXTitle("p_{T} (GeV)")
        h.SetLineColor(icut+1)
        h.drawOption_="hist"
        if norm: h.Scale(1./h.Integral())
        m.add(h, cut)
    m.Draw()
    aux.save(name+"EMHTcuts")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    for icut, cut in enumerate(["1", "ht>700", "pt/emht<.2", ".2<pt/emht"]):
        h = aux.createHistoFromTree(t, "pt", cut, ptBinning)
        h.SetXTitle("p_{T} (GeV)")
        h.SetLineColor(icut+1)
        h.drawOption_="hist"
        if norm: h.Scale(1./h.Integral())
        m.add(h, cut)
    m.Draw()
    aux.save(name+"OtherCuts")

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    for icut, cut in enumerate(["1", "ht>700"]):
        h = aux.createHistoFromTree(t, "pt/emht", cut, aux.frange(0,0.8,0.01))
        h.SetXTitle("p_{T}/EMH_{T}")
        h.SetLineColor(icut+1)
        h.drawOption_="hist"
        if norm: h.Scale(1./h.Integral())
        m.add(h, cut)
    h = aux.createHistoFromTree(tSig, "pt/emht", "1", aux.frange(0,0.8,0.01))
    h.SetLineColor(ROOT.kMagenta)
    h.SetXTitle("p_{T}/EMH_{T}")
    h.drawOption_="hist"
    if norm: h.Scale(1./h.Integral())
    m.add(h, "Signal")
    m.Draw()
    aux.save(name+"PtOverEMHT")

def correlation2d(x="pt", xBins=(100,100,200), y="emht", yBins=(130,700,2000), treeDir="tr", cut="1", drawMiddle=False):
    name = "correlation2d_{}_{}_vs_{}".format(treeDir,y,x)
    tree = ROOT.TChain(treeDir+"/simpleTree")
    for f in data.files: tree.Add(f)
    hname = aux.randomName()
    result = ROOT.TH2F(hname, ";{};{};Events".format(x,y), *list(xBins+yBins))
    tree.Draw("{}:{}>>{}".format(y,x,hname), cut, "goff")
    style.style2d()
    c = ROOT.TCanvas()
    result.Draw("colz")
    if drawMiddle:
        l = ROOT.TLine()
        l.SetLineColor(ROOT.kBlue)
        if x=="pt" and y=="emht":
            l.DrawLine(350,700,1000,2000)
        else:
            l.DrawLine(result.GetXaxis().GetXmin(), result.GetYaxis().GetXmin(), result.GetXaxis().GetXmax(), result.GetYaxis().GetXmax())
    aux.save(name)
    style.defaultStyle()

def studyCorrelation(dataset, x, xName, xBins, y, yName, yBins, dirname="tr", cut="1"):
    name = aux.randomName()
    h2 = ROOT.TH2F(name, ";{};{};a.u.".format(xName,yName), len(xBins)-1, array.array("d", xBins), len(yBins)-1, array.array("d", yBins))
    h2.Sumw2()
    tree = ROOT.TChain(dirname+"/simpleTree")
    for f in dataset.files: tree.Add(f)
    tree.Draw("{}:{}>>{}".format(y,x,name), "weight*{}".format(cut), "goff")
    #aux.appendFlowBin2d(h2)

    m = multiplot.Multiplot()
    ROOT.gStyle.SetPalette(55)
    for h in aux.getProjections(h2, scale=True):
        h.Scale(1., "width")
        m.add(h, h.GetName())
    c = ROOT.TCanvas()
    m.minimum = m.getMinimum()/2
    m.maximum = m.getMaximum()*1.5
    m.Draw()
    l = aux.Label()
    sName = "correlationStudy_{}_{}_{}".format(dirname, x, y)
    aux.save(sName, log=False)
    c.SetLogy()
    ROOT.gPad.SaveAs("plots/{}_log.pdf".format(sName))
    c.SetLogy(False)
    px = h2.ProfileX()
    px.Draw()
    aux.save(sName+"_profileX", log=False)
    py = h2.ProfileY()
    py.Draw()
    aux.save(sName+"_profileY", log=False)

def studyCorrelations():
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "abs(eta)", "|#eta|", aux.frange(0,2.6,.25), "tr_jControl")
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "abs(phi)", "|#eta|", aux.frange(0,2.6,.25), "tr_jControl")
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "abs(eta)", "|#eta|", aux.frange(0,1.6,.25))
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "abs(eta)", "|#eta|", aux.frange(1.5,2.6,.25), "tr_ee")
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "pt", "p_{T}", [100,150,200,300,400,500,600])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "pt", "p_{T}", [100,150,200,300,400,500,600], "tr_ee")
    studyCorrelation(data, "abs(eta)", "|#eta|", [0,.5,1,1.5], "pt", "p_{T}", [100,150,200,300,400,500,600])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "sie", "#sigma_{i#etai#eta}", [0, 0.008, 0.0087, 0.009, 0.0095, 0.01, 0.0103])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "sie", "#sigma_{i#etai#eta}", [0, 0.02, 0.022, 0.023, 0.024, 0.025, 0.026, 0.0275], "tr_ee")
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "sip", "#sigma_{i#phii#phi}", [0, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.02])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "sip", "#sigma_{i#phii#phi}", [0, 0.025, 0.025, 0.0275, 0.03, 0.04, 0.05], "tr_ee")
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "cIso", "I_{#pm}", [0, 1e-6, 1, 2, 3])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "nIso", "I_{0}", [0, 1e-6, 1, 2, 4, 6])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "pIso", "I_{#gamma}", [0, 1e-6, 1, 2, 3])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "cIsoWorst", "worst I_{#pm}", [0, 1e-6, 20, 50, 100, 200, 500])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "hoe", "H/E", [0, 1e-6, 0.003, 0.005 ])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "r9", "r9", [0, .8, .9, .925, .95, 1])
    studyCorrelation(data, "met", "#it{E}_{T}^{miss} (GeV)", range(0,200,10)+[200, 250, 300, 350], "vtx", "vertex multiplicity", [10, 15, 20, 25, 30])
    studyCorrelation(data, "pt", "#it{p}_{T} (GeV)", range(100,400,10)+[400, 500, 600, 700, 800, 900, 1000], "sie", "#sigma_{i#etai#eta}", [0, 0.008, 0.0087, 0.009, 0.0095, 0.01, 0.0103])
    studyCorrelation(data, "njet", "jet multiplicity", range(1,10), "sie", "#sigma_{i#etai#eta}", [0, 0.008, 0.0087, 0.009, 0.0095, 0.01, 0.0103])
    studyCorrelation(data, "sie", "#sigma_{i#etai#eta}", aux.frange(0.007, 0.0106, 0.0001), "pt", "p_{T}", [100,150,200,300,400,500,600])



def checkMetCorrelation(var, dname, dataset, dirname="tr", cut="1"):
    import array
    metBins = range(0,200,10)+[200, 250, 300, 350]
    if var == "sie":
        varBins = [0, 0.008, 0.0087, 0.009, 0.0095, 0.01, 0.0103]
        vName = "#sigma_{i#etai#eta}"
        if "ee" in dirname:
            varBins = [0, 0.02, 0.022, 0.023, 0.024, 0.025, 0.026, 0.0275]
    elif var == "sip":
        varBins = [0, 0.009, 0.01, 0.011, 0.012, 0.013, 0.014, 0.015, 0.02]
        vName = "#sigma_{i#phii#phi}"
        if "ee" in dirname:
            varBins = [0, 0.025, 0.025, 0.0275, 0.03, 0.04, 0.05]
    elif var == "cIso":
        varBins = [0, 1e-6, 1, 2, 3]
        vName = "I_{#pm}"
    elif var == "nIso":
        varBins = [0, 1e-6, 1, 2, 4, 6]
        vName = "I_{0}"
    elif var == "pIso":
        varBins = [0, 1e-6, 1, 2, 3]
        vName = "I_{#gamma}"
    elif var == "r9":
        varBins = [0, .8, .9, .925, .95, 1]
        vName = "r9"
    elif var == "hoe":
        varBins = [0, 1e-6, 0.003, 0.005 ]
        vName = "H/E"
    elif var == "cIsoWorst":
        varBins = [0, 1e-6, 20, 50, 100, 200, 500 ]
        vName = "worst I_{#pm}"


    name = aux.randomName()
    h2 = ROOT.TH2F(name, ";#it{{E}}_{{T}}^{{miss}} (GeV);{};a.u.".format(vName), len(metBins)-1, array.array("d", metBins), len(varBins)-1, array.array("d", varBins))
    h2.Sumw2()
    tree = ROOT.TChain(dirname+"/simpleTree")
    for f in dataset.files: tree.Add(f)
    tree.Draw("{}:met>>{}".format(var,name), "weight*{}".format(cut), "goff")
    aux.appendFlowBin2d(h2)

    m = multiplot.Multiplot()
    ROOT.gStyle.SetPalette(55)
    for h in aux.getProjections(h2, scale=True):
        h.Scale(1.,"width")
        m.add(h, h.GetName())
    c = ROOT.TCanvas()
    m.minimum = 1e-5
    m.Draw()
    l = aux.Label(info=dname)
    sName = "metCorrolation_{}_{}_{}".format(dname,dirname, var)
    aux.save(sName, log=False)
    c.SetLogy()
    ROOT.gPad.SaveAs("plots/{}_log.pdf".format(sName))

def runAreas():
    dataB= Dataset("SinglePhoton_Run2016B-23Sep2016-v3", 0, ROOT.kBlack )
    dataB.label = "RunB"
    dataHtB = Dataset("JetHT_Run2016B-23Sep2016-v3", 0, ROOT.kBlack )
    finalDistributionLowerLumi("dataB", dataB, dataHtB, lumi=5.93369235121)
    finalDistributionLowerLumi("dataB_2000emht", dataB, dataHtB, cut="2000<emht", lumi=5.93369235121)
    finalDistributionLowerLumi("dataB_emht1000", dataB, dataHtB, cut="emht<1000", lumi=5.93369235121)

    dataC= Dataset("SinglePhoton_Run2016C-23Sep2016-v1", 0, ROOT.kBlack )
    dataC.label = "RunC"
    dataHtC = Dataset("JetHT_Run2016C-23Sep2016-v1", 0, ROOT.kBlack )
    finalDistributionLowerLumi("dataC", dataC, dataHtC, lumi=2.64596808309)

    dataD= Dataset("SinglePhoton_Run2016D-23Sep2016-v1", 0, ROOT.kBlack )
    dataD.label = "RunD"
    dataHtD = Dataset("JetHT_Run2016D-23Sep2016-v1", 0, ROOT.kBlack )
    finalDistributionLowerLumi("dataD", dataD, dataHtD, lumi=4.35344881055)

    dataE= Dataset("SinglePhoton_Run2016E-23Sep2016-v1", 0, ROOT.kBlack )
    dataE.label = "RunE"
    dataHtE = Dataset("JetHT_Run2016E-23Sep2016-v1", 0, ROOT.kBlack )
    finalDistributionLowerLumi("dataE", dataE, dataHtE, lumi=4.04973213414)

    dataF= Dataset("SinglePhoton_Run2016F-23Sep2016-v1", 0, ROOT.kBlack )
    dataF.label = "RunF"
    dataHtF = Dataset("JetHT_Run2016F-23Sep2016-v1", 0, ROOT.kBlack )
    finalDistributionLowerLumi("dataF", dataF, dataHtF, lumi=3.16008842022)

    dataG= Dataset("SinglePhoton_Run2016G-23Sep2016-v1", 0, ROOT.kBlack )
    dataG.label = "RunG"
    dataHtG = Dataset("JetHT_Run2016G-23Sep2016-v1", 0, ROOT.kBlack )
    finalDistributionLowerLumi("dataG", dataG, dataHtG, lumi=7.55445363513)

    dataH= Dataset("SinglePhoton_Run2016H-PromptReco-v2", 0, ROOT.kBlack )
    dataH.label = "RunH-v2"
    dataHtH = Dataset("JetHT_Run2016H-PromptReco-v2", 0, ROOT.kBlack )
    finalDistributionLowerLumi("dataH", dataH, dataHtH, lumi=8.54503959046)

def ht800Trigger():
    dataBG = Dataset("SinglePhoton_Run2016B-23Sep2016-v3", 0, ROOT.kBlack ) \
        + Dataset("SinglePhoton_Run2016C-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("SinglePhoton_Run2016D-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("SinglePhoton_Run2016E-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("SinglePhoton_Run2016F-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("SinglePhoton_Run2016G-23Sep2016-v1", 0, ROOT.kBlack )
    dataBG.label = "RunB-G"
    dataHtBG = Dataset("JetHT_Run2016B-23Sep2016-v3", 0, ROOT.kBlack ) \
        + Dataset("JetHT_Run2016C-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("JetHT_Run2016D-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("JetHT_Run2016E-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("JetHT_Run2016F-23Sep2016-v1", 0, ROOT.kBlack ) \
        + Dataset("JetHT_Run2016G-23Sep2016-v1", 0, ROOT.kBlack )
    finalDistributionTwoTrees("data_ht800trigger", dataHtBG, dataHtBG, "tr_ht900/simpleTree", "tr_jControl_ht900/simpleTree", "1", 27.697383434340001)
    #finalDistributionLowerLumi("dataBG_900emht", dataBG, dataHtBG, cut="emht>900", lumi=27.6973)






if __name__ == "__main__":
    allMC = gjets+qcd+zg+wg+ttg+wjets+ttjets_ht+znunu
    allMC.label = "MC mix"
    #finalDistributionEta("mc_eta_0", allMC, allMC, "tr/simpleTree")
    #finalDistributionEta("mc_eta_1", allMC, allMC, "tr_eta1/simpleTree")
    #finalDistributionEta("mc_eta_2", allMC, allMC, "tr_eta2/simpleTree")
    #finalDistributionEta("mc_eta_3", allMC, allMC, "tr_eta3/simpleTree")
    #finalDistributionEta("mc_eta_4", allMC, allMC, "tr_eta4/simpleTree")
    #finalDistributionEta("mc_eta_5", allMC, allMC, "tr_eta5/simpleTree")

    #finalDistributionEta("eta_all", data, dataHt, "tr/simpleTree")
    #finalDistributionEta("eta_1", data, dataHt, "tr_eta1/simpleTree")
    #finalDistributionEta("eta_2", data, dataHt, "tr_eta2/simpleTree")
    #finalDistributionEta("eta_3", data, dataHt, "tr_eta3/simpleTree")
    #finalDistributionEta("eta_4", data, dataHt, "tr_eta4/simpleTree")
    #finalDistributionEta("eta_5", data, dataHt, "tr_eta5/simpleTree")
    #finalDistributionEta("tight_eb", data, dataHt, "tr_tight/simpleTree")
    #finalDistributionEta("tight_ee", data, dataHt, "tr_ee_tight/simpleTree")
    #finalDistribution("data", data, dataHt)
    #finalDistributionMediumPU("mediumPU", data, dataHt)
    #finalDistributionHTbins("htbins", data, dataHt)
    #finalDistributionEta("tight", data, dataHt, "tr_tight/simpleTree")
    #finalDistributionTwoCuts("onlyPhotonCut_100pt", data, dataHt, cut="100<pt&&pt<120")
    #finalDistributionTwoCuts("onlyPhotonCut_120pt", data, dataHt, cut="120<pt&&pt<140")
    #finalDistributionTwoCuts("onlyPhotonCut_140pt", data, dataHt, cut="140<pt&&pt<160")
    #finalDistributionTwoCuts("onlyPhotonCut_160pt", data, dataHt, cut="160<pt&&pt<180")
    #finalDistributionTwoCuts("onlyPhotonCut_180pt", data, dataHt, cut="180<pt&&pt<200")
    #finalDistributionTwoCuts("onlyPhotonCut_200pt", data, dataHt, cut="200<pt&&pt<220")
    #finalDistributionTwoCuts("onlyPhotonCut_220pt", data, dataHt, cut="220<pt&&pt<240")
    #finalDistributionTwoCuts("onlyPhotonCut_240pt", data, dataHt, cut="240<pt&&pt<260")
    #finalDistributionTwoCuts("onlyPhotonCut_260pt", data, dataHt, cut="260<pt&&pt<280")
    #finalDistributionTwoCuts("onlyPhotonCut_280pt", data, dataHt, cut="280<pt&&pt<300")
    #finalDistributionTwoCuts("onlyPhotonCut_300pt", data, dataHt, cut="300<pt&&pt<320")
    #finalDistributionTwoCuts("onlyPhotonCut_320pt", data, dataHt, cut="320<pt&&pt<340")
    #finalDistributionTwoCuts("onlyPhotonCut_340pt", data, dataHt, cut="340<pt&&pt<360")
    #finalDistributionTwoCuts("onlyPhotonCut_360pt", data, dataHt, cut="360<pt&&pt<380")
    #finalDistributionTwoCuts("onlyPhotonCut_380pt", data, dataHt, cut="380<pt&&pt<400")
    #finalDistributionTwoCuts("onlyPhotonCut_400pt", data, dataHt, cut="400<pt&&pt<420")
    #finalDistributionTwoCuts("onlyPhotonCut_420pt", data, dataHt, cut="420<pt&&pt<440")
    #finalDistributionTwoCuts("onlyPhotonCut_440pt", data, dataHt, cut="440<pt&&pt<460")
    #finalDistributionTwoCuts("onlyPhotonCut_460pt", data, dataHt, cut="460<pt&&pt<480")
    #finalDistributionTwoCuts("onlyPhotonCut_480pt", data, dataHt, cut="480<pt&&pt<500")
    #finalDistributionTwoCuts("onlyPhotonCut_500pt", data, dataHt, cut="500<pt&&pt<520")
    #finalDistributionTwoCuts("onlyPhotonCut_520pt", data, dataHt, cut="520<pt&&pt<540")
    #finalDistributionTwoCuts("onlyPhotonCut_540pt", data, dataHt, cut="540<pt&&pt<560")
    #finalDistributionTwoCuts("onlyPhotonCut_560pt", data, dataHt, cut="560<pt&&pt<580")
    #finalDistributionTwoCuts("onlyPhotonCut_580pt", data, dataHt, cut="580<pt&&pt<600")
    #finalDistributionTwoCuts("onlyPhotonCut_600pt", data, dataHt, cut="600<pt&&pt<620")
    #finalDistributionTwoCuts("onlyPhotonCut_620ppppt", data, dataHt, cut="620<pt")
    #finalDistributionTwoCuts("onlyPhotonCut_ht700", data, dataHt, cut="700<ht")
    #finalDistributionTwoCuts("onlyPhotonCut_ht800", data, dataHt, cut="800<ht")
    #finalDistributionTwoCuts("onlyPhotonCut_ht1300", data, dataHt, cut="1300<ht")
    #finalDistributionTwoCuts("onlyPhotonCut_ht1000", data, dataHt, cut="1000<ht")
    #finalDistributionTwoCuts("onlyPhotonCut_ptOverEmht0p2", data, dataHt, cut="pt/emht<0.2")
    #finalDistributionTwoCuts("onlyPhotonCut_0p2ptOverEmht", data, dataHt, cut="0.2<pt/emht")
    #finalDistributionTwoCuts("onlyPhotonCut_ee_ptOverEmht0p2", data, dataHt, "tr_ee/simpleTree", cut="pt/emht<0.2")
    #finalDistributionTwoCuts("onlyPhotonCut_ee_0p2ptOverEmht", data, dataHt, "tr_ee/simpleTree", cut="0.2<pt/emht")
    #finalDistributionTwoCuts("onlyPhotonCut_ptOverEmht0p2_eta1", data, dataHt, cut="pt/emht<0.2&&abs(eta)<.5")
    #finalDistributionTwoCuts("onlyPhotonCut_ptOverEmht0p2_eta3", data, dataHt, cut="pt/emht<0.2&&1<abs(eta)")
    #finalDistributionTwoCuts("onlyPhotonCut_0p2ptOverEmht_eta1", data, dataHt, cut="0.2<pt/emht&&abs(eta)<.5")
    #finalDistributionTwoCuts("onlyPhotonCut_0p2ptOverEmht_eta3", data, dataHt, cut="0.2<pt/emht&&1<abs(eta)")
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
    #finalDistributionEta("data_noLep", data, dataHt, "tr_noLep/simpleTree")

    #comparisons()
    #comparisons("n_jet")
    #runAreas()

    #photonPtAnalysis()
    #photonPtAnalysis(False)
    #photonPtAnalysis(treename="tr_eta1/simpleTree")
    #photonPtAnalysis(treename="tr_eta2/simpleTree")
    #photonPtAnalysis(treename="tr_eta3/simpleTree")
    #photonPtAnalysis(treename="tr_eta4/simpleTree")
    #photonPtAnalysis(treename="tr_eta5/simpleTree")

    #finalDistributionTwoCuts("veryTight", data, dataHt, "tr_veryTight/simpleTree")
    #finalDistributionTwoCuts("veryTight_eta1", data, dataHt, "tr_veryTight/simpleTree", cut="abs(eta)<.5")
    #finalDistributionTwoCuts("veryTight_eta2", data, dataHt, "tr_veryTight/simpleTree", cut=".5<abs(eta)&&abs(eta)<1")
    #finalDistributionTwoCuts("veryTight_eta3", data, dataHt, "tr_veryTight/simpleTree", cut="1<abs(eta)&&abs(eta)<1.5")

    #finalDistributionTwoCuts("onlyPhotonCut_eta1", data, dataHt, cut="abs(eta)<.25")
    #finalDistributionTwoCuts("onlyPhotonCut_eta2", data, dataHt, cut="0.25<abs(eta)&&abs(eta)<.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta3", data, dataHt, cut="0.5<abs(eta)&&abs(eta)<.75")
    #finalDistributionTwoCuts("onlyPhotonCut_eta4", data, dataHt, cut="0.75<abs(eta)&&abs(eta)<1")
    #finalDistributionTwoCuts("onlyPhotonCut_eta5", data, dataHt, cut="1<abs(eta)&&abs(eta)<1.25")
    #finalDistributionTwoCuts("onlyPhotonCut_eta6", data, dataHt, cut="1.25<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta7", data, dataHt, "tr_ee/simpleTree", cut="1.5<abs(eta)&&abs(eta)<1.75")
    #finalDistributionTwoCuts("onlyPhotonCut_eta8", data, dataHt, "tr_ee/simpleTree", cut="1.75<abs(eta)&&abs(eta)<2")
    #finalDistributionTwoCuts("onlyPhotonCut_eta9", data, dataHt, "tr_ee/simpleTree", cut="2<abs(eta)&&abs(eta)<2.25")
    #finalDistributionTwoCuts("onlyPhotonCut_eta10", data, dataHt, "tr_ee/simpleTree", cut="2.25<abs(eta)&&abs(eta)<2.5")

    #finalDistributionTwoCuts("onlyPhotonCut_onsided_etaPlus", data, dataHt, cut="eta>0")
    #finalDistributionTwoCuts("onlyPhotonCut_onsided_etaMinus", data, dataHt, cut="eta<0")
    #finalDistributionTwoCuts("onlyPhotonCut_phiPlus", data, dataHt, cut="phi<0")
    #finalDistributionTwoCuts("onlyPhotonCut_phiMinus", data, dataHt, cut="phi>0")
    #finalDistributionTwoCuts("onlyPhotonCut_phi0", data, dataHt, cut="0.25<phi&&phi<0.75")
    #finalDistributionTwoCuts("onlyPhotonCut_phi1", data, dataHt, cut="-5<phi&&phi<-2")
    #finalDistributionTwoCuts("onlyPhotonCut_phi2", data, dataHt, cut="-2<phi&&phi<-1")
    #finalDistributionTwoCuts("onlyPhotonCut_phi3", data, dataHt, cut="-1<phi&&phi<0")
    #finalDistributionTwoCuts("onlyPhotonCut_phi4", data, dataHt, cut="-0<phi&&phi<1")
    #finalDistributionTwoCuts("onlyPhotonCut_phi5", data, dataHt, cut="1<phi&&phi<2")
    #finalDistributionTwoCuts("onlyPhotonCut_phi6", data, dataHt, cut="2<phi&&phi<5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta1phi1", data, dataHt, cut="-5<phi&&phi<-1 && -1.5<eta&&eta<-0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta1phi2", data, dataHt, cut="-1<phi&&phi<1 && -1.5<eta&&eta<-0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta1phi3", data, dataHt, cut=" 1<phi&&phi<5 && -1.5<eta&&eta<-0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta2phi1", data, dataHt, cut="-5<phi&&phi<-1 && -.5<eta&&eta<.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta2phi2", data, dataHt, cut="-1<phi&&phi<1 && -.5<eta&&eta<.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta2phi3", data, dataHt, cut=" 1<phi&&phi<5 && -.5<eta&&eta<.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta3phi1", data, dataHt, cut="-5<phi&&phi<-1 && .5<eta&&eta<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta3phi2", data, dataHt, cut="-1<phi&&phi<1 && .5<eta&&eta<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_eta3phi3", data, dataHt, cut="1<phi&&phi<5 && .5<eta&&eta<1.5")

    #for njet in range(1,9): finalDistributionTwoCuts("onlyPhotonCut_{}jet".format(njet), data, dataHt, cut="njet=={}".format(njet))
    #for njet in range(1,9): finalDistributionTwoCuts("data_{}jetBoth".format(njet), data, dataHt, cut="njet=={}".format(njet), cutPre="njet=={}".format(njet))
    #for njet in range(1,9): finalDistributionTwoCuts("data_{}jetMinusPred".format(njet), data, dataHt, cut="njet=={}".format(njet), cutPre="njet=={}".format(njet-1))
    #finalDistributionTwoCuts("onlyPhotonCut_dphiPar", data, dataHt, cut="abs(dphi)<0.785")
    #finalDistributionTwoCuts("onlyPhotonCut_dphiR", data, dataHt, cut="0.758<dphi&&dphi<2.356")
    #finalDistributionTwoCuts("onlyPhotonCut_dphiL", data, dataHt, cut="-2.356<dphi&&dphi<-0.758")
    #finalDistributionTwoCuts("onlyPhotonCut_dphiAnti", data, dataHt, cut="2.356<abs(dphi)")

    #finalDistributionTwoCuts("onlyPhotonCut_sie1", data, dataHt, cut="sie<0.008")
    #finalDistributionTwoCuts("onlyPhotonCut_sie2", data, dataHt, cut="0.008<sie&&sie<0.009")
    #finalDistributionTwoCuts("onlyPhotonCut_sie3", data, dataHt, cut="0.009<sie&&sie<0.0095")
    #finalDistributionTwoCuts("onlyPhotonCut_sie4", data, dataHt, cut="0.0095<sie&&sie<0.01")
    #finalDistributionTwoCuts("onlyPhotonCut_sie5", data, dataHt, cut="0.01<sie")

    #finalDistributionTwoCuts("onlyPhotonCut_sie1_eta1", data, dataHt, cut="sie<0.008&&abs(eta)<0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie2_eta1", data, dataHt, cut="0.008<sie&&sie<0.009&&abs(eta)<0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie3_eta1", data, dataHt, cut="0.009<sie&&sie<0.0095&&abs(eta)<0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie4_eta1", data, dataHt, cut="0.0095<sie&&sie<0.01&&abs(eta)<0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie5_eta1", data, dataHt, cut="0.01<sie&&abs(eta)<0.5")

    #finalDistributionTwoCuts("onlyPhotonCut_sie1_eta2", data, dataHt, cut="sie<0.008&&0.5<abs(eta)&&abs(eta)<1.0")
    #finalDistributionTwoCuts("onlyPhotonCut_sie2_eta2", data, dataHt, cut="0.008<sie&&sie<0.009&&0.5<abs(eta)&&abs(eta)<1.0")
    #finalDistributionTwoCuts("onlyPhotonCut_sie3_eta2", data, dataHt, cut="0.009<sie&&sie<0.0095&&0.5<abs(eta)&&abs(eta)<1.0")
    #finalDistributionTwoCuts("onlyPhotonCut_sie4_eta2", data, dataHt, cut="0.0095<sie&&sie<0.01&&0.5<abs(eta)&&abs(eta)<1.0")
    #finalDistributionTwoCuts("onlyPhotonCut_sie5_eta2", data, dataHt, cut="0.01<sie&&0.5<abs(eta)&&abs(eta)<1.0")

    #finalDistributionTwoCuts("onlyPhotonCut_sie1_eta3", data, dataHt, cut="sie<0.008&&1<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie2_eta3", data, dataHt, cut="0.008<sie&&sie<0.009&&1<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie3_eta3", data, dataHt, cut="0.009<sie&&sie<0.0095&&1<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie4_eta3", data, dataHt, cut="0.0095<sie&&sie<0.01&&1<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie5_eta3", data, dataHt, cut="0.01<sie&&1<abs(eta)&&abs(eta)<1.5")

    #checkMetCorrelation("sie", "data", data)
    #checkMetCorrelation("sip", "data", data)
    #checkMetCorrelation("cIso", "data", data)
    #checkMetCorrelation("nIso", "data", data)
    #checkMetCorrelation("pIso", "data", data)
    #checkMetCorrelation("r9", "data", data)
    #checkMetCorrelation("hoe", "data", data)
    #checkMetCorrelation("cIsoWorst", "data", data)
    #checkMetCorrelation("sie", "data", data, "tr_ee")
    #checkMetCorrelation("sip", "data", data, "tr_ee")
    #checkMetCorrelation("cIso", "data", data, "tr_ee")
    #checkMetCorrelation("nIso", "data", data, "tr_ee")
    #checkMetCorrelation("pIso", "data", data, "tr_ee")
    #checkMetCorrelation("r9", "data", data, "tr_ee")
    #checkMetCorrelation("hoe", "data", data, "tr_ee")
    #checkMetCorrelation("cIsoWorst", "data", data, "tr_ee")

    #finalDistributionTwoCuts("onlyPhotonCut_veryTightCopy", data, dataHt, cut="0.006<sie&&sie<0.0095&&0.005<sip&&sip<0.012&&hoe<0.005&&r9>.8&&cIso<0.005&&nIso<0.97&&pIso<0.5&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNoSie", data, dataHt, cut="0.005<sip&&sip<0.012&&hoe<0.005&&r9>.8&&cIso<0.005&&nIso<0.97&&pIso<0.5&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNoSip", data, dataHt, cut="0.006<sie&&sie<0.0095&&hoe<0.005&&r9>.8&&cIso<0.005&&nIso<0.97&&pIso<0.5&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNoHoe", data, dataHt, cut="0.006<sie&&sie<0.0095&&0.005<sip&&sip<0.012&&r9>.8&&cIso<0.005&&nIso<0.97&&pIso<0.5&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNor9", data, dataHt, cut="0.006<sie&&sie<0.0095&&0.005<sip&&sip<0.012&&hoe<0.005&&cIso<0.005&&nIso<0.97&&pIso<0.5&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNocIso", data, dataHt, cut="0.006<sie&&sie<0.0095&&0.005<sip&&sip<0.012&&hoe<0.005&&r9>.8&&nIso<0.97&&pIso<0.5&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNonIso", data, dataHt, cut="0.006<sie&&sie<0.0095&&0.005<sip&&sip<0.012&&hoe<0.005&&r9>.8&&cIso<0.005&&pIso<0.5&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNopIso", data, dataHt, cut="0.006<sie&&sie<0.0095&&0.005<sip&&sip<0.012&&hoe<0.005&&r9>.8&&cIso<0.005&&nIso<0.97&&cIsoWorst<50")
    #finalDistributionTwoCuts("onlyPhotonCut_veryTightNocIsoWorst", data, dataHt, cut="0.006<sie&&sie<0.0095&&0.005<sip&&sip<0.012&&hoe<0.005&&r9>.8&&cIso<0.005&&nIso<0.97&&pIso<0.5")
    #finalDistributionTwoCuts("onlyPhotonCut_noCut", data, dataHt, cut="1")
    #finalDistribution("data_sieCut", data, dataHt, cut="(0.006<sie&&sie<0.0095)||sie<1e-6")

    #ht800Trigger()
    #finalDistribution("data_vtx10", data, dataHt, cut="vtx<10")
    #finalDistribution("data_10vtx15", data, dataHt, cut="10<vtx&&vtx<15")
    #finalDistribution("data_15vtx20", data, dataHt, cut="15<vtx&&vtx<20")
    #finalDistribution("data_20vtx25", data, dataHt, cut="20<vtx&&vtx<25")
    #finalDistribution("data_25vtx30", data, dataHt, cut="25<vtx&&vtx<30")
    #finalDistribution("data_25vtx", data, dataHt, cut="25<vtx")
    #finalDistributionRawMet("data", data, dataHt)
    studyCorrelations()
    #qcdClosure("gqcd", gjets+qcd)
    #qcdClosure("qcd", qcd)
    #finalDistributionTwoCuts("onlyJetCut_1jet_eta1", data, dataHt, cut="1", cutPre="abs(eta)<0.5")
    #finalDistributionTwoCuts("onlyJetCut_1jet_eta2", data, dataHt, cut="1", cutPre="0.5<abs(eta)&&abs(eta)<1")
    #finalDistributionTwoCuts("onlyJetCut_1jet_eta3", data, dataHt, cut="1", cutPre="1<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("onlyJetCut_2jet_eta1", data, dataHt, cut="1", cutPre="abs(phi)<0.5")
    #finalDistributionTwoCuts("onlyJetCut_2jet_eta2", data, dataHt, cut="1", cutPre="0.5<abs(phi)&&abs(phi)<1")
    #finalDistributionTwoCuts("onlyJetCut_2jet_eta3", data, dataHt, cut="1", cutPre="1<abs(phi)&&abs(phi)<1.5")

    #finalDistributionTwoCuts("bothEtaCuts_1jet_eta1", data, dataHt, cut="abs(eta)<.5", cutPre="abs(eta)<0.5")
    #finalDistributionTwoCuts("bothEtaCuts_1jet_eta2", data, dataHt, cut=".5<abs(eta)&&abs(eta)<1", cutPre="0.5<abs(eta)&&abs(eta)<1")
    #finalDistributionTwoCuts("bothEtaCuts_1jet_eta3", data, dataHt, cut="1<abs(eta)&&abs(eta)<1.5", cutPre="1<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("bothEtaCuts_2jet_eta1", data, dataHt, cut="abs(eta)<.5", cutPre="abs(phi)<0.5")
    #finalDistributionTwoCuts("bothEtaCuts_2jet_eta2", data, dataHt, cut=".5<abs(eta)&&abs(eta)<1", cutPre="0.5<abs(phi)&&abs(phi)<1")
    #finalDistributionTwoCuts("bothEtaCuts_2jet_eta3", data, dataHt, cut="1<abs(eta)&&abs(eta)<1.5", cutPre="1<abs(phi)&&abs(phi)<1.5")

    #finalDistributionTwoCuts("photonEtaCuts_1jet_eta1", data, dataHt, cut="abs(eta)<.5", cutPre="1")
    #finalDistributionTwoCuts("photonEtaCuts_1jet_eta2", data, dataHt, cut=".5<abs(eta)&&abs(eta)<1", cutPre="1")
    #finalDistributionTwoCuts("photonEtaCuts_1jet_eta3", data, dataHt, cut="1<abs(eta)&&abs(eta)<1.5", cutPre="1")

    #finalDistributionTwoCuts("onlyPhotonCut_sieAll_eta3", data, dataHt, cut="1<abs(eta)&&abs(eta)<1.5")
    #finalDistributionTwoCuts("onlyPhotonCut_sie1to4_eta3", data, dataHt, cut="sie<.009&&1<abs(eta)&&abs(eta)<1.5")


    #qcdClosure("mc_onlyPhotonCut_100pt", gjets+qcd, cut="100<pt&&pt<120")
    #qcdClosure("mc_onlyPhotonCut_120pt", gjets+qcd, cut="120<pt&&pt<140")
    #qcdClosure("mc_onlyPhotonCut_140pt", gjets+qcd, cut="140<pt&&pt<160")
    #qcdClosure("mc_onlyPhotonCut_160pt", gjets+qcd, cut="160<pt&&pt<180")
    #qcdClosure("mc_onlyPhotonCut_180pt", gjets+qcd, cut="180<pt&&pt<200")
    #qcdClosure("mc_onlyPhotonCut_200pt", gjets+qcd, cut="200<pt&&pt<220")
    #qcdClosure("mc_onlyPhotonCut_220pt", gjets+qcd, cut="220<pt&&pt<240")
    #qcdClosure("mc_onlyPhotonCut_240pt", gjets+qcd, cut="240<pt&&pt<260")
    #qcdClosure("mc_onlyPhotonCut_260pt", gjets+qcd, cut="260<pt&&pt<280")
    #qcdClosure("mc_onlyPhotonCut_280pt", gjets+qcd, cut="280<pt&&pt<300")
    #qcdClosure("mc_onlyPhotonCut_300pt", gjets+qcd, cut="300<pt&&pt<320")
    #qcdClosure("mc_onlyPhotonCut_320pt", gjets+qcd, cut="320<pt&&pt<340")
    #qcdClosure("mc_onlyPhotonCut_340pt", gjets+qcd, cut="340<pt&&pt<360")
    #qcdClosure("mc_onlyPhotonCut_360pt", gjets+qcd, cut="360<pt&&pt<380")
    #qcdClosure("mc_onlyPhotonCut_380pt", gjets+qcd, cut="380<pt&&pt<400")
    #qcdClosure("mc_onlyPhotonCut_400pt", gjets+qcd, cut="400<pt&&pt<420")
    #qcdClosure("mc_onlyPhotonCut_420pt", gjets+qcd, cut="420<pt&&pt<440")
    #qcdClosure("mc_onlyPhotonCut_440pt", gjets+qcd, cut="440<pt&&pt<460")
    #qcdClosure("mc_onlyPhotonCut_460pt", gjets+qcd, cut="460<pt&&pt<480")
    #qcdClosure("mc_onlyPhotonCut_480pt", gjets+qcd, cut="480<pt&&pt<500")
    #qcdClosure("mc_onlyPhotonCut_500pt", gjets+qcd, cut="500<pt&&pt<520")
    #qcdClosure("mc_onlyPhotonCut_520pt", gjets+qcd, cut="520<pt&&pt<540")
    #qcdClosure("mc_onlyPhotonCut_540pt", gjets+qcd, cut="540<pt&&pt<560")
    #qcdClosure("mc_onlyPhotonCut_560pt", gjets+qcd, cut="560<pt&&pt<580")
    #qcdClosure("mc_onlyPhotonCut_580pt", gjets+qcd, cut="580<pt&&pt<600")
    #qcdClosure("mc_onlyPhotonCut_600pt", gjets+qcd, cut="600<pt&&pt<620")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta1", gjets+qcd, cut="abs(eta)<.25")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta2", gjets+qcd, cut=".25<abs(eta)&&abs(eta)<0.5")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta3", gjets+qcd, cut=".5<abs(eta)&&abs(eta)<.75")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta4", gjets+qcd, cut=".75<abs(eta)&&abs(eta)<1")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta5", gjets+qcd, cut="1<abs(eta)&&abs(eta)<1.25")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta6", gjets+qcd, cut="1.25<abs(eta)&&abs(eta)<1.5")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta7", gjets+qcd, "tr_ee/simpleTree", cut="1.5<abs(eta)&&abs(eta)<1.75")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta8", gjets+qcd, "tr_ee/simpleTree", cut="1.75<abs(eta)&&abs(eta)<2")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta9", gjets+qcd, "tr_ee/simpleTree", cut="2<abs(eta)&&abs(eta)<2.25")
    #qcdClosureTwoCuts("mc_onlyPhotonCut_eta10", gjets+qcd, "tr_ee/simpleTree", cut="2.25<abs(eta)&&abs(eta)<2.5")

    #for scale in [1.]: finalDistributionEGscale("egScaled_{}_100pt".format(int(scale*100)), data, dataHt, cut="100<pt*{0}&&pt*{0}<120".format(scale), egScale=scale)
    #finalDistributionTwoCuts("onlyPhotonCut_100pt", data, dataHt, cut="100<pt&&pt<120")
    #finalDistributionTwoCuts("onlyPhotonCut_120pt", data, dataHt, cut="120<pt&&pt<140")
    #finalDistributionTwoCuts("onlyPhotonCut_140pt", data, dataHt, cut="140<pt&&pt<160")
    #finalDistributionTwoCuts("onlyPhotonCut_160pt", data, dataHt, cut="160<pt&&pt<180")
    #finalDistributionTwoCuts("onlyPhotonCut_180pt", data, dataHt, cut="180<pt&&pt<200")
    #finalDistributionTwoCuts("onlyPhotonCut_200pt", data, dataHt, cut="200<pt&&pt<220")
    #finalDistributionTwoCuts("onlyPhotonCut_220pt", data, dataHt, cut="220<pt&&pt<240")
    #finalDistributionTwoCuts("onlyPhotonCut_240pt", data, dataHt, cut="240<pt&&pt<260")
    #finalDistributionTwoCuts("onlyPhotonCut_260pt", data, dataHt, cut="260<pt&&pt<280")
    #finalDistributionTwoCuts("onlyPhotonCut_300pt", data, dataHt, cut="300<pt&&pt<320")
    #for scale in [1.01,1.02,1.03,1.04,1.05,1.06]: finalDistributionEGscale("egScaled_{}_300pt".format(int(scale*100)), data, dataHt, cut="300<pt*{0}&&pt*{0}<320".format(scale), egScale=scale)
    #for scale in [1.,1.01,1.02,1.03,1.04,1.05,1.06]: finalDistributionEGscale("egScaled_{}_300pt".format(int(scale*100)), data, dataHt, cut="300<pt*{0}&&pt*{0}<320".format(scale), egScale=scale)
    #for egScale in [0.99, 1., 1.01, 1.04]:
    #    correlation2d("met", (200,0,200), "sqrt(met**2+(pt*{0})**2+2*met*(pt*{0})*cos(dphi))".format(1-egScale), (200,0,200), drawMiddle=True)
    #correlation2d("pt", (200,100,1000), "emht", (200,700,2000), drawMiddle=True)

