#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def ewkClosure(dataset, name, samplename, binning, binningName, dirDir="tr_genWZ11", preDir="tr_eControl"):

    f = 0.015 # fakeRate
    ef = 0.3 # relative uncertainy
    if "Weighted" in preDir: f=1.

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    hdir = aux.stdHist(dataset, dirDir+"/"+name, binning)
    aux.drawOpt(hdir, "data")
    m.add(hdir, "e#rightarrow#gamma")

    hpre = aux.stdHist(dataset, preDir+"/"+name, binning)
    hpre.SetLineColor(ROOT.kRed)
    hpre.drawOption_ = "hist"
    hpre.Scale(f)
    m.add(hpre, "")

    hsys = aux.getSysHisto(hpre, ef)
    aux.drawOpt(hsys, "sysUnc")
    m.add( hsys, "Prediction" )

    m.leg.SetHeader(preDir.replace("tr_eControl_",""))
    m.leg.SetX1(0.5)
    m.leg.SetY1(0.8)
    m.Draw()

    r = ratio.Ratio("e#rightarrow#gamma/Pred", hdir, hpre, hsys)
    r.draw(.5,1.5)

    l = aux.Label(sim=True,info=dataset.label)
    if binningName: binningName = "_"+binningName
    aux.save( "ewkClosure_{}_{}_vs_{}_{}{}".format(samplename,dirDir,preDir,name,binningName ) )

def ewkClosures(datasetName, dataset):
    names = aux.getObjectNames( dataset.files[0], "tr_eControl", [ROOT.TH1F] )
    names = ["met", "g_pt"]
    for name in names:
        for binningName, binning in aux.getBinningsFromName( name ).iteritems():
            ewkClosure(dataset, name, datasetName, binning, binningName)
            #ewkClosure(dataset, name, datasetName, binning, binningName, preDir="tr_eControl_etaWeighted")
            #ewkClosure(dataset, name, datasetName, binning, binningName, preDir="tr_eControl_ptWeighted")
            #ewkClosure(dataset, name, datasetName, binning, binningName, preDir="tr_eControl_ptFitWeighted")
            #ewkClosure(dataset, name, datasetName, binning, binningName, preDir="tr_eControl_vtxWeighted")

def ewkClosuresTree(datasetName, dataset, dirDir, preDir, variable, cut, binning):

    f = 0.015 # fakeRate
    ef = 0.3 # relative uncertainy
    if "Weighted" in preDir: f=1.

    can = ROOT.TCanvas()
    m = multiplot.Multiplot()
    hdir = aux.createHistoFromDatasetTree(dataset, variable, "1*({})".format(cut), binning, dirDir+"/simpleTree")
    aux.drawOpt(hdir, "data")
    m.add(hdir, "e#rightarrow#gamma")

    hpre = aux.createHistoFromDatasetTree(dataset, variable, "1*({})".format(cut), binning, preDir+"/simpleTree")
    hpre.SetLineColor(ROOT.kRed)
    hpre.drawOption_ = "hist"
    hpre.Scale(f)
    m.add(hpre, "")

    hsys = aux.getSysHisto(hpre, ef)
    aux.drawOpt(hsys, "sysUnc")
    m.add( hsys, "Prediction" )

    m.leg.SetHeader(preDir.replace("tr_eControl_",""))
    m.leg.SetX1(0.5)
    m.leg.SetY1(0.8)
    m.Draw()

    r = ratio.Ratio("e#rightarrow#gamma/Pred", hdir, hpre, hsys)
    r.draw(.5,1.5)

    l = aux.Label(sim=True,info=dataset.label)
    saveName = "ewkClosureTree_{}_{}_vs_{}_{}".format(datasetName, dirDir, preDir, variable, len(binning))
    aux.save(saveName)



if __name__ == "__main__":
    #ewkClosuresTree("tt", ttjets, "tr_genWZ11", "tr_eControl", "met", "1", range(0,200,10)+[200, 300, 400, 500, 600])
    ewkClosures("ewk", wjets+ttjets)
    ewkClosures("ewk_ht", wjets+ttjets_ht)
    ewkClosures("wjets", wjets)
    #ewkClosures("tt_ht", ttjets_ht)
    #ewkClosures("tt", ttjets)
    #ewkClosures("tt_nlo", ttjets_nlo)
    #ewkClosures("tt_ht0", ttjets0)
    #ewkClosures("tt_ht600", ttjets600)
    #ewkClosures("tt_ht800", ttjets800)
    #ewkClosures("tt_ht1200", ttjets1200)
    #ewkClosures("tt_ht2500", ttjets2500)
