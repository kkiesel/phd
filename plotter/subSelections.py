#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

colors = [rwth.blue, rwth.green, rwth.red, rwth.lila, rwth.magenta, rwth.petrol, rwth.cyan, rwth.orange, rwth.blue75] + range(50)

labels = {
    "tr_true_GUDSCB": "#gamma_{gen} from g/q",
    "tr_true_other": "#gamma_{gen} not from g/q",
    "tr_true_pi0": "#pi^{0}_{gen}",
    "tr_unmatched": "no #gamma_{gen}",
    "tr_genE": "e_{gen}",
    "tr_noGenE": "no e_{gen}",
}


def subSelection(sampleName, dataset, totalDir, subDirs, name, binning, binningName):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    hTot = aux.stdHist(dataset, totalDir+"/"+name, binning)

    for iSubDir, subDir in enumerate(subDirs):
        hSub = aux.stdHist(dataset, subDir+"/"+name, binning)
        if not hSub: continue
        hSub.Divide(hTot)
        hSub.SetLineColor(colors[iSubDir])
        hSub.SetYTitle("Fraction")
        m.addStack(hSub, labels[subDir] if subDir in labels else subDir)
    m.Draw()

    l = aux.Label(sim=dataset!=data, info=dataset.label)

    if binningName: binningName = "_"+binningName
    aux.save("subSelections_{}_{}{}".format(sampleName, name, binningName), log=False)


def subSelections(sampleName, dataset):
    names = aux.getObjectNames(dataset.files[0], "tr", [ROOT.TH1F])
    names = ["met"]
    for name in names:
        for binningName, binning in aux.getBinningsFromName(name).iteritems():
            if binningName != "3": continue
            #subSelection(sampleName+"_electrons", dataset, "tr", ["tr_genE", "tr_noGenE"], name, binning, binningName)
            subSelection(sampleName+"_electrons", dataset, "tr", ["tr_gen0","tr_gen-11", "tr_gen11", "tr_gen-13", "tr_gen13", "tr_gen-22", "tr_gen22"], name, binning, binningName)
            subSelection(sampleName+"_electronsWZ", dataset, "tr", ["tr_genWZ0", "tr_genWZ11", "tr_genWZ13", "tr_genWZ15","tr_genWZ12", "tr_genWZ14", "tr_genWZ16"], name, binning, binningName)
            subSelection(sampleName, dataset, "tr", ["tr_unmatched", "tr_true_pi0", "tr_true_other", "tr_true_GUDSCB"], name, binning, binningName)


if __name__ == "__main__":
    subSelections("wjets", wjets)
    subSelections("tt", ttjets)
    subSelections("gjets", gjets)
    subSelections("qcd", qcd)
    subSelections("ttg", ttg)
    subSelections("wg", wg)
    subSelections("znunu", znunu)
    subSelections("zg", zg)


