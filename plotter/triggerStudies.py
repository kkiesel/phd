#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def getRootMap(filename, dictname="triggerStudies/rawEff_vs_run"):
    inFile = ROOT.TFile.Open(filename)
    ROOT.gInterpreter.GenerateDictionary("map<int,pair<int,int>>", "map")
    myMap = ROOT.MakeNullPointer(ROOT.map("int,pair<int,int>"))
    inFile.GetObject(dictname, myMap)
    return myMap

def rootMap2pyDict(rootMap):
    pyDict = {}
    for x,(y,z) in rootMap: pyDict[x] = (y,z)
    return pyDict

def getMergedDict(filenames):
    m = {}
    for f in filenames: m.update(rootMap2pyDict(getRootMap(f)))
    return m

def drawEfficiencyVsRun(dataset, saveName=""):
    map = getMergedDict(dataset.files)
    od = collections.OrderedDict(sorted(map.items()))
    hrat = ROOT.TH1F("rat", "#gamma90,HT500 trigger efficiency measured with #gamma90 baseline trigger for H_{T}>600;;#varepsilon", len(map), 0, len(map))
    hpas = ROOT.TH1F("hpas", "title", len(map), 0, len(map))
    htot = ROOT.TH1F("htot", "title", len(map), 0, len(map))
    for i, (run, (p,t)) in enumerate(od.iteritems()):
        hpas.SetBinContent(i+1, p)
        htot.SetBinContent(i+1, t)
        hrat.GetXaxis().SetBinLabel(i+1, str(run))
        hrat.SetBinContent(i+1, 1.*p/t)
    for h in hpas,htot: h.Sumw2(False)
    eff = ROOT.TEfficiency(hpas, htot)
    c = ROOT.TCanvas("","",2400,800)
    c.SetLeftMargin(0.05)
    hrat.SetMaximum(1)
    hrat.SetMinimum(0)
    hrat.Draw("axis")
    eff.Draw("p same Z")
    aux.save("triggerEfficiencyVsRun"+saveName, log=False)

def efficiency(dataset, name, savename="", binning=None, binningName=""):
    c = ROOT.TCanvas()
    c.SetLogy(0)
    eff = dataset.getHist(name)
    if eff.UsesWeights(): eff.SetStatisticOption(ROOT.TEfficiency.kFNormal)

    h_pas = eff.GetPassedHistogram()
    h_tot = eff.GetTotalHistogram()

    if name.endswith("_preScaled"):
        eff2 = dataset.getHist(name.replace("_preScaled", ""))
        h_tot = eff2.GetTotalHistogram()

    if binning:
        h_pas = aux.rebin(h_pas, binning, False)
        h_tot = aux.rebin(h_tot, binning, False)

    if name.endswith("_preScaled"):
        ratio = h_pas.Clone(aux.randomName())
        ratio.Divide(h_tot)
        gr = ROOT.TGraphAsymmErrors(ratio)
        gr.SetTitle(";{};prescaled #varepsilon".format(eff.CreateGraph().GetHistogram().GetXaxis().GetTitle()))
        gr.SetLineColor(1)
        gr.Draw("ap")
    else:
        x = h_pas.Clone(aux.randomName())
        y = h_tot.Clone(aux.randomName())
        eff = ROOT.TEfficiency(h_pas, h_tot)
        h_pas, h_tot = x, y
        eff.Draw()
        ROOT.gPad.Update()
        gr = eff.GetPaintedGraph()

    gr.GetYaxis().SetRangeUser(0., 1.1)
    if name.endswith("eff_hlt_ht_ct") or name.endswith("eff_hlt_met_ct"):
        gr.GetYaxis().SetRangeUser(0., 0.1)


    if name.endswith("eff_hlt_pt") or name.endswith("eff_hlt_pt_endcap"):
        cutValue = 100
    elif name.endswith("eff_hlt_ht") \
            or name.endswith("eff_hlt_ht_ct") \
            or name.endswith("eff_hlt_ht_ct_preScaled"):
        cutValue = 700
    elif name.endswith("eff_hlt_ele_pt"):
        cutValue = 30
    else:
        cutValue = 0

    if cutValue or True:
        bin = h_pas.FindFixBin(cutValue)
        passed = int(h_pas.Integral(bin, -1))
        total = int(h_tot.Integral(bin, -1))
        if not total: return
        e = 1.*passed/total
        if passed<=total:
            conf = ROOT.TEfficiency().GetConfidenceLevel()
            e_up = ROOT.TEfficiency.ClopperPearson(total, passed, conf, True)
            e_dn = ROOT.TEfficiency.ClopperPearson(total, passed, conf, False)
        else:
            passed, epassed = aux.integralAndError(h_pas,bin,-1)
            total, etotal = aux.integralAndError(h_pas,bin,-1)
            ee = e * math.sqrt((epassed/passed)**2 + (etotal/total)**2)
            e_up = e + ee
            e_dn = e - ee
        eLabel = ROOT.TLatex(0.7, .15, "#varepsilon = {:.1f}^{{#plus{:.1f}}}_{{#minus{:.1f}}}%".format(100*e, 100*(e_up-e),100*(e-e_dn)))
        eLabel.SetNDC()
        eLabel.Draw()

        # graphical representation
        l = ROOT.TLine()
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kRed)
        xmin = gr.GetHistogram().GetXaxis().GetXmin()
        xmax = gr.GetHistogram().GetXaxis().GetXmax()
        l.DrawLine(max(cutValue,xmin), e, xmax, e)
        l.DrawLine(max(cutValue,xmin), e_up, xmax, e_up)
        l.DrawLine(max(cutValue,xmin), e_dn, xmax, e_dn)

        if cutValue > xmin:
            # cut line
            l.SetLineStyle(2)
            ymin = gr.GetYaxis().GetXmin()
            ymax = gr.GetYaxis().GetXmax()
            l.DrawLine(cutValue, ymin, cutValue, ymax)

    saveLumi = aux.intLumi
    if "B-F" in dataset.label:
        aux.intLumi = 20.101e3 # brilcalc lumi -b "STABLE BEAMS" --normtag=/afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -u /fb -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt --end 278808 # Up to run F
    l = aux.Label(sim="Data" not in dataset.label, info=dataset.label)
    aux.intLumi = saveLumi
    name = "_"+name.split("/")[-1]
    if binningName: binningName = "_"+binningName
    aux.save("efficiency_"+savename+name+binningName, log=False)

    if False:
        h_tot.SetLineColor(2)
        h_tot.Draw("hist")
        h_pas.Draw("same e*")
        aux.save("efficiency_"+savename+name+"_raw")

def efficiencies(dataset, savename=""):
    names = ["triggerStudies/"+x for x in aux.getObjectNames(dataset.files[0], "triggerStudies", [ROOT.TEfficiency]) ]
    for name in names:
        efficiency(dataset, name, savename)
        if name.endswith("eff_hlt_pt") or name.endswith("eff_hlt_pt_endcap"):
            efficiency(dataset, name, savename, range(0,80,8) + range(80,108,4) + range(108,300,12) + [300,400,500, 1000], "1")
        if name.endswith("eff_hlt_ht"):
            efficiency(dataset, name, savename, range(0,1001,40) + range(1000,1500,2000), "1")
        if name.endswith("eff_hlt_ht_ct") or name.endswith("eff_hlt_ht_ct_preScaled"):
            efficiency(dataset, name, savename, range(450,1001,10), "1")
        if name.endswith("eff_hlt_ht_ct2") or name.endswith("eff_hlt_ht_ct2_preScaled"):
            efficiency(dataset, name, savename, range(500,1001,50), "1")
        if "_met_" in name:
            efficiency(dataset, name, savename, range(0, 100, 5)+range(100,151,10), "1")
        if name.endswith("eff_hlt_nIso"):
            efficiency(dataset, name, savename, [0, .2, .4, .8, 1, 1.2, 1.4, 1.6, 2, 2.2, 2.4, 2.6, 2.8, 3, 4, 5, 6, 7, 8, 9, 10, 15], "1")
        if name.endswith("eff_hlt_r9"):
            efficiency(dataset, name, savename, [.4,.5,.6,.7,.8]+aux.frange(0.8, 1.1, 0.01), "1")
        if name.endswith("eff_hlt_sie"):
            efficiency(dataset, name, savename, [.0045,.007,.0075]+aux.frange(0.0075,0.011,5e-5), "1")

def main():
    #drawEfficiencyVsRun(dataHt)
    efficiencies(data, "singlePhoton")
    efficiencies(dataHt, "jetHt")
    #dataHtBtoF = Dataset("JetHT_Run2016B-PromptReco-v2", 0, ROOT.kBlack ) \
    #+ Dataset("JetHT_Run2016C-PromptReco-v2", 0, ROOT.kBlack ) \
    #+ Dataset("JetHT_Run2016D-PromptReco-v2", 0, ROOT.kBlack ) \
    #+ Dataset("JetHT_Run2016E-PromptReco-v2", 0, ROOT.kBlack ) \
    #+ Dataset("JetHT_Run2016F-PromptReco-v1", 0, ROOT.kBlack )
    #dataHtBtoF.label = "Data (JetHt B-F)"
    #efficiencies(dataHtBtoF, "jetHt_BtoF")

if __name__ == "__main__":
    main()
