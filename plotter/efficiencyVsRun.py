#!/usr/bin/env python2
# -*- coding: utf-8 -*-

def getRootMap( filename, dictname="triggerStudies/rawEff_vs_run" ):
    inFile = ROOT.TFile.Open( filename )
    ROOT.gInterpreter.GenerateDictionary("map<int,pair<int,int>>","map")
    myMap = ROOT.MakeNullPointer( ROOT.map("int,pair<int,int>") )
    inFile.GetObject(dictname,myMap)
    return myMap

def rootMap2pyDict( rootMap ):
    pyDict = {}
    for x,(y,z) in rootMap: pyDict[x] = (y,z)
    return pyDict

def getMergedDict( filenames ):
    m = rootMap2pyDict( getRootMap( filenames[0] ) )
    for f in filenames[1:]: m.update( rootMap2pyDict( getRootMap( f ) ) )
    return m

def drawEfficiencyVsRun(dataset, saveName=""):
    map = getMergedDict(dataset.files)

    od = collections.OrderedDict(sorted(map.items()))

    hrat = ROOT.TH1F("rat", "#gamma90,HT500 trigger efficiency measured with #gamma90 baseline trigger for H_{T}>600;;#varepsilon", len(map), 0, len(map) )
    hpas = ROOT.TH1F("hpas", "title", len(map), 0, len(map) )
    htot = ROOT.TH1F("htot", "title", len(map), 0, len(map) )

    for i, (run, (p,t)) in enumerate( od.iteritems() ):
        hpas.SetBinContent( i+1, p )
        htot.SetBinContent( i+1, t )
        hrat.GetXaxis().SetBinLabel( i+1, str(run) )
        hrat.SetBinContent( i+1, 1.*p/t )

    for h in hpas,htot: h.Sumw2(False)

    eff = ROOT.TEfficiency( hpas, htot )

    c = ROOT.TCanvas("","",2400,800)
    c.SetLeftMargin(0.05)
    hrat.SetMaximum(1)
    hrat.SetMinimum(0)
    hrat.Draw("axis")
    eff.Draw("p same Z")

    aux.save("triggerEfficiencyVsRun"+saveName)


def main():
    drawEfficiencyVsRun(dataHt)


if __name__ == "__main__":
    main()
