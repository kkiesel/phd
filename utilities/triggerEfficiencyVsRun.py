#!/usr/bin/env python2
import argparse
import ROOT
import pickle

ROOT.gSystem.Load("../../tools/templates/TreeParticles.so")
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

def write():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile')
    args = parser.parse_args()

    ch = ROOT.TChain("TreeWriter/eventTree")
    ch.AddFile( args.inputFile )

    runs = {}
    susp = {}

    nEvents = ch.GetEntries()
    for ievent, event in enumerate( ch ):
        #if ievent>100: break
        if not ievent%100000: print "{:.1%}".format(1.*ievent/nEvents)

        if not event.HLT_Photon90_v: continue
        ht = 0
        for j in event.jets:
            if j.p.Pt() > 40 and abs(j.p.Eta())<3: ht += j.p.Pt()
        if ht < 600: continue

        if event.runNo not in runs: runs[event.runNo] = [0,0]
        runs[event.runNo][0] += 1
        if event.HLT_Photon90_CaloIdL_PFHT500_v: runs[event.runNo][1] += 1

        if event.runNo == 259637:
            if event.lumNo not in susp: susp[event.lumNo] = [0,0]
            susp[event.lumNo][0] += 1
            if event.HLT_Photon90_CaloIdL_PFHT500_v: susp[event.lumNo][1] += 1

    pickle.dump( runs, open( "runInfo.pkl", "wb" ) )
    pickle.dump( susp, open( "lumInfo.pkl", "wb" ) )

def drawEfficiency( map, savename ):
    sortedKeys = map.keys()
    sortedKeys.sort()

    hrat = ROOT.TH1F("rat", "#gamma90,HT500 trigger efficiency measured with #gamma90 baseline trigger for H_{T}>600;;#varepsilon", len(map), 0, len(map) )
    htot = ROOT.TH1F("htot", "title", len(map), 0, len(map) )
    hpas = ROOT.TH1F("hpas", "title", len(map), 0, len(map) )

    for ik,k in enumerate(sortedKeys):
        htot.SetBinContent( ik+1, map[k][1] )
        hpas.SetBinContent( ik+1, map[k][0] )

        hrat.GetXaxis().SetBinLabel( ik+1, str(k) )
        hrat.SetBinContent( ik+1, 1.*map[k][1]/map[k][0] )

    for h in htot,hpas:
        h.Sumw2()

    eff = ROOT.TEfficiency( htot,hpas )
    hrat.Draw("hist")
    eff.Draw("p same")

    ROOT.gPad.SaveAs(savename)



def read():
    runs = pickle.load( open( "runInfo.pkl", "rb" ) )
    drawEfficiency( runs, "eff_ht_vs_runs.pdf" )

    susp = pickle.load( open( "lumInfo.pkl", "rb" ) )
    drawEfficiency( susp, "eff_ht_vs_lumis.pdf" )


if __name__ == "__main__":
    #write()
    read()
