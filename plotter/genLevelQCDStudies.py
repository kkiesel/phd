#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
if sys.version_info[:2] == (2,6):
    print "Initialize correct python version first!"
    sys.exit()

import ConfigParser
import ROOT
import math
import argparse
import re
from random import randint
from sys import maxint

# private libs
import ratio
import style
import multiplot

import auxiliary as aux




def singleProcessEta():
    jFile = "../../CMSSW/genLevelChecks/CMSSW_7_4_14/src/GenWriter/GenWriter/histogrammer/QCD_merged_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_hists.root"
    gFile = "../../CMSSW/genLevelChecks/CMSSW_7_4_14/src/GenWriter/GenWriter/histogrammer/GJets_merged-TuneCUETP8M1_13TeV-madgraphMLM-pythia8_hists.root"
    c = ROOT.TCanvas()
    hg = aux.getFromFile( gFile, "etaG_gq->qG" )
    hj = aux.getFromFile( jFile, "etag1_gq->gq" )

    hj.SetLineColor( ROOT.kRed )

    for h in hg,hj:
        h.Rebin(4, "width")
        h.drawOption_="hist"
        h.Scale( 1./h.Integral() )
        h.SetTitle(";|#eta^{gen}|;normalized entries")
    hg.Scale( 2 ) # temporary fix: not abs eta taken

    m = multiplot.Multiplot()
    m.leg.SetY2(0.75)
    m.add( hg, "#gamma  gq#rightarrowq#gamma" )
    m.add( hj, "g  gq#rightarrowqg" )
    m.Draw()
    aux.save("genLevelQCDStudies_test" )



def simpleCompare(filename, name, settings, norm=True):
    file = "../../CMSSW/genLevelChecks/CMSSW_7_4_14/src/GenWriter/GenWriter/histogrammer/hists_nTuple_all100.root"
    file = "../../CMSSW/genLevelChecks/CMSSW_7_4_14/src/GenWriter/GenWriter/histogrammer/hists_slimmedTreeAN.root"
    colors = [ ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+3, ROOT.kMagenta ]

    c = ROOT.TCanvas()
    m = multiplot.Multiplot()
    for iS, (hname, leg) in enumerate( settings ):
        h = aux.getFromFile( filename, hname )
        if not h: continue
        if norm: h.Scale(1./h.Integral(0,-1))
#        if "pt" in name:
#            h = aux.rebin( h, range(100,200,20)+range(200,400,50)+[400, 600])
        h.drawOption_="hist e"
        h.SetLineColor(colors[iS])
        m.add( h, leg )
    m.Draw()
    aux.save("genLevelQCDStudies_"+name )
    c.SetLogy()
    aux.save("genLevelQCDStudies_"+name+"_log" )


#singleProcessEta()
#simpleCompare("emht", [("gjet_emht", "#gamma+Jet"),("qcd_emht", "QCD")])
#simpleCompare("nPart", [("gjet_nPart", "#gamma+Jet"),("qcd_nPart", "QCD")])
#simpleCompare("o2_eta_gjet", [("gjet_2_etaG", "#gamma"),("gjet_2_etaJ1","j")])
#simpleCompare("o2_eta_qcd", [("qcd_2_etaJ1", "j1"),("qcd_2_etaJ2","j2")])
#simpleCompare("o2_pt_gjet_qcd", [("gjet_2_ptG", "#gamma+Jet"),("qcd_2_ptJ1","QCD")])
#simpleCompare("o2_ht", [("gjet_2_emht", "EMH_{T}, #gamma+Jet"),("gjet_2_ht","H_{T}, #gamma+Jet"),("qcd_2_emht","(EM)H_{T}, QCD")])
#simpleCompare("o3_qcd_eta", [("qcd_3_etaJ1","j1"),("qcd_3_etaJ2","j2"),("qcd_3_etaJ3","j3")])
#simpleCompare("o3_qcd_pt", [("qcd_3_ptJ1","j1"),("qcd_3_ptJ2","j2"),("qcd_3_ptJ3","j3")])
#simpleCompare("o3_gjet_nPos", [("gjet_3_posG","")])
#simpleCompare("o3_gjet_pt1", [("gjet_3_ptJ1","j1"),("gjet_3_ptJ2","j2"),("gjet_3_ptG","#gamma")])
#simpleCompare("o3_gjet_pt_o", [("gjet_3_ptO1","obj1"),("gjet_3_ptO2","obj2"),("gjet_3_ptO3","obj3")])
#simpleCompare("o3_gjet_pt_lead", [("gjet_3_ptO1","obj1"),("gjet_3_ptO1G","obj1, #gamma"),("gjet_3_ptO1J","obj1, j")])
#simpleCompare("o3_gjet_pt_midd", [("gjet_3_ptO2","obj2"),("gjet_3_ptO2G","obj2, #gamma"),("gjet_3_ptO2J","obj2, j")])
#simpleCompare("o3_gjet_pt_trail", [("gjet_3_ptO3","obj3"),("gjet_3_ptO3G","obj3, #gamma"),("gjet_3_ptO3J","obj3, j")])
#simpleCompare("o3_gqcd_pt_lead", [("gjet_3_ptO1","obj1"),("gjet_3_ptO1G","obj1, #gamma"),("gjet_3_ptO1J","obj1, j"),("qcd_3_ptJ1","QCD obj1")])
#simpleCompare("o3_gqcd_pt_mid", [("gjet_3_ptO2","obj2"),("gjet_3_ptO2G","obj2, #gamma"),("gjet_3_ptO2J","obj2, j"),("qcd_3_ptJ2","QCD obj2")])
#simpleCompare("o3_gqcd_pt_trail", [("gjet_3_ptO3","obj3"),("gjet_3_ptO3G","obj3, #gamma"),("gjet_3_ptO3J","obj3, j"),("qcd_3_ptJ3","QCD obj3")])
#simpleCompare("o4_gjet_nPos", [("gjet_4_posG","4 objects"), ("gjet_5_posG", "5 objects")])
#simpleCompare("o4_gqcd_pt_lead", [("gjet_4_ptO1","obj1"),("gjet_4_ptO1G","obj1, #gamma"),("gjet_4_ptO1J","obj1, j"),("qcd_4_ptJ1","QCD obj1")])
#simpleCompare("o4_gqcd_pt_lead2", [("gjet_4_ptO2","obj2"),("gjet_4_ptO2G","obj2, #gamma"),("gjet_4_ptO2J","obj2, j"),("qcd_4_ptJ2","QCD obj2")])
#simpleCompare("o4_gqcd_pt_trail", [("gjet_4_ptO4","obj4"),("gjet_4_ptO4G","obj4, #gamma"),("gjet_4_ptO4J","obj4, j"),("qcd_4_ptJ4","QCD obj4")])
#simpleCompare("emht_nPart", [("gjet_emhtCut_nPart", "#gamma+Jet"),("qcd_emhtCut_nPart", "QCD"),("gjet_nPart", "#gamma+Jet (EMH_{T} inclusive)"),("qcd_nPart", "QCD (EMH_{T} inclusive)")])
#simpleCompare("emht_nPos", [("gjet_emhtCut_3_posG","3 objects"), ("gjet_emhtCut_4_posG", "4 objects"), ("gjet_emhtCut_5_posG", "5 objects")])

#simpleCompare("gjet_nPos", [("gjet_{}_posG".format(i),"{} objects".format(i)) for i in range(2,6)])

f1 = "../../CMSSW/genLevelChecks/CMSSW_7_4_14/src/GenWriter/GenWriter/histogrammer/hists_nTuple_all100.root"
f2 = "../../CMSSW/genLevelChecks/CMSSW_7_4_14/src/GenWriter/GenWriter/histogrammer/hists_slimmedTreeAN.root"

#for nObjects in range(2,6):
#    for object in range(1,nObjects+1):
#        simpleCompare(f2,"gqcdAN_n{}_o{}_pt".format(nObjects,object), [("gjet_{}_ptO{}".format(nObjects,object),"obj{}".format(object)),("gjet_{}_ptO{}G".format(nObjects,object),"obj{}, #gamma".format(object)),("gjet_{}_ptO{}J".format(nObjects,object),"obj{}, j".format(object)),("qcd_{}_ptJ{}".format(nObjects,object),"QCD obj{}".format(object))] )

def getProcessN(file):
    f = ROOT.TFile(file)
    hists = dict( [(key.GetName(),key.ReadObj()) for key in f.GetListOfKeys()] )
    nEntries = dict( [(name[:-5], h.Integral(0,-1)) for name, h in hists.iteritems() if name.endswith("_emht") and "TO" in name] )
    nEntriesG = dict( [(k[5:],v) for k,v in nEntries.iteritems() if k.startswith("gjet_") ] )
    nEntriesQ = dict( [(k[4:],v) for k,v in nEntries.iteritems() if k.startswith("qcd_") ] )
    import collections
    #nEntriesG = sorted(nEntriesG.items(), key=lambda x: len(x[0])*1e6+x[1], reverse=True)
    nEntriesG = sorted(nEntriesG.items(), key=lambda x: x[1], reverse=True)
    #nEntriesQ = sorted(nEntriesQ.items(), key=lambda x: len(x[0])*1e6+x[1], reverse=True)
    nEntriesQ = sorted(nEntriesQ.items(), key=lambda x: x[1], reverse=True)

    # number for 2/fb
    for name, n in nEntriesG: print "γJet", name, "   \t", 2000*n
    for name, n in nEntriesQ: print "QCD ", name, "   \t", 2000*n

#getProcessN(f1)
#getProcessN(f2)

#simpleCompare(f1, "gjet_nPos", [("gjet_{}_posG".format(i),"{} objects".format(i)) for i in range(2,6)])
#simpleCompare(f2, "gjetAN_nPos", [("gjet_{}_posG".format(i),"{} objects".format(i)) for i in range(2,6)])
#simpleCompare(f1, "gjet_nPart", [("gjet_nPart","#gamma+Jet"),("qcd_nPart","QCD")], norm=False)
#simpleCompare(f2, "gjetAN_nPart", [("gjet_nPart","#gamma+Jet"),("qcd_nPart","QCD")], norm=False)

