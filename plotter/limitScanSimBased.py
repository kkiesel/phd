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
from string import Template

# private libs
import ratio
import style
import multiplot
from rwthColors import rwth
import limitTools

import auxiliary as aux
from datasets import *

def getMetYields( dataset ):
    h = dataset.getHist("tr/met")
    h = aux.rebin( h, binEdges, False )
    cont = [ h.GetBinContent(i) for i in range(1,h.GetNbinsX()+2) ]
    return cont


def getRgraphs( files ):

    binEdges = [ 150, 200, 250, 300, 370, 450 ]

    dataCardTemplate = Template("""
imax 6
jmax 1
kmax 2
bin bin0 bin1 bin2 bin3 bin4 bin5
observation ${observation}

bin bin0 bin1 bin2 bin3 bin4 bin5 bin0 bin1 bin2 bin3 bin4 bin5
process bg bg bg bg bg bg sig sig sig sig sig sig
process 1 1 1 1 1 1 0 0 0 0 0 0
rate ${bgRate} ${sigRate}
signalUnc lnN - - - - - - ${sigUnc}
bgUnc lnN ${bgUnc} - - - - - -
    """)

    dataCardTemplate = Template(dataCardTemplate.safe_substitute( {"sigUnc": "1.15 "*6, "bgUnc": "1.5 "*6 } ))

    bgYields = getMetYields( gjets+ qcd+ ttjets+ ttg+ wjets+wg_mg+zg_130+znunu )
    dataCardTemplate = Template( dataCardTemplate.safe_substitute( {"bgRate": " ".join([str(i) for i in bgYields]), "observation": " ".join([str(int(i)) for i in bgYields]) } ) )

    defaultGr = ROOT.TGraph2D(len(files))
    graphs = dict( (x,defaultGr.Clone(x)) for x in ["obs","exp","exp1up","exp1dn","exp2up","exp2dn"] )
    for ifile, file in enumerate(files):
        p = guessSignalPoint(file)
        pointName = file.split("/")[-1][:-11]
        sigSet = signal[pointName]

        sigYields = getMetYields( sigSet )

        with open("tmpDataCard.txt", "w") as f:
            f.write( dataCardTemplate.safe_substitute( {"sigRate": " ".join([str(i) for i in sigYields]) } ) )
        rInfo = limitTools.infosFromDatacard("tmpDataCard.txt")
        for name, gr in graphs.iteritems():
            graphs[name].SetPoint(ifile, p[0], p[1], rInfo[name] )

    return graphs



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()

    scanName = limitTools.guessScanName(args.files[0])

    graphs = getRgraphs( args.files )
    f = ROOT.TFile( "tmp/%s_graphs2d.root"%scanName,"recreate")
    for name, gr in graphs.iteritems():
        gr.Write(name)
    f.Close()


    f = ROOT.TFile("tmp/%s_graphs2d.root"%scanName)
    fout = ROOT.TFile( "tmp/%_graphs1d.root"%scanName, "recreate" )

    for name in "obs","exp","exp1up","exp1dn","exp2up","exp2dn":
        gr2d = f.Get(name)
        # get histogram
        if name == "obs":
            #h = ROOT.TH2F("","", 14, 950, 2050, 19, 50, 1950 )
            #h = gr2d.GetHistogram()
            #gr2d.SetHistogram(h)
            h = gr2d.Project("xy")
            fout.cd()
            h.Write("obs_hist")
        gr = getContour(gr2d)
        fout.cd()
        gr.Write(name)

    subprocess.call(["python", "smsPlotter/python/makeSMSplots.py", "smsPlotter/config/SUS15xxx/%s_SUS15xxx.cfg"%scanName, "plots/%s_limits_"%scanName])




