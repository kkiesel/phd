#!/usr/bin/env python2
# -*- coding: utf-8 -*-

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

def packLists( inList ):
    inList.sort()

    outList = []
    tmpList = []
    for i in inList:
        if not tmpList: pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFile', default="/user/kiesel/nTuples/V03_2/SinglePhoton_nTuple.root" )
    args = parser.parse_args()

    ch = ROOT.TChain("TreeWriter/eventTree")
    ch.AddFile( args.inputFile )

    out = {}
    for e in ch:
        run = str( e.runNo )
        if run not in out:
            out[run] = []
        if e.lumNo not in out[run]:
            out[run].append( e.lumNo )

    runs = out.keys()
    runs.sort()
    print runs

    #print out

if __name__ == "__main__":
    main()

