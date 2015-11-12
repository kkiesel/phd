#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import ROOT
import argparse

def packLists( inList ):
    inList.sort()

    outList = []
    tmpList = []
    for i in inList:
        if not tmpList: pass


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputFiles', default=["/user/kiesel/nTuples/V03_2/SinglePhoton_nTuple.root"], nargs="+" )
    args = parser.parse_args()

    ch = ROOT.TChain("TreeWriter/eventTree")
    for f in args.inputFiles:
        ch.AddFile( f )

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

