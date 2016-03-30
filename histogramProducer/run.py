#!/usr/bin/env python2
import ROOT
import argparse
import os

def run(infile="", selector="HistogramProducer.cc"):
    # load libraries
    ROOT.gSystem.Load("pluginTreeWriterTreeWriterAuto.so")
    lib = "AutoDict_map_int_pair_int_int____cxx.so"
    if not os.path.isfile(lib):
        ROOT.gInterpreter.GenerateDictionary("map<int,pair<int,int> >", "map")
    ROOT.gSystem.Load("AutoDict_map_int_pair_int_int____cxx.so")

    if infile:
        ch = ROOT.TChain("TreeWriter/eventTree")
        ch.AddFile(infile)
        ch.Process(selector+"+")
    else:
        ROOT.TSelector.GetSelector(selector+"++")
        print "Compiled TSelector"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', default="", nargs="?")

    args = parser.parse_args()

    run(args.file)
