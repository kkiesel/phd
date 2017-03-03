#!/usr/bin/env python2
import suppressor
with suppressor.suppress_stdout_stderr():
    import ROOT
import ROOT
import argparse
import os

def run(infile="", selector="HistogramProducer.cc", ext=False):
    # load libraries
    ROOT.gSystem.Load("pluginTreeWriterTreeWriterAuto.so")
    lib = "AutoDict_map_int_pair_int_int____cxx.so"
    if not os.path.isfile(lib):
        ROOT.gInterpreter.GenerateDictionary("map<int,pair<int,int> >", "map")
    ROOT.gSystem.Load("AutoDict_map_int_pair_int_int____cxx.so")

    if infile:
        ch = ROOT.TChain("TreeWriter/eventTree")
        ch.AddFile(infile)
        extName = infile.replace("_nTuple", "_ext_nTuple")
        if ext and os.path.isfile(extName):
            print "Add file", extName
            ch.AddFile(extName)
        ch.Process(selector+"+")
    elif ROOT.TSelector.GetSelector(selector+"++"):
        print "Compiled TSelector"
    else:
        raise Exception ('TSelector could not be compiled!')

def runExt(infile="", selector="HistogramProducer.cc"):
    # wrapper
    run(infile, selector, True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', default="", nargs="?")
    parser.add_argument('--ext', action='store_true')
    parser.add_argument('--signal', action='store_true')

    args = parser.parse_args()
    signalScans = ["SMS-T5Wg_nTuple.root", "SMS-T6Wg_nTuple.root", "SMS-T5Wg_mGo2150To2500_nTuple.root"]
    if os.path.basename(args.file) in signalScans or args.signal:
        run(args.file, "SignalScan.cc")
    else:
        run(args.file, ext=args.ext)
