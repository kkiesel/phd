#!/usr/bin/env python2
#import suppressor
#with suppressor.suppress_stdout_stderr():
#    import ROOT
import ROOT
import argparse
import os

def run(infile="", selector="FakeRateSelector.cc"):
    # load libraries
    ROOT.gSystem.Load("pluginTreeWriterTreeWriterAuto.so")

    if infile:
        ch = ROOT.TChain("TreeWriter/eventTree")
        ch.AddFile(infile)
        ch.Process(selector+"+")
    elif ROOT.TSelector.GetSelector(selector+"++"):
        print "Compiled TSelector"
    else:
        raise Exception ('TSelector could not be compiled!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', default="", nargs="?")

    args = parser.parse_args()

    run(args.file)
