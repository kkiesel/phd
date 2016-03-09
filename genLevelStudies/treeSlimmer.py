#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import ROOT
import collections
import copy
import argparse

ROOT.gSystem.Load("templates/MyGenParticle.so")

def getOutName( inName, appendix="_nTuple_all100" ):
    return inName.split("/")[-1].replace("_nTuple",appendix)


parser = argparse.ArgumentParser()
parser.add_argument('inFile', default="../GJets_HT-merged-TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" )
args = parser.parse_args()



chain = ROOT.TChain("GenWriter/genParticleTree")
chain.AddFile( args.inFile )

f = ROOT.TFile( getOutName(args.inFile), "recreate")
newTree = chain.CloneTree(0)

for ievent, event in enumerate(chain):

    # only for gjets
    if event.g.Pt()>0.01 and (event.g.Pt() < 100 or abs(event.g.Eta())>1.4442): continue

    fill=True
    for j in event.jets:
        if j.p.Pt() < 100 or abs(j.p.Eta())>1.4442:
            fill=False
            break
    if fill: newTree.Fill()

newTree.AutoSave()

