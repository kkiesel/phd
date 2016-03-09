#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import ROOT
import collections
import copy
import argparse

ROOT.gSystem.Load("templates/MyGenParticle.so")

def idToStr( pdgId ):
    if abs(pdgId) in range(1,7): return "q"
    if pdgId == 22: return "γ"
    if pdgId == 21: return "g"
    return str(pdgId)

def getProcessString( event ):
    init = ''.join(sorted([idToStr(event.idMother1),idToStr(event.idMother2)]))
    jetIds = ''
    for j in event.jets: jetIds+=idToStr(j.id)
    final = ''.join(sorted(jetIds))
    if event.g.Px():
        final += "G"
    return "{}->{}".format(init,final)

def getOutName( inName ):
    return inName.split("/")[-1].replace("_nTuple.root","_slimmedTreeAN.root")


parser = argparse.ArgumentParser()
parser.add_argument('inFile', default="../GJets_HT-merged-TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" )
args = parser.parse_args()



chain = ROOT.TChain("GenWriter/genParticleTree")
chain.AddFile( args.inFile )

f = ROOT.TFile( getOutName(args.inFile), "recreate")
newTree = chain.CloneTree(0)

if "GJets" in args.inFile:
    for event in chain:
        if event.g.Pt() < 100 or abs(event.g.Eta())>1.4442: continue
        emht = event.g.Pt()
        for j in event.jets:
            pt = j.p.Pt()
            if pt>40 and abs(j.p.Eta())<3: emht += pt
        if emht>700: newTree.Fill()
else:
    for event in chain:
        emht = 0
        maxJetPt = 0
        for j in event.jets:
            pt = j.p.Pt()
            if pt > maxJetPt: maxJetPt = pt
            if pt>40 and abs(j.p.Eta())<3: emht += pt
        if emht>700 and maxJetPt>100: newTree.Fill()


newTree.AutoSave()

