#!/usr/bin/env python2
# import ROOT in batch mode

import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
from math import exp
from math import sqrt
import numpy as n
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

genP, genPlabel = Handle("vector<reco::GenParticle>"), "prunedGenParticles"
genJ, genJlabel = Handle("vector<reco::GenJet>"), "slimmedGenJets"

filenames = [
    ["/store/user/kiesel/13TeV/miniaod/SMS-T5gg/SMS-T5gg_mGluino-2000_mNeutralino-200/151109_163629/0000/T5gg_1.root"]
]

import subprocess

out = subprocess.check_output( "srmls -recursion_depth=4 srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/kiesel/13TeV/miniaod/SMS-T5gg", shell=True )

lines = out.split("\n")
lines = [ l.strip()[41:] for l in lines if l.endswith(".root") ]

lines = [ [l] for l in lines if l.endswith("_1.root") ] # only 1/4 of statistics

filenames = lines

for filename in filenames:
    prefix = "root://xrootd-cms.infn.it//"
    events = Events([ prefix+x for x in filename])

    outName = filename[0].split("/")[7] + ".root"
    f = ROOT.TFile(outName, "recreate")
    t = ROOT.TTree("signalTree", "")

    ht = n.zeros(1, dtype=float)
    pt = n.zeros(1, dtype=float)
    met = n.zeros(1, dtype=float)
    eta = n.zeros(1, dtype=float)
    nJet = n.zeros(1, dtype=float)

    t.Branch('ht', ht, 'ht/D')
    t.Branch('pt', pt, 'pt/D')
    t.Branch('met', met, 'met/D')
    t.Branch('eta', eta, 'eta/D')
    t.Branch('nJet', nJet, 'nJet/D')


    for iev,event in enumerate(events):
        #if iev > 100: break # do not use too much events
        event.getByLabel(genPlabel, genP )
        event.getByLabel(genJlabel, genJ )

        ht[0] = 0
        pt[0] = 0
        eta[0] = 0
        met[0] = 0
        nJet[0] = 0
        gravitinos = []

        for p in genP.product():
            if p.status() == 1 and p.pdgId() == 22 and pt[0] < p.pt():
                pt[0] = p.pt()
                eta[0] = p.eta()
            if p.status() == 1 and p.pdgId() == 1000039:
                gravitinos.append(p)
        for j in genJ.product():
            if j.pt() < 40 or abs(j.eta())>3: continue
            ht[0] += j.pt()
            nJet[0] += 1

        if len(gravitinos) != 2: print "there have to be exactly two gravitinos, not", len(gravitinos)

        met[0] = n.sqrt( (gravitinos[0].px()+gravitinos[1].px())**2 + (gravitinos[0].py()+gravitinos[1].py())**2 )

        t.Fill()



    f.Write()
    f.Close()
