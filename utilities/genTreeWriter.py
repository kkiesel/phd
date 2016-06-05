#!/usr/bin/env python2
# dataset = /ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/RunIIWinter15wmLHE-MCRUN2_71_V1-v1/LHE
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

import numpy as n

# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

lheeventproduct, label = Handle("LHEEventProduct"), "externalLHEProducer"

# open file (you can use 'edmFileUtil -d /store/whatever.root' to get the physical file name)
events = Events("root://xrootd-cms.infn.it//store/mc/RunIIWinter15wmLHE/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/LHE/MCRUN2_71_V1-v1/50000/9A7D08C2-E172-E511-A43B-0025905A60F4.root")

f = ROOT.TFile("test.root", "recreate")
t = ROOT.TTree("genTree", "Simple gen tree")

p3 = ROOT.TVector3()
p4 = ROOT.TVector3()
p5 = ROOT.TVector3()
p6 = ROOT.TVector3()
p7 = ROOT.TVector3()
t.Branch( 'p3', p3 )
t.Branch( 'p4', p4 )
t.Branch( 'p5', p5 )
t.Branch( 'p6', p5 )
t.Branch( 'p7', p7 )

for iev,event in enumerate(events):
    #if iev >= 10: break
    event.getByLabel(label, lheeventproduct)

    info = lheeventproduct.product().hepeup()

    p3.SetXYZ(0,0,0)
    p4.SetXYZ(0,0,0)
    p5.SetXYZ(0,0,0)
    p6.SetXYZ(0,0,0)
    p7.SetXYZ(0,0,0)

    for status, id, p in zip(info.ISTUP,info.IDUP,info.PUP):
        if status != 1: continue

        if abs(id) in [12,14,16]: # neutrino
            if p3.X()<1e-6: p3.SetXYZ(p[0],p[1],p[2])
            else: p4.SetXYZ(p[0],p[1],p[2]) # often not filled, why?
        elif id == 22:
            p5.SetXYZ(p[0],p[1],p[2])
        elif id == 21 or abs(id) in range(1,7):
            if p6.X()<1e-6: p6.SetXYZ(p[0],p[1],p[2]) # p5 and p6 same?
            else: p7.SetXYZ(p[0],p[1],p[2])
        else:
            print "unknown id", id

    t.Fill()

f.Write()
f.Close()
