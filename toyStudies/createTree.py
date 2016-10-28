#!/usr/bin/env python

import numpy as n
import ROOT
from math import *

rand = ROOT.TRandom2()


def getRandomMetResmeared():
    o1 = ROOT.TVector3()
    o1.SetPtEtaPhi(100*rand.Gaus(1,0.02)*rand.Gaus(1,0.2), 0, pi/4)
    o2 = ROOT.TVector3()
    o2.SetPtEtaPhi(200*rand.Gaus(1, .2), 0, pi/2)
    o3 = ROOT.TVector3()
    o3.SetPtEtaPhi(200*rand.Gaus(1, .2), 0, 3*pi/2)
    o4 = ROOT.TVector3()
    o4.SetPtEtaPhi(100*rand.Gaus(1, .2), 0, 5*pi/4)
    s = o1 + o2 + o3 + o4
    return s.Pt()

def getRandomMet(res=.2):
    o1 = ROOT.TVector3()
    o1.SetPtEtaPhi(100*rand.Gaus(1,res),0,pi/4)
    o2 = ROOT.TVector3()
    o2.SetPtEtaPhi(200*rand.Gaus(1, .2),0,pi/2)
    o3 = ROOT.TVector3()
    o3.SetPtEtaPhi(200*rand.Gaus(1, .2),0,3*pi/2)
    o4 = ROOT.TVector3()
    o4.SetPtEtaPhi(100*rand.Gaus(1, .2),0,5*pi/4)
    s = o1 + o2 + o3 + o4
    return s.Pt()

f = ROOT.TFile("toy_tree.root", "recreate")
t1 = ROOT.TTree("gam", "")
t2 = ROOT.TTree("jet", "")
t3 = ROOT.TTree("gamjet", "")

met1 = n.zeros(1, dtype=float)
met2 = n.zeros(1, dtype=float)
met3 = n.zeros(1, dtype=float)

t1.Branch('met', met1, 'met/D')
t2.Branch('met', met2, 'met/D')
t3.Branch('met', met3, 'met/D')

for i in xrange(2000000):
    met1[0] = getRandomMet(0.02)
    t1.Fill()
    met2[0] = getRandomMet()
    t2.Fill()
    met3[0] = getRandomMetResmeared()
    t3.Fill()

f.Write()
f.Close()


