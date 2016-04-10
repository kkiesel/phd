#!/usr/bin/env python2

import re
import argparse
import subprocess
import ROOT

def infoFromOut(out):
    infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }
    for line in out.split("\n"):
        if line.startswith("Observed Limit: r < "): infos["obs"]   = float(line.split("<")[1])
        if line.startswith("Expected  2.5%: r < "): infos["exp2dn"] = float(line.split("<")[1])
        if line.startswith("Expected 16.0%: r < "): infos["exp1dn"] = float(line.split("<")[1])
        if line.startswith("Expected 50.0%: r < "): infos["exp"]    = float(line.split("<")[1])
        if line.startswith("Expected 84.0%: r < "): infos["exp1up"] = float(line.split("<")[1])
        if line.startswith("Expected 97.5%: r < "): infos["exp2up"] = float(line.split("<")[1])
    return infos

def guessSignalPoint(name):
    m = re.match(".*_(\d+)_(\d+).*",name)
    if m:
        return int(m.group(1)), int(m.group(2))
    else:
        print "could not determine signal point for", name
        return 0, 0

def guessScanName( name ):
    short = "unknown"
    if "T5Wg" in name: short = "T5Wg"
    if "T5gg" in name: short = "T5gg"
    return short

def infosFromDatacard(name):
    return infoFromOut(subprocess.check_output( ["combine", "-M", "Asymptotic", name], stderr=subprocess.STDOUT ))


def getContour( gr2d ):
    gr2d.Draw()
    ROOT.gPad.Update()
    contoursN = [(c,c.GetN()) for c in gr2d.GetContourList(1.) ]
    contoursN = sorted( contoursN, key=lambda x: x[1] )
    if contoursN: return contoursN[-1][0]
    print "Could not find contour"


