#!/usr/bin/env python2

import re
import subprocess
import ROOT
import optparse
import os
from Datacard import Datacard


def infoFromOut(out):
    infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }
    for line in out.split("\n"):
        if line.startswith("NLL at global minimum of asimov:"): infos["rMinNLL"] = float(re.match(".*\(r = (.*)\).*",line).group(1))
        if line.startswith("Observed Limit: r < "): infos["obs"]   = float(line.split("<")[1])
        if line.startswith("Expected  2.5%: r < "): infos["exp2dn"] = float(line.split("<")[1])
        if line.startswith("Expected 16.0%: r < "): infos["exp1dn"] = float(line.split("<")[1])
        if line.startswith("Expected 50.0%: r < "): infos["exp"]    = float(line.split("<")[1])
        if line.startswith("Expected 84.0%: r < "): infos["exp1up"] = float(line.split("<")[1])
        if line.startswith("Expected 97.5%: r < "): infos["exp2up"] = float(line.split("<")[1])
    if infos["rMinNLL"] == 2:
        infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }
    return infos

def significanceFromOut(out):
    sig = -10
    for line in out.split("\n"):
        if line.startswith("Significance:"): sig = float(line.split(":")[1])
    return sig

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
    return infoFromOut(callCombine(name))

def callCombine(name):
    nameLimit = name + ".limit"
    if not os.path.isfile(nameLimit) or os.path.getmtime(name)>os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            out = subprocess.check_output(["combine", "-M", "Asymptotic", name], stderr=subprocess.STDOUT)
            f.write(out)
    with open(nameLimit) as f:
        out = f.read()
    return out

def callCombineSignificance(name):
    nameLimit = name + ".signif.limit"
    if not os.path.isfile(nameLimit) or os.path.getmtime(name)>os.path.getmtime(nameLimit):
        with open(nameLimit, "w+") as f:
            out = subprocess.check_output(["combine", "-M", "ProfileLikelihood", "--uncapped", "1", "--rMin", "-2", name], stderr=subprocess.STDOUT)
            f.write(out)
    with open(nameLimit) as f:
        out = f.read()
    return out

def getContour( gr2d ):
    gr2d.Draw()
    ROOT.gPad.Update()
    contoursN = [(c,c.GetN()) for c in gr2d.GetContourList(1.) ]
    contoursN = sorted( contoursN, key=lambda x: x[1] )
    if contoursN: return contoursN[-1][0]
    print "Could not find contour"

class MyDatacard(Datacard):
    def __init__(self, dc):
        if isinstance(dc, Datacard):
            self.bins = dc.bins
            self.obs = dc.obs
            self.processes = dc.processes
            self.signals = dc.signals
            self.isSignal = dc.isSignal
            self.keyline = dc.keyline
            self.exp = dc.exp
            self.systs = dc.systs
            self.shapeMap = dc.shapeMap
            self.hasShape = dc.hasShape
            self.flatParamNuisances = dc.flatParamNuisances
            self.rateParams = dc.rateParams
            self.rateParamsOrder = dc.rateParamsOrder
        elif dc == "photonHt":
            self.bins = []
            self.obs = {}
            self.processes = ['signal', 'gqcd', 'ele', 'zg', 'wg', 'ttg']
            self.signals = ['signal']
            self.isSignal = {'wg': False, 'signal': True, 'zg': False, 'gqcd': False, 'ele': False, 'ttg': False}
            self.keyline = []
            self.exp = {}
            self.systs = [(x, False, 'lnN', [], {}) for x in "lumi", "jec", "pdf", "gqcdSyst", "eleSyst", "wgSyst", "zgSyst", "ttgSyst"]
            self.shapeMap = {}
            self.hasShape = False
            self.flatParamNuisances = {}
            self.rateParams = {}
            self.rateParamsOrder = []

    def __str__(self):
        print "Class MyDatacard"
        print "bins               ", self.bins
        print "obs                ", self.obs
        print "processes          ", self.processes
        print "signals            ", self.signals
        print "isSignal           ", self.isSignal
        print "keyline            ", self.keyline
        print "exp                ", self.exp
        print "systs              ", self.systs
        #print "shapeMap           ", self.shapeMap
        #print "hasShape           ", self.hasShape
        #print "flatParamNuicances ", self.flatParamNuisances
        #print "rateParams         ", self.rateParams
        #print "rateParamsOrder    ", self.rateParamsOrder
        return ""

    def write(self, filename=""):
        maxInfoLen = max( [len(line[0])+7 for line in self.systs])
        out = ""
        out += "\nimax "+str(len(self.bins))
        out += "\njmax *"
        out += "\nkmax *"
        out += "\n\nbin         " + ("{:<15}"*len(self.bins)).format(*self.bins)
        out += "\nobservation " + ("{:<15}"*len(self.bins)).format(*[str(int(self.obs[x])) for x in self.bins])
        keylineInv0, keylineInv1, keylineInv2 = zip(*self.keyline)
        out += "\n\nbin".ljust(maxInfoLen) + ("{:<15}"*len(keylineInv0)).format(*keylineInv0)
        out += "\nprocess".ljust(maxInfoLen) + ("{:<15}"*len(keylineInv1)).format(*keylineInv1)
        counter = 1
        processNumbers = {}
        for a, b in self.isSignal.iteritems():
            if b: processNumbers[a] = 0
            else:
                processNumbers[a] = counter
                counter += 1
        out += "\nprocess".ljust(maxInfoLen) + ("{:<15}"*len(keylineInv1)).format(*[str(processNumbers[x]) for x in keylineInv1]) # name of processes has to be the same
        out += "\nrate".ljust(maxInfoLen) + ("{:<15}"*len(keylineInv0)).format(*[str(self.exp[keylineInv0[i]][keylineInv1[i]]) for i in range(len(keylineInv0))])
        for line in self.systs:
            out += "\n{} {} ".format(line[0],line[2]).ljust(maxInfoLen) + ("{:<15}"*len(keylineInv0)).format(*[str(line[4][keylineInv0[i]][keylineInv1[i]]) for i in range(len(keylineInv0))])
        if filename:
            with open(filename, "wb") as f:
                f.write(out)
                print "Writing to file:", filename
        else:
            print out

    def addBin(self, name, obs, gqcd, gqcdStat, gqcdSyst, ele, eleStat, eleSyst, zg, zgStat, zgSyst, wg, wgStat, wgSyst, ttg, ttgStat, ttgSyst, signal, signalStat, signalSyst):
        self.bins = list(set(self.bins+[name]))
        self.obs[name] = obs
        self.keyline.extend([(name, 'signal', True), (name, 'gqcd', False), (name, 'ele', False), (name, 'zg', False), (name, 'wg', False), (name, "ttg", False)])
        self.exp[name] = {"signal": signal, "gqcd": gqcd, "ele": ele,
                            "wg": wg, "zg": zg, "ttg": ttg}
        for line in self.systs:
            if   line[0] == "lumi": line[4][name] = {'signal': 1.062, 'gqcd': 1, 'ele': 1, 'wg': 1.062, 'zg': 1.062, 'ttg': 1.062}
            elif line[0] == "jec": line[4][name] = {'signal': 1.11, 'gqcd': 1, 'ele': 1, 'wg': 1.11, 'zg': 1.11, 'ttg': 1.11}
            elif line[0] == "pdf": line[4][name] = {'signal': 1.11, 'gqcd': 1, 'ele': 1, 'wg': 1.11, 'zg': 1.11, 'ttg': 1.11}
            elif line[0] == "gqcdSyst": line[4][name] = {'signal': 1, 'gqcd': gqcdSyst, 'ele': 1, 'wg': 1, 'zg': 1, 'ttg': 1}
            elif line[0] == "eleSyst": line[4][name] = {'signal': 1, 'gqcd': 1, 'ele': eleSyst, 'wg': 1, 'zg': 1, 'ttg': 1}
            elif line[0] == "wgSyst": line[4][name] = {'signal': 1, 'gqcd': 1, 'ele': 1, 'wg': wgSyst, 'zg': 1, 'ttg': 1}
            elif line[0] == "zgSyst": line[4][name] = {'signal': 1, 'gqcd': 1, 'ele': 1, 'wg': 1, 'zg': zgSyst, 'ttg': 1}
            elif line[0] == "ttgSyst": line[4][name] = {'signal': 1, 'gqcd': 1, 'ele': 1, 'wg': 1, 'zg': 1, 'ttg': ttgSyst}
            else: line[4][name] = {'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}

        self.systs.append(('qcdStat_'+name, False, 'lnN', [], dict(zip(self.bins, [{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]*(len(self.bins)-1)+[{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': gqcdStat, 'ele': 1, 'ttg': 1}]))))
        self.systs.append(('eleStat_'+name, False, 'lnN', [], dict(zip(self.bins, [{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]*(len(self.bins)-1)+[{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': eleStat, 'ttg': 1}]))))
        self.systs.append(('wgStat_'+name, False, 'lnN', [], dict(zip(self.bins, [{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]*(len(self.bins)-1)+[{'wg': wgStat, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]))))
        self.systs.append(('zgStat_'+name, False, 'lnN', [], dict(zip(self.bins, [{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]*(len(self.bins)-1)+[{'wg': 1, 'signal': 1, 'zg': zgStat, 'gqcd': 1, 'ele': 1, 'ttg': 1}]))))
        self.systs.append(('ttgStat_'+name, False, 'lnN', [], dict(zip(self.bins, [{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]*(len(self.bins)-1)+[{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': ttgStat}]))))
        self.systs.append(('signalStat_'+name, False, 'lnN', [], dict(zip(self.bins, [{'wg': 1, 'signal': 1, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]*(len(self.bins)-1)+[{'wg': 1, 'signal': signalStat, 'zg': 1, 'gqcd': 1, 'ele': 1, 'ttg': 1}]))))

    def newSignal(self, infos):
        if sorted(infos.keys()) != sorted(self.bins): print "Error, not all bins replaced:", infos.keys(), self.bins
        for b, v in infos.iteritems():
            self.exp[b]["signal"] = v[0]
            lines = [i for i, j in enumerate(self.systs) if j[0] == "signalStat_"+b]
            if not len(lines) != 1: "Error: Statistical uncertainty not found"
            self.systs[lines[0]][4][b]["signal"] = v[1]

    def limit(self):
        self.write("/tmp/tmpDataCard.txt")
        return infosFromDatacard("/tmp/tmpDataCard.txt")

    def limitFast(self):
        infos = { "obs":0, "exp":0, "exp1up":0, "exp1dn":0, "exp2up":0, "exp2dn":0 }

        for bin in self.bins:
            obs = self.obs[bin]
            bkg = sum([b for a,b in self.exp[bin].iteritems() if a is not "signal"])
            signal = self.exp[bin]["signal"]
            if not signal: continue
            err = ROOT.TMath.Sqrt(bkg)
            r = signal/abs(err - abs(obs - bkg))
            infos["obs"] = max(infos["obs"], r)
            r_exp = signal/err
            if r_exp > infos["exp"]:
                infos["exp"] = r_exp
                infos["exp1up"] = (signal+err)/err
                infos["exp2up"] = (signal+2*err)/err
                infos["exp1dn"] = (signal-err)/err
                infos["exp2dn"] = (signal-2*err)/err
        return infos

if False:
    inFileName = "limitCalculations/observation_v1.txt"
    import DatacardParser
    options, b = DatacardParser.addDatacardParserOptions(optparse.OptionParser())
    dc = MyDatacard(DatacardParser.parseCard(file(inFileName), options))

    dc.newSignal({
        "bin25": (12,1.1),
        "bin26": (1, 1.1),
        "bin27": (3, 1.2)
    })

    print dc
    dc.write()
    #print dc.limit()


