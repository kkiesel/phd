#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from include import *

def drawSameHistogram(sampleNames, name, bkg=[], additional=[], binning=None, binningName="", scale=False):
    can = ROOT.TCanvas()
    m = multiplot.Multiplot()

    dataIntegral = 0
    for d in additional:
        if d == data or d == dataHt:
            h = d.getHist(name)
            dataIntegral = h.Integral(0,-1)
    bkgHist = None
    for d in bkg:
        h = d.getHist(name)
        if bkgHist: bkgHist.Add(h)
        else: bkgHist = h
    bkgIntegral = bkgHist.Integral(0,-1)

    for d in bkg[-1::-1]:
        h = d.getHist(name)
        if not h: continue
        if not h.Integral(): continue
        if scale: h.Scale(dataIntegral/bkgIntegral)
        if binning: h = aux.rebin( h, binning )

        aux.appendFlowBin( h )
        h.SetYTitle( aux.getYAxisTitle( h ) )
        m.addStack( h, d.label )

    dataHist = None
    for d in additional:
        h = d.getHist( name )
        if not h: continue
        if not h.Integral(): continue
        if binning: h = aux.rebin( h, binning )
        aux.appendFlowBin( h )

        if h.GetLineColor() == ROOT.kBlack: # data
            h.drawOption_ = "ep"
            h.SetMarkerStyle(20)
            # disable errors for data, so that ErrorOption is working
            # TODO: kPoisson also for rebinned and scaled histograms
            if not binning: h.Sumw2(False)
            h.SetBinErrorOption( ROOT.TH1.kPoisson )
            dataHist = h
        else:
            h.drawOption_ = "hist"
            h.SetLineWidth(3)

        m.add( h, d.label )

    m.sortStackByIntegral()
    if m.Draw():

        # ratio
        hsm = m.hists[0].GetStack().Last()
        if dataHist:
            r = ratio.Ratio( "Data/SM", dataHist, hsm )
            rMean = dataHist.Integral()/hsm.Integral()
            if rMean > 2 or rMean < 0.25: # for jcontrol with prescaled data
                r.draw(rMean/2,1.5*rMean)
            else:
                r.draw(.5,1.5)

        l = aux.Label(sim=data not in additional)

        if binningName: binningName = "_"+binningName
        name = name.replace("/","__")
        saveName = "sameHistograms_{}_{}{}".format(sampleNames, name, binningName )
        aux.save( saveName )

        if "emht" in name and dataHist:
            c = ROOT.TCanvas()
            myRatio = dataHist.Clone()
            myRatio.Divide(hsm)
            myRatio.SetMaximum(1.5)
            myRatio.SetMinimum(0.5)
            myRatio.Draw()
            aux.save(saveName+"_ratio", endings=[".root"], log=False)


def drawSameHistograms( sampleNames="test", stack=[], additional=[] ):
    file = stack[0].files[0] if stack else additional[0].files[0]
    names = aux.getObjectNames( file, "tr", [ROOT.TH1F] )
    dirs = [d for d in aux.getDirNames(file) if "tr_" in d and not "gen" in d and not "true" in d]

    if data in additional:
        if "genMatch" in names: names.remove("genMatch")
        if "genHt" in names: names.remove("genHt")

    names = ["met"]
    dirs = ["tr_jControl", "tr_eControl", "tr"]

    for name in names:
        for binningName, binning in aux.getBinningsFromName( name ).iteritems():
            #drawSameHistogram( sampleNames, "tr/"+name, stack, additional, binning, binningName )
            for directory in dirs:
                if "jControl" in directory:
                    thisAdditional = [dataHt] + [x for x in additional if x is not data]
                    drawSameHistogram( sampleNames, directory+"/"+name, stack, thisAdditional, binning, binningName )
                else:
                    drawSameHistogram( sampleNames, directory+"/"+name, stack, additional, binning, binningName )


if __name__ == "__main__":
    #transitions()
    #drawSameHistograms( "gqcd_data", [gjets, qcd], additional=[data])
    drawSameHistograms( "mc_data", [gjets, qcd, ttjets, ttg, wjets, wg, zg, znunu], additional=[data])
    drawSameHistograms( "mc", [gjets, qcd, ttjets, ttg, wjets, wg, zg, znunu], additional=[])

